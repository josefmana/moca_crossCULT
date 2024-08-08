# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list packages to use
pkgs <- c( "tidyverse", "dplyr", # data wrangling
           "ggplot2", "patchwork", # plotting
           "rcompanion", # Vargha and Delaney’s A calculation
           "brms", "tidybayes" # Bayesian IRT model fitting and summaries
           )

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# prepare a folder for tables, figures, models, and sessions info
sapply( c("tabs","figs","mods","sess"), function(i) if( !dir.exists(i) ) dir.create(i) )

# prepare colorblind palette
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )


# ---- data set prep ----

# read out the structure of the MoCA test
struct <- read.csv( "data/moca_struct.csv", sep = "," )

# read the raw data
d0 <- list(
  dem = read.csv( "data/moca_crossCULT_data_demo.csv", sep = ";", dec = "," ),
  ctr = read.csv( "data/moca_crossCULT_data_ctr.csv", sep = ";" ) %>% mutate( id = as.character(id) ),
  exp = read.csv( "data/moca_crossCULT_data_exp.csv", sep = ";" )
)

# re-code "na" to NA, then change the variables to numeric
for ( i in names(d0)[2:3] ) {
  d0[[i]][ d0[[i]] == "na" ] <- NA # changing "na" to NA
  for ( j in names(d0[[i]])[-1] ) if( is.character( d0[[i]][,j] ) ) d0[[i]][,j] <- as.integer( d0[[i]][,j] ) # transforming to numeric
}

# re-code abstraction from 99 to 0 (99 means a concrete instead of an abstract answer)
for ( i in c("ctr","exp") ) {
  for ( j in paste0("abstraction_", 1:2) ) d0[[i]][,j] <- ifelse( d0[[i]][,j] == 99, 0, d0[[i]][,j] )
}

# change coding of the fluency score of the patient with shifted row
d0$ctr <- d0$ctr %>% mutate( fluency = ifelse( fluency > 1, 1, fluency ) )

# add item sum scores wherever missing
for ( i in c("ctr","exp") ) {
  for ( j in c("clock","naming","digits","repeating","abstraction","frecall") ) {
    # for each participant (row hence "rowSums") sum (hence "rowSums") correct answers for items in "j"
    d0[[i]][,j] <- rowSums( d0[[i]][ , grepl( paste0( j, "_" ), colnames( d0[[i]] ) ) ] )
  }
}

# for the item-level analysis data set keep only the items listed in struct
d1 <- d0$ctr[ , colnames(d0$ctr) %in% c("id",struct$item) ] %>%
  
  # add an index of the group and flag subjects to be included
  add_column( grp = "ctrl", .after = 1 ) %>%
  mutate( included = ifelse( id %in% d0$dem[ d0$dem$exp.id != "", "ctr.id"], 1, 0 ), .after = 2 ) %>%
  
  # bind the new experimental subjects below the old normative control ones
  rbind.data.frame(
    d0$exp[ , colnames(d0$exp) %in% c("id",struct$item) ] %>%
      add_column( grp = "exp", .after = 1 ) %>%
      mutate( included = ifelse( id %in% d0$dem[ d0$dem$exp.id != "", "exp.id"], 1, 0 ), .after = 2 )
  ) %>%
  
  # add columns for demographic variables
  add_column( age_y = NA, .after = 2 ) %>%
  add_column( edu_y = NA, .after = 3 ) %>%
  add_column( sex = NA, .after = 4 )

# add demographic variables via loops
for ( i in d1$id ) {
  
  if ( i %in% d0$dem$ctr.id ) d1[ d1$id == i , c("age_y","edu_y","sex") ] <- d0$dem[ d0$dem$ctr.id == i, paste0( "ctr.",c("age_y","edu_y","sex") ) ]
  if ( i %in% d0$dem$exp.id ) d1[ d1$id == i , c("age_y","edu_y","sex") ] <- d0$dem[ d0$dem$exp.id == i, paste0( "exp.",c("age_y","edu_y","sex") ) ]
  if ( i %in% d0$dem$excl.id ) d1[ d1$id == i , c("age_y","edu_y","sex") ] <- d0$dem[ d0$dem$excl.id == i, paste0( "excl.",c("age_y","edu_y","sex") ) ]

}

# prepare a data set with total MoCA scores for surface-level analysis
d2 <- d1 %>%
  # re-code "subtraction" from raw correct answers to MoCA points and calculate the sum score
  mutate( subtraction = case_when( subtraction %in% c(5,4) ~ 3, subtraction %in% c(3,2) ~ 2, subtraction == 1 ~ 1, subtraction == 0 ~ 0 ) ) %>%
  mutate( sum_score = rowSums( .[ , which(colnames(.) == "trail"):ncol(.) ] ) )
  

# switch the item-level data set (d1) to long format for IRT
# (note that this step comes after creating d2 because d2 was built upon a wide format d1)
d1 <- d1 %>% pivot_longer( cols = 7:ncol(.), names_to = "item", values_to = "score" ) %>%
  # make the item name to an ordered factor for straighforward plotting
  mutate( item = factor( item, levels = struct$item, ordered = T) )

# prepare a data set with all memory scores for a task-level analysis
d3 <- d0$ctr[ , colnames(d0$ctr) %in% c("id",with( d0, names(ctr)[ grepl( "recall_|recog", names(ctr) ) ] ) ) ] %>%
  
  # add an index of the group and flag subjects to be included
  mutate( included = ifelse( id %in% d0$dem[ d0$dem$exp.id != "", "ctr.id"], 1, 0 ), .after = 1 ) %>%
  pivot_longer( cols = 3:ncol(.), names_to = c("type","item"), names_sep = "_", values_to = "score" ) %>%
  add_column( grp = "ctrl", .after = 1 ) %>%
  
  # bind the new experimental subjects below the old normative control ones
  rbind.data.frame(
    d0$exp[ , colnames(d0$exp) %in% c("id",with( d0, names(exp)[ grepl( "recall_|recog", names(exp) ) ] ) ) ] %>%
      
      # do the same as above but with d0$exp instead
      mutate( included = ifelse( id %in% d0$dem[ d0$dem$exp.id != "", "exp.id"], 1, 0 ), .after = 1 ) %>%
      pivot_longer( cols = 3:ncol(.), names_to = c("type","item"), names_sep = "_", values_to = "score" ) %>%
      add_column( grp = "exp", .after = 1 )
    
  ) %>%
  
  # add order of each item in the list
  mutate( order = case_when( item == "face" ~ "i1",
                             item == "velvet" ~ "i2",
                             item == "church" ~ "i3",
                             item %in% c("daisy","rye") ~ "i4",
                             item %in% c("red","salt") ~ "i5"
                             ), .before = "score" )


# ---- descriptive stats ----

# loop through control and experimental groups
t1 <- lapply( unique(d2$grp), function(i)
  
  # loop through the continuous variables
  sapply( names(d2)[c(3,4,18,19)], function(j)
    
    # extract a vector of values of each outcome then compute descriptive stats
    na.omit( d2[ with(d2, grp == i & included == 1 ) , j ] ) %>% ( function(x = .) {
      c( length(x), # number of observations
         round( mean(x), 2 ) %>% sprintf( "%.2f", . ), # mean
         round( sd(x), 2 ) %>% sprintf( "%.2f", . ), # standard deviation
         round( median(x), 2 ) %>% sprintf( "%.2f", . ) ) # median
      }
    )
  ) %>%
    
    # some book-keeping to prepare for the finish line
    t() %>% `colnames<-`( c("N", "M", "SD", "Md") ) %>%
    as.data.frame() %>%
    rownames_to_column( var = "Variable" ) %>%
    mutate( grp = i, .before = 1 )
  
  # pull SA and non-SA tables together
) %>% do.call( rbind.data.frame, . ) %>%
  pivot_wider( values_from = c("N","M","SD","Md"),
               names_from = "grp",
               names_glue = "{grp}_{.value}"
  ) %>%
  
  # reorder columns such that they make a sense
  relocate( 1,2,4,6,8,3,5,7,9 )

# add sex numbers
t1[5, ] <- c( "sex_m",
              
              # non-sex way to extract the numbers
              table( d2[ with( d2, included == 1 & grp == "ctrl" ), "sex" ] ) %>% as.data.frame() %>% select(Freq) %>%
                mutate( N = paste0( Freq, " (", ( 100 * Freq / sum(Freq) ) %>% round(2) %>% sprintf( "%.2f", . ), "%)" ) ) %>%
                select(N) %>% slice(1),
              # blank cells
              rep("-",3),
              
              # the same for experimental group
              table( d2[ with( d2, included == 1 & grp == "exp" ), "sex" ] ) %>% as.data.frame() %>% select(Freq) %>%
                mutate( N = paste0( Freq, " (", ( 100 * Freq / sum(Freq) ) %>% round(2) %>% sprintf( "%.2f", . ), "%)" ) ) %>%
                select(N) %>% slice(1),
              # blank cells
              rep("-",3)
              )

# save as .csv
write.table( t1[ c(1,2,5,3,4), ], "tabs/t01_sample_description.csv", sep = ",", quote = F, row.names = F )


# ---- surface-level analysis ----

# plot all raw frequencies into a grid of histograms
d2 %>%
  filter( included == 1 ) %>% # keep included subjects only
  pivot_longer( cols = 7:(ncol(.)-1), names_to = "item", values_to = "score" ) %>%
  mutate( label = sapply( item, function(i) struct[ struct$item == i , "label"] ) %>% factor( levels = struct$label, ordered = T ),
          Group = case_when( grp == "ctrl" ~ "MoCA", grp == "exp" ~ "MoCA-WLE" ) %>% factor( levels = c("MoCA","MoCA-WLE"), ordered = T )
          ) %>%
  
  # plotting proper
  ggplot( aes(x = score, fill = Group ) ) + # colors by group
  geom_histogram( binwidth = .3, position = position_dodge(.3) ) + # make the histograms
  scale_x_continuous( breaks = seq(0,6,1), labels = seq(0,6,1) ) + # only integers on abscissa
  scale_fill_manual( values = cbPal[c(5,6)] ) + # colorblind friendly color palette
  labs( x = "Score", y = "Count" ) + # axes names
  facet_wrap( ~ label, nrow = 4, ncol = 3, scales = "free" ) + # make a grid
  theme_bw( base_size = 16 ) + theme( legend.position = "bottom" ) # final touches

# save asi Fig 1 (for now)
ggsave( "figs/fig1_item_scores_frequencies.jpg", dpi = 300, width = 12, height = 13.3 )

# next do NHST analysis on the item scores
t2 <- data.frame( `0` = rep(NA, length(struct$item) ), `1` = NA, `2` = NA, `3` = NA, `4` = NA, `5` = NA, `6` = NA, # frequencies
                  W = NA, p = NA, VDA  = NA, # test statistic, p-value and effect size,
                  row.names = struct$item )

# fill-in frequencies in a format N_ctr/N_exp wherever appropriate or "-" wherever not
for ( i in rownames(t2) ) {
  
  # fill-in frequencies in a format N_ctr/N_exp wherever appropriate or "-" wherever not
  for ( j in 1:7 ) t2[i,j] <- tryCatch(
    table( d2[ d2$included == 1, i ], d2[ d2$included == 1, "grp" ] )[as.character(j-1),] %>% paste( collapse = "/" ),
    error = function(e) "-"
  )
  
  # add stats
  t2[ i , "W" ] <- wilcox.test( as.formula( paste0( i, " ~ grp" ) ), data = d2[ d2$included == 1, ], paired = F )$statistic %>% round(1) %>% sprintf( "%.1f", . )
  t2[ i , "p" ] <- wilcox.test( as.formula( paste0( i, " ~ grp" ) ), data = d2[ d2$included == 1, ], paired = F )$p.value %>% round(3) %>% sprintf( "%.3f", . )
  
  # before the next step set seed for exact replocation of bootstapped CIs
  set.seed(87542)
  
  # fill-in the effect siue (Vargha and Delaney’s A)
  t2[ i, "VDA" ] <- paste0(
    vda( as.formula( paste0( i, " ~ grp" ) ), data = d2[ d2$included == 1, ], ci = F, conf = .95 ) %>% round(2) %>% sprintf( "%.2f", . ), " [",
    vda( as.formula( paste0( i, " ~ grp" ) ), data = d2[ d2$included == 1, ], ci = T, conf = .95 )[,2:3] %>% round(2) %>% sprintf( "%.2f", . ) %>% paste( collapse = ", " ), "]"
  )

}

# save the table to .csv
write.table( t2 %>% add_column( Item = struct$label, .before = 1 ), "tabs/t1_nhst_results.csv", sep = ";", quote = F, row.names = F )

# for in-text reporting calculate the same for the MoCA total score
set.seed(87542)
print(
  paste0(
    "MoCA total score: W = ", wilcox.test( sum_score ~ grp, data = d2[ d2$included == 1, ], paired = F )$statistic %>% round(1) %>% sprintf( "%.1f", . ),
    ", p = ", wilcox.test( sum_score ~ grp, data = d2[ d2$included == 1, ], paired = F )$p.value %>% round(3) %>% sprintf( "%.3f", . ),
    ", VDA = ", vda( sum_score ~ grp, data = d2[ d2$included == 1, ], ci = F, conf = .95 ) %>% round(2) %>% sprintf( "%.2f", . ),
    " (95% CI [", vda( sum_score ~ grp, data = d2[ d2$included == 1, ], ci = T, conf = .95 )[,2:3] %>% round(2) %>% sprintf( "%.2f", . ) %>% paste( collapse = ", " ), "])"
  )
)


t# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/moca_crossCULT_stat_modelling.txt" )
