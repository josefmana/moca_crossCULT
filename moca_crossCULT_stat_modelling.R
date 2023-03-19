# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list packages to use
pkgs <- c( "tidyverse", "dplyr", # data wrangling
           "ggplot2", "patchwork", # plotting
           "brms", "tidybayes"
           )

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# prepare a folder for tables, figures, models, and sessions info
sapply( c("tabs","figs","mods","sess"), function(i) if( !dir.exists(i) ) dir.create(i) )


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

# add sum scores wherever missing
for ( i in c("ctr","exp") ) {
  for ( j in c("clock","naming","abstraction","frecall") ) d0[[i]][,j] <- rowSums( d0[[i]][ , grepl( paste0(j,"_"), colnames(d0[[i]]) ) ] )
}

# keep only the items list in struct
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

# prepare a data set with MoCA scores
d2 <- d1 %>% mutate(
  # re-code "subtraction" from raw correct answers to MoCA points and calculate the sum score
  subtraction = case_when( subtraction %in% c(5,4) ~ 3, subtraction %in% c(3,2) ~ 2, subtraction == 1 ~ 1, subtraction == 0 ~ 0 ),
  sum_score = rowSums( .[ , which(colnames(.) == "trail"):ncol(.) ] )
)

# prepare a data set with all memory scores
d3 <- d0$ctr[ , colnames(d0$ctr) %in% c("id",with( d0, names(ctr)[ grepl( "recall_|recog", names(ctr) ) ] ) ) ] %>%
  
  # add an index of the group and flag subjects to be included
  mutate( included = ifelse( id %in% d0$dem[ d0$dem$exp.id != "", "ctr.id"], 1, 0 ), .after = 1 ) %>%
  pivot_longer( cols = 3:ncol(.), names_to = c("type","item"), names_sep = "_", values_to = "score" ) %>%
  add_column( grp = "ctrl", .after = 1 ) %>%
  
  # bind the new experimental subjects below the old normative control ones
  rbind.data.frame(
    d0$exp[ , colnames(d0$exp) %in% c("id",with( d0, names(exp)[ grepl( "recall_|recog", names(exp) ) ] ) ) ] %>%
      
      # do the same as above with d0$ctr
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


# ---- session info ----

# write the sessionInfo() into a .txt file
capture.output( sessionInfo(), file = "sess/moca_crossCULT_stat_modelling.txt" )
