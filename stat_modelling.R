# This is a script that analyses MoCA data via a series of M-U tests and multilevel Bayesian binomial model
#
# The analysis works for all but the 'DESCRIPTIVES' section because we do not share continuous demography data
# for privacy reasons to achiev k-anonymity (on k ≥ 5 level)


rm( list = ls() ) # clear environment
options( mc.cores = parallel::detectCores() ) # use all parallel CPU cores

# load packages
library(here)
library(tidyverse)
library(patchwork)
library(rcompanion)
library(lavaan)
library(psych)
library(brms)
library(priorsense)
library(tidybayes)

# create folders for models, figures, tables and sessions to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present
sapply( c("figures", "tables"), function(i) if( !dir.exists(i) ) dir.create(i) )

# prepare colorblind palette
cbPal <- c( "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7" )


# DATA SET ----

# read out the structure of the MoCA test
struct <- read.csv( here("data","moca_struct.csv"), sep = "," )

# read the raw data
d0 <- list(
  #dem = read.csv( here("data","_data_demo.csv"), sep = ";", dec = "," ), # hashtag this line as I do not share full demography variables for privacy reasons
  dem = read.csv( here("data","data_demo.csv"), sep = ";", dec = "," ),
  ctr = read.csv( here("data","data_ctr.csv"), sep = ";" ) %>% mutate( id = as.character(id) ),
  exp = read.csv( here("data","data_exp.csv"), sep = ";" )
)

# re-code "na" to NA, then change the variables to numeric
for ( i in c("ctr","exp") ) {
  
  d0[[i]][ d0[[i]] == "na" ] <- NA # changing "na" to NA
  for ( j in names(d0[[i]])[-1] ) {
    if( is.character( d0[[i]][ ,j] ) ) d0[[i]][ ,j] <- as.integer( d0[[i]][ ,j] ) # transforming to numeric
  }
}

# re-code abstraction from 99 to 0 (99 means a concrete instead of an abstract answer)
for ( i in c("ctr","exp") ) for ( j in paste0("abstraction_", 1:2) ) d0[[i]][,j] <- ifelse( d0[[i]][,j] == 99, 0, d0[[i]][,j] )

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
  mutate(
    subtraction = case_when(
      subtraction %in% c(5,4) ~ 3,
      subtraction %in% c(3,2) ~ 2,
      subtraction == 1 ~ 1,
      subtraction == 0 ~ 0
    )
  ) %>%
  mutate( sum_score = rowSums( .[ , which(colnames(.) == "trail"):ncol(.) ] ) )

# switch the item-level data set (d1) to long format for IRT
# (note that this step comes after creating d2 because d2 was built upon a wide format d1)
d1 <- d1 %>%
  
  pivot_longer( cols = 7:ncol(.), names_to = "item", values_to = "score" ) %>%
  mutate(
    item = factor(
      unlist( sapply( 1:nrow(.), function(i) struct$label[struct$item == item[i]] ), use.names = F),
      levels = struct$label,
      ordered = T
    ),
    max = unlist(
      sapply( 1:nrow(.), function(i) struct$max_score[struct$label == item[i]] ),
      use.names = F
    )
  )


# DESCRIPTIVES ----

# loop through control and experimental groups
t1 <- lapply(
  
  unique(d2$grp),
  function(i)
    
    # loop through the continuous variables
    sapply(
      
      names(d2)[c(3,4,18,19)],
      function(j)
        
        # extract a vector of values of each outcome then compute descriptive stats
        na.omit( d2[ with(d2, grp == i & included == 1 ) , j ] ) %>%
        (function(x = .) {
          c( length(x), # number of observations
             round( mean(x), 2 ) %>% sprintf( "%.2f", . ), # mean
             round( sd(x), 2 ) %>% sprintf( "%.2f", . ), # standard deviation
             round( median(x), 2 ) %>% sprintf( "%.2f", . ) # median
             )
        }
      
    )
  ) %>%
    
    # some book-keeping to prepare for the finish line
    t() %>%
    `colnames<-`( c("N", "M", "SD", "Md") ) %>%
    as.data.frame() %>%
    rownames_to_column( var = "Variable" ) %>%
    mutate( grp = i, .before = 1 )
  
# pull SA and non-SA tables together
) %>%
  
  do.call( rbind.data.frame, . ) %>%
  pivot_wider(
    values_from = c("N","M","SD","Md"),
    names_from = "grp",
    names_glue = "{grp}_{.value}"
  ) %>%
  
  # reorder columns such that they make a sense
  relocate(1,2,4,6,8,3,5,7,9)

# add comparison of means/chi-square results
t1 <- t1 %>%
  
  left_join(
    
    sapply(
      t1$Variable,
      function(y)
        unlist(
          t.test( as.formula( paste0(y," ~ grp") ), data = subset(d2, included == 1) )[c("statistic","parameter","p.value")],
          use.names = F
        ) %>%
        `names<-`( c("Stat","df","p_value") ) %>%
        c(Variable = y, Test = "t-test", . )
    ) %>%
      t() %>%
      as.data.frame() %>%
      mutate(
        Stat = sprintf( "%.3f", round( as.numeric(Stat), 3 ) ),
        df = sprintf( "%.2f", round( as.numeric(df), 2 ) ),
        p_value = sprintf( "%.3f", round( as.numeric(p_value), 3 ) )
      ),
    
    by = "Variable"
    
  )

# add sex numbers
t1[5, ] <- c(
  
  # Variable
  "sex_m",
  
  # non-sex way to extract the numbers
  table( ( d2[ with( d2, included == 1 & grp == "ctrl" ), ] %>% mutate(sex = gsub( " ", "", gsub("ž","z",sex) ) ) )$sex ) %>%
    as.data.frame() %>%
    select(Freq) %>%
    mutate( N = paste0( Freq, " (", ( 100 * Freq / sum(Freq) ) %>% round(2) %>% sprintf( "%.2f", . ), "%)" ) ) %>%
    select(N) %>%
    slice(1),
  
  # blank cells
  rep("-",3),
  
  # the same for experimental group
  table( ( d2[ with( d2, included == 1 & grp == "ctrl" ), ] %>% mutate(sex = gsub( " ", "", gsub("ž","z",sex) ) ) )$sex ) %>%
    as.data.frame() %>%
    select(Freq) %>%
    mutate( N = paste0(Freq, " (", ( 100 * Freq / sum(Freq) ) %>% round(2) %>% sprintf( "%.2f", . ), "%)") ) %>%
    select(N) %>%
    slice(1),
  
  # blank cells
  rep("-",3),
  
  # statistics
  # stupidly long code because I could not care less now about making it tidy and the data are formatted pityfully
  "Chi-square",
  sprintf( "%.3f", round(chisq.test( table( subset(d2, included == 1)[ , c("sex","grp") ] %>% mutate(sex = gsub( " ", "", gsub("ž","z",sex) ) ) ) )$statistic, 3) ),
  sprintf( "%.0f", round(chisq.test( table( subset(d2, included == 1)[ , c("sex","grp") ] %>% mutate(sex = gsub( " ", "", gsub("ž","z",sex) ) ) ) )$parameter, 0) ),
  sprintf( "%.3f", round(chisq.test( table( subset(d2, included == 1)[ , c("sex","grp") ] %>% mutate(sex = gsub( " ", "", gsub("ž","z",sex) ) ) ) )$p.value, 3) )
  
)


# save as .csv
write.table(
  x = t1[ c(1,2,5,3,4), ],
  file = here("tables","sample_description.csv"),
  sep = ",",
  quote = F,
  row.names = F
)


# SURFACE-LEVEL ANALYSIS ----

# do NHST analysis on the item scores
t2 <- data.frame(
  
  # frequencies
  `0` = rep(NA, length(struct$item) ),
  `1` = NA,
  `2` = NA,
  `3` = NA,
  `4` = NA,
  `5` = NA,
  `6` = NA,
  
  # test statistic, p-value and effect size
  W = NA,
  p = NA,
  VDA  = NA,
  
  # rownames
  row.names = struct$item
  
)

# fill-in frequencies in a format N_ctr/N_exp wherever appropriate or "-" wherever not
for ( i in rownames(t2) ) {
  
  # fill-in frequencies in a format N_ctr/N_exp wherever appropriate or "-" wherever not
  for ( j in 1:7 ) t2[i,j] <- tryCatch(

    table( d2[ d2$included == 1, i ], d2[ d2$included == 1, "grp" ] )[as.character(j-1),] %>%
      paste( collapse = "/" ),
    error = function(e) "-"

  )
  
  # add stats
  t2[ i , "W" ] <- wilcox.test( as.formula( paste0( i, " ~ grp" ) ), data = d2[ d2$included == 1, ], paired = F )$statistic %>% round(1) %>% sprintf( "%.1f", . )
  t2[ i , "p" ] <- wilcox.test( as.formula( paste0( i, " ~ grp" ) ), data = d2[ d2$included == 1, ], paired = F )$p.value %>% round(3) %>% sprintf( "%.3f", . )
  
  # before the next step set seed for exact replocation of bootstapped CIs
  set.seed(87542)
  
  # fill-in the effect size (Vargha and Delaney’s A)
  t2[ i, "VDA" ] <- paste0(

    vda( as.formula( paste0( i, " ~ grp" ) ), data = d2[ d2$included == 1, ], ci = F, conf = .95 ) %>% round(2) %>% sprintf( "%.2f", . ), " [",
    vda( as.formula( paste0( i, " ~ grp" ) ), data = d2[ d2$included == 1, ], ci = T, conf = .95 )[ ,2:3] %>% round(2) %>% sprintf( "%.2f", . ) %>% paste( collapse = ", " ), "]"

  )

}

# save the table to .csv
write.table(
  x = add_column(t2, Domain = struct$label, .before = 1),
  file = here("tables","nhst_results.csv"),
  sep = ";",
  quote = F,
  row.names = F
)

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


# keep only included patients in data sets used in next sections
d2 <- subset(d2, included == 1)
d1 <- subset(d1, id %in% d2$id)


# INTERNAL CONSISTENCY ----

# confirmatory factor analysis
form0 <- 'moca =~ trail + cube + clock + naming + digits + tapping + subtraction + repeating + fluency + abstraction + frecall + orientation'
summary(cfa(form0, data = d2), standardized = T, fit.measures = T) # not very good

# internal consistency measures
#alpha( d2[ , struct$item] )
omega(d2[ , struct$item]) # Alpha = .67, G6 = .69, Omega Hier = .48, Omega H asymp = .65, Omega T = .73


# ITEM-RESPONSE THEORY ANALYSIS ----

# set-up formulas
form1 <- list(
  base = bf( score | trials(max) ~ 1 + (1 | item) + (1 | id) ),
  dif = bf( score | trials(max) ~ 0 + grp + (0 + grp | item) + (1 | id) )
)

# set-up weakly informative priors
prior1 <- list(
  
  base = c(
    prior("normal(0,3)", class = "Intercept"),
    prior("normal(0,3)", class = "sd")
  ),
  
  dif = c(
    prior("normal(0,3)", class = "b"),
    prior("normal(0,3)", class = "sd"),
    prior("lkj(2)", class = "cor")
  )
  
)

# fit them
fit1 <- lapply(
  
  setNames( names(prior1), names(prior1) ),
  function(i)
    
    brm( formula = form1[[i]],
         prior = prior1[[i]],
         family = binomial(link = "logit"),
         data = d1,
         seed = 87542,
         control = list(adapt_delta = .9),
         save_pars = save_pars(all = T)
         )
 
)

# PSIS-LOO compare the models
with( fit1, loo(base, dif, moment_match = T) ) # the difference in ELPD is circa 1/2 MCSE in favour of 'base', keeping 'dif' for next analysis

# powerscaling prior-sensitivity
print(powerscale_sensitivity(fit1$dif), n = 390) # item parameters a-ok
powerscale_plot_quantities(fit1$dif, variable = "sd_item__grpctrl") # prior-data conflicte detecting, does not look bad though

# item-specific posterior predictions
pp_check( fit1$dif, type = "bars_grouped", group = "item", ndraws = NULL) +
  labs(
    x = "Score",
    title = "Posterior predictive check of the Multilevel Binomial regression with Differential Item Functioning",
    subtitle = "Bars represent counts of each score as observed in the sample,\npoints and whiskers represent posterior medians and 95% posterior probability intervals as predicted by the model"
  ) +
  scale_x_continuous(breaks = 0:6, labels = 0:6) +
  theme_grey() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5)
  )

# save it
ggsave(
  plot = last_plot(),
  filename = here("figures","ppc.jpg"),
  dpi = 300,
  width = 10,
  height = 9
)

## item parameters ----

# item easiness parameters
fig1a <- fit1$dif %>%
  
  spread_draws(r_item[item,group]) %>%
  median_qi(r_item = -r_item, .width = .95) %>%
  mutate( item =  gsub(".", " ", item, fixed = T) ) %>%
  arrange( match( item, rev(struct$label) ) ) %>%
  mutate(
    item = factor(
      x = item,
      levels = rev(struct$label),
      ordered = T
    ),
    `Group: ` = factor(
      x = case_when(group == "grpctrl" ~ "MoCA", group == "grpexp" ~ "MoCA-WLE")
    )
  ) %>%
  
  ggplot(aes(y = item, x = r_item, xmin = .lower, xmax = .upper, colour = `Group: `) ) +
  geom_vline(xintercept = 0, colour = "red") +
  geom_point(position = position_dodge(width = .4), size = 5) +
  geom_linerange(position = position_dodge(width = .4), linewidth = 1.5) +
  labs(y = NULL, x = "Item difficulty (logit scale)") +
  scale_colour_manual( values = cbPal[c(2,6)] ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# item differences
fig1b <- fit1$dif %>%
  
  spread_draws(r_item[item,group]) %>%
  pivot_wider(names_from = "group", values_from = "r_item") %>%
  median_qi(diff = -grpexp - (-grpctrl), .width = .95) %>%
  mutate( item =  gsub(".", " ", item, fixed = T) ) %>%
  arrange( match( item, rev(struct$label) ) ) %>%
  mutate(
    item = factor(
      x = item,
      levels = rev(struct$label),
      ordered = T
    ),
    target = if_else(item == "Free Recall", T, F)
  ) %>%
  
  ggplot() +
  aes(y = item, x = diff, xmin = .lower, xmax = .upper, colour = target) +
  geom_vline(xintercept = 0, colour = "red") +
  geom_point(size = 5) +
  geom_linerange(linewidth = 1.5) +
  scale_colour_manual( values = c("grey67","black") ) +
  labs(y = NULL, x = "Item difficulty difference (MoCA-WLE-minus-MoCA)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

# organise them
(fig1a | fig1b) + plot_layout(axes = "collect")

# save the results
ggsave(
  plot = last_plot(),
  filename = here("figures","item_parameters.jpg"),
  dpi = 300,
  width = 9,
  height = 10
)

# extract the item difference of Free Recall for abstract
fit1$dif %>%
  
  spread_draws(r_item[domain,group]) %>%
  pivot_wider(names_from = "group", values_from = "r_item") %>%
  median_qi(diff = -grpexp - (-grpctrl), .width = .95) %>%
  filter(domain == "Free.Recall") %>%
  mutate( across( where(is.numeric), ~round(.x,2) ) )

## person parameters ----

# plot it
fit1$dif %>%
  
  spread_draws(b_grpctrl, b_grpexp, r_id[id, ]) %>%
  median_qi( ability = if_else(id %in% subset(d2, grp == "exp")$id, b_grpexp + r_id, b_grpctrl + r_id ) ) %>%
  mutate(`Group: ` = if_else(id %in% subset(d2, grp == "exp")$id, "MoCA-WLE", "MoCA") ) %>%
  arrange(ability) %>%
  mutate( id = seq_len(n() ) ) %>%
  
  ggplot(aes(y = id, x = ability, xmin = .lower, xmax = .upper, colour = `Group: `) ) +
  geom_point(size = 3) +
  geom_linerange(linewidth = 2, alpha = .66) +
  scale_colour_manual( values = cbPal[c(2,6)] ) +
  labs(
    y = "Person number (sorted)",
    x = "Person ability (logit scale)",
    title = "Person ability estimates of the Multilevel Binomial regression with Differential Item Functioning",
    subtitle = "Points and whiskers represent posterior medians and 95% posterior probability intervals of peesons' ability parameters\nsorted from the lowest (bottom) to the highest (top)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid = element_blank(),
    plot.title = element_text(hjust = .5, face = "bold"),
    plot.subtitle = element_text(hjust = .5)
  )

ggsave(
  plot = last_plot(),
  filename = here("figures","person_parameters.jpg"),
  dpi = 300,
  width = 12,
  height = 9
)


