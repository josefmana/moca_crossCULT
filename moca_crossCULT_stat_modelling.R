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
d <- list(
  dem = read.csv( "data/moca_crossCULT_data_demo.csv", sep = ";" ),
  ctr = read.csv( "data/moca_crossCULT_data_ctrl.csv", sep = ";" ),
  exp = read.csv( "data/moca_crossCULT_data_exp.csv", sep = ";" )
)

# for further analyses keep only paired subjects (i.e., those ctrl subjects that have been paired with an exp subject)
d_nanok <- list( incl = d$dem[ d$dem$exp.id != "", ], excl = d$dem[ d$dem$exp.id == "", ] )

# add sum scores wherever missing
for ( i in c("ctr","exp") ) {
  for ( j in c("clock","naming","abstraction","free_recall") ) d[[i]][,j] <- rowSums( d[[i]][ , grepl( paste0(j,"_"), colnames(d[[i]]) ) ] )
}

# keep only the items list in struct
df <- d$ctr[ , colnames(d$ctr) %in% c("id",struct$item) ] %>%
  # add an index of the group and flag subjects to be included
  add_column( grp = "ctrl", .after = 1 ) %>%
  mutate( included = ifelse( id %in% d$dem[ d$dem$exp.id != "", "ctr.id"], 1, 0 ), .after = 2 ) %>%
  # bind the new experimental subjects below the old normative control ones
  rbind.data.frame(
    d$exp[ , colnames(d$exp) %in% c("id",struct$item) ] %>%
      add_column( grp = "exp", .after = 1 ) %>%
      mutate( included = ifelse( id %in% d$dem[ d$dem$exp.id != "", "exp.id"], 1, 0 ), .after = 2 )
  ) %>%
  # add columns for demographic variables
  add_column( age_y = NA, .after = 2 ) %>%
  add_column( edu_y = NA, .after = 3 ) %>%
  add_column( sex = NA, .after = 4 )

# add demographic variables via loops
for ( i in df$id ) {
  
  if ( i %in% d$dem$ctr.id ) df[ df$id == i , c("age_y","edu_y","sex") ] <- d$dem[ d$dem$ctr.id == i, paste0( "ctr.",c("age_y","edu_y","sex") ) ]
  if ( i %in% d$dem$exp.id ) df[ df$id == i , c("age_y","edu_y","sex") ] <- d$dem[ d$dem$exp.id == i, paste0( "exp.",c("age_y","edu_y","sex") ) ]
  if ( i %in% d$dem$excl.id ) df[ df$id == i , c("age_y","edu_y","sex") ] <- d$dem[ d$dem$excl.id == i, paste0( "excl.",c("age_y","edu_y","sex") ) ]

}
