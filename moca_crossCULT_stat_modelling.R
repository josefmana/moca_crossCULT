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

# add sum scores wherever missing
for ( i in c("ctr","exp") ) {
  for ( j in c("clock","naming","abstraction","free_recall") ) d[[i]][,j] <- rowSums( d[[i]][ , grepl( paste0(j,"_"), colnames(d[[i]]) ) ] )
}

# keep only the items list in struct
df <- d$ctr[ , colnames(d$ctr) %in% c("id",struct$item) ] %>%
  mutate( grp = "ctrl" ) %>%
  rbind.data.frame( d$exp[ , colnames(d$exp) %in% c("id",struct$item) ] %>% mutate( grp = "exp" ) )
