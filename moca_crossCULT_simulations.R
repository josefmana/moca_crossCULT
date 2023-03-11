# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list packages to be used
pkgs <- c("dplyr", "tidyverse", # data wrangling
          "ggplot2", "patchwork", # plotting
          "mvtnorm" # multivariate normal distribution
          )

# load or install each of the packages as needed
for ( i in pkgs ) {
  if ( i %in% rownames( installed.packages() ) == F ) install.packages(i) # install if it ain't installed yet
  if ( i %in% names( sessionInfo()$otherPkgs ) == F ) library( i , character.only = T ) # load if it ain't loaded yet
}

# create folders for models, figures, tables and sessions to store results and sessions info in
# prints TRUE and creates the folder if it was not present, prints NULL if the folder was already present
sapply( c("mods", "figs", "tabs", "sess"), function(i) if( !dir.exists(i) ) dir.create(i) )


# ---- simulate latent scores ----

# list all score types
scores <- c( "tmt_b", "cube", paste( "clock", c("contour","numbers","hands"), sep = "_" ), # visuospatial/executive
             paste( "naming", c("lion","rhino","camel"), sep = "_" ), # naming
             paste( "digit_span", c("forward","backward"), sep = "_" ), "tapping", "subtraction", # attention
             paste( "sentence", c("john","cat"), sep = "_" ), "fluency", # language
             paste( "similarity", c("vehicles","measurements"), sep = "_" ), # abstraction
             paste( "memory", c("face","velvet","church","daisy","red"), sep = "_" ), # delayed recall
             paste( "orientation", c("date","month","year","day","place","city"), sep = "_" ) # orientation
             )
