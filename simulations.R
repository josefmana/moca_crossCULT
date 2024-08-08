# set working directory (works only in RStudio)
setwd( dirname(rstudioapi::getSourceEditorContext()$path) )

# list packages to be used
pkgs <- c("dplyr", "tidyverse", # data wrangling
          "ggplot2", "patchwork", # plotting
          "brms", "tidybayes", # stat modelling
          "MASS", "MBESS" # data-generating
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

# read out the structure of the MoCA test
struct <- read.csv( "data/moca_struct.csv", sep = "," )


# ---- generate synthetic observations ---

generate_data <- function( N = 1e5, m = 1.2, sd = 0.5, rho = .5 ) {
  
  # prepare a vector of means and SDs on latent scale and their implications for observed probabilities
  true <- struct %>% mutate(
    # add means and SDs in latent standard normal space
    M = m, SD = sd,
    # calculate to be observed probabilities on raw scales via standard normal cumulative function (i.e., inverse probit)
    prob = pnorm(M), E_score = max_score * prob
  )
  
  # prepare a correlation matrixr egarding MoCA subtasks
  cor <- with( struct, matrix( nrow = length(item), ncol = length(item), dimnames = list( item, item ) ) )
  
  # fill-in with Pearson's correlations
  for( i in rownames(cor) ) {
    for ( j in colnames(cor) ) {
      # r = 1 on the diagonal
      if (i == j) cor[i,j] <- 1
      # identical correlation for the other pairs
      else cor[i,j] <- rho
    }
  }
  
  # simulated a data set of latent scores
  lat <- mvrnorm( N, true$M, cor2cov(cor,true$SD) ) %>% pnorm()
  
  # generate observed values from the latent ones
  obs <- sapply( 1:N,
                 function(i) sapply( true$item, function(j) rbinom( 1, with( true, max_score[item == j] ), lat[i,j] ) )
  ) %>% t() %>% as.data.frame()
  
  # return the observations
  return(obs)
  
}

# ---- do some simulations ----

# data
d <- lapply( c(0,1.2), function(i)
  lapply( c( 0,0.3,0.5,0.8), function(j)
    generate_data( m = as.numeric(i), sd = .5, rho = as.numeric(j) )
  ) )

# grid for plotting
par( mfrow = c(4,2) )

# histograms
for ( i in 1:4 ) {
  for ( j in 1:2 )
    hist( d[[j]][[i]] %>%
            mutate( subtraction = case_when( subtraction %in% c(5,4) ~ 3,
                                             subtraction %in% c(3,2) ~ 2,
                                             subtraction == 1 ~ 1,
                                             subtraction == 0 ~ 0 )
                    ) %>% rowSums(),
          
          xlim = c(0,30), xlab = "MoCA total score", breaks = 20,
          main = paste0( "M = ", c(0,1.2)[j],", rho = ", c(0,0.3,0.5,0.8)[i] )
    )
}
