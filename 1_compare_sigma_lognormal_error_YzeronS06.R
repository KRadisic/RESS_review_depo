##########################################################
#####    GENERATE RAIN PERTURBATIONS W LOGNORMAL DISTRIB
#####    COMPARISON FOR VARIOUS SIGMA

## Path and libraries
rm(list=ls())
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_functions.R")

## Load Yzeron rains
setwd('~/code/6_GSA_stochastic_article/data/raw/FORCING_DATA/new_forcing_data_Yzeron_S06_hourly/')

## Choose one rain as true rain
true_rain_idx <- 53
original_rain_hourly <- read.csv(paste('new_forcing_data_hourly',true_rain_idx,'.in', sep=""), sep ="\t")[1:6,]
head(original_rain_hourly)

## Create dataframe where each row is a error realization, each column the hourly timestep
R <- 500 # number of error realizations

###############################################
####    ERROR DISTRIBUTION CHOICE AND SAMPLING
###############################################
## Sample error from LogNormal distribution, for each timestep
list_of_plots <- list()
z <- 1

for (stddev in c(0.20,0.05,0.04,0.03,0.01)){
  mat_rain_error_YzeronS06 <- t(matrix(original_rain_hourly$Rain..cm.h.))
  for (error_idx in seq(1,R)){
    epsilon <- exp(rnorm(mean = 0, sd = stddev, n=6))
    epsilon
    mat_rain_error_YzeronS06 <- rbind(t(original_rain_hourly$Rain..cm.h.) * epsilon, mat_rain_error_YzeronS06)
    # put the original last, otherwise can't be seen on the plot...
  }
  mat_rain_error_YzeronS06
  ## Format dataframe for plotting
  df_rain_error_YzeronS06 <- data.frame(mat_rain_error_YzeronS06)
  colnames(df_rain_error_YzeronS06) <- c(1,2,3,4,5,6)
  df_rain_error_YzeronS06$error_realization <- row.names(df_rain_error_YzeronS06)
  df_rain_error_YzeronS06$orig_rain <- factor(c(rep('error', R), 'original'), levels = c('error','original'))
  head(df_rain_error_YzeronS06)
  
  dfdf_rain_error_YzeronS06 <- melt(df_rain_error_YzeronS06)
  head(dfdf_rain_error_YzeronS06)
  
  ## Plot the one rain with the newly samples errors (transform in mm only for plot)
  gg_rain_error_YzeronS06 <- ggplot(dfdf_rain_error_YzeronS06, 
                                    aes(x = variable, y = value*10, group = orig_rain)) + 
    geom_point(aes(shape = orig_rain, color = orig_rain)) + theme_bw() + 
    ylab('Rain [mm/h]') + xlab('Time [h]') +
    scale_shape_manual(name = "Rain", values = c(error = 4, original = 19))+
    scale_color_manual(name = "Rain", values = c(error = "#c4c0b3", original = "red"))+
    ggtitle(paste("x * exp N(0, sd=",stddev,")",sep="")) + 
    theme(legend.position="none") +ylim(c(0,25))+
    theme(legend.position="left")
  list_of_plots[[z]] <- gg_rain_error_YzeronS06
  z <- z+1
}

## ggarr  nge the list of plots
gg_rain_error_YzeronS06_varying_sigma <- do.call(grid.arrange, c(list_of_plots, ncol = 5))
gg_rain_error_YzeronS06_varying_sigma
  
