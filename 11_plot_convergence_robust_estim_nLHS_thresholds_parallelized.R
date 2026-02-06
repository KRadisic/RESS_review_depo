###################################################################
####   EVALUATE EMULATOR IN NEW POINTS TO GET ROBUST MINIMIZERS
###################################################################
rm(list=ls())
set.seed(111222)

## Libraries
library('lhs')
library('EQL')

## Functions : 
# isoprobabilistic transform z to x
IsoProb_z_to_x <- function(z, df_prior){
  z*df_prior$sd[1:5] + df_prior$mean[1:5]
}
# isoprobabilistic transform x to z
isoprobTransfGaussian <- function(x,mean,stddev){
  return((x-mean)/stddev)
}
# evaluate PCE in given point.
eval_surrogate_OLS_in_point_x <- function(x, df_prior, A_multidx, coeff_PCE_fixed_rain) {
  # evaluates the PCE in a fixed 6-dimensional point x
  Z <- mapply(isoprobTransfGaussian, x, df_prior$mean, df_prior$sd)
  #print('evaluate OLS normalized hermite')
  Psi_alpha = c()
  for (idx_alpha in seq(1,dim(A_multidx)[1])){
    alpha <- A_multidx[idx_alpha, ]
    psi_alpha_Z <- mapply(hermite, Z, alpha, prob = TRUE)
    psi_alpha_Z[alpha == 2] <- psi_alpha_Z[alpha == 2]/sqrt(2) # normalize dirty
    psi_alpha_Z[alpha == 3] <- psi_alpha_Z[alpha == 3]/sqrt(2*3) # normalize dirty
    psi_alpha_Z[alpha == 4] <- psi_alpha_Z[alpha == 4]/sqrt(2*3*4) # normalize dirty
    
    Psi_alpha_Z <- prod(psi_alpha_Z)
    Psi_alpha <- cbind(Psi_alpha, Psi_alpha_Z)
  }    
  # exp des x card of basis
  colnames(Psi_alpha) <- paste(A_multidx[seq(1,dim(A_multidx)[1]),1], 
                               A_multidx[seq(1,dim(A_multidx)[1]),2],
                               A_multidx[seq(1,dim(A_multidx)[1]),3], 
                               A_multidx[seq(1,dim(A_multidx)[1]),4],
                               A_multidx[seq(1,dim(A_multidx)[1]),5],
                               A_multidx[seq(1,dim(A_multidx)[1]),6])
  
  J_surr_OLS <- Psi_alpha%*%coeff_PCE_fixed_rain
  return(J_surr_OLS)
}

## sources (for libraries and functions)
#source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")
#source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_functions.R")

## paths local/ cluster change accordingly only working directory
# working_directory <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/"
working_directory <- "/lustre/radisick/CALIB_FWK"
excursion_set_results <- paste(working_directory, "/data/processed/excursion_sets/", sep="")


## constants
vector_input_indices = c(63,68,72,99,71,65); # indices of PESH parameters
n_params <- 5 # Number of parameters to be optimized
Nomega <- 500 # number of trajectiories used
str_beta <- '005' # beta regularization as a character.

## Sizes of the space-filling design in the parameter space; 0 at the beginning
# N_LHS <- c(0, 20)#, 30, 50, 100)#, 150, 200, 300, 400, 500, 750, 1000, 1500, 2000, 
# #2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000) 

## load data
setwd(paste(working_directory, "/scripts_PCE/metamodels/", sep=""))
A_multidx <- read.csv(paste('union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb',str_beta,'.csv',sep=""), header = TRUE)
prior_distrib <- read.csv('prior_params_wo_bds145.csv',header = FALSE)
load(file = paste('c_alpha_KDE_new_beta', str_beta, '.RData', sep = "")) # c_alpha_KDE_new

## Dataframe priors
df_prior <- data.frame(param = c('th9', 'mn10', 'th10', 'th13', 'thr10', 'hg10'),
                       mean = prior_distrib[vector_input_indices,1],
                       sd =   prior_distrib[vector_input_indices,2])
c_alpha_KDE_new_temp <- c_alpha_KDE_new
c_alpha_KDE_new <- c_alpha_KDE_new_temp[1:Nomega,]
c_alpha_KDE_new_multidx_r <- t(c_alpha_KDE_new)


# prevu :2h00 pour 1000x1000 simulations
# 50..5000 x 500 = 14 650 000 simu = 30h.
# 50...6000 x 500 = 40h utorak popodne -> 8 ujutru cetvrtak.
# 50..4500 48h
# Loop on the number of candidate points, starts from 2 to compensate the 0 in the beginning
dfdf_dataframe <- c()
#for (idx_n_lhs in seq(2,length(N_LHS))){
## Initialize sets
#n_lhs <-  #N_LHS[idx_n_lhs]
## read from terminal argument (parallel version)
n_lhs <- commandArgs(trailingOnly = TRUE)[1]
print(paste("size of n_lhs",n_lhs,sep=""))
n_lhs <- as.integer(n_lhs)

print(paste('n_lhs = ', n_lhs))

# generate LHS in the parameter space
z_LHS_quantiles <- randomLHS(n_lhs, n_params)
z_LHS <- qnorm(z_LHS_quantiles, mean = 0, sd = 1)

# Loop on ALL trajectories
df <- c()
for (rain_idx in seq(1:length(c_alpha_KDE_new_multidx_r[1,]))){
  print(paste('omega idx = ', rain_idx))
  
  # PCE coefficients corresponding to this trajectory
  c_alpha <- c_alpha_KDE_new_multidx_r[,rain_idx]
  
  # Loop on the space-filling design in the parameter space
  for (idx_theta in seq(1, n_lhs)){
    
    # Fix one point from the unitary LHS  # Take its value in the parameter space
    zk <- z_LHS[idx_theta,]; theta_k <- c(IsoProb_z_to_x(zk,df_prior = df_prior),-3.50)
    
    # Evaluate the trajectory in this point
    value_J <- eval_surrogate_OLS_in_point_x(theta_k, df_prior, A_multidx, c_alpha_KDE_new_multidx_r[,rain_idx])
    df <- rbind(df, c(theta_k, rain_idx, idx_theta, value_J))
  }
}
df_dataframe <- data.frame(df)
colnames(df_dataframe) <- c(df_prior$param, 'omega_idx', 'idx_theta', 'value_J')

## Add column for the cardinality of X candidates
df_dataframe$n_candidates <- n_lhs
#head(df_dataframe)

# Save (intermediary) results
setwd(excursion_set_results)
save(df_dataframe, file = paste("df_interm_MM_realiz_Nomega",Nomega,"_nLHS",n_lhs,".RData", sep = ""))

## Concatenate with previous cardinalities of the parameter space
#dfdf_dataframe <- rbind(dfdf_dataframe, df_dataframe)
#}

# Save results
#setwd(excursion_set_results)
#save(dfdf_dataframe, file = paste("df_MM_realiz_Nomega",Nomega,"_nLHS",n_lhs,".RData", sep = ""))
