#############################################
#####     SEE REALIZATIONS OF MM STOCHASTIC 
#################################################
rm(list=ls())

working_directory <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/"
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_functions.R")
setwd(paste(working_directory, "/scripts_PCE/metamodels/", sep=""))
prior_distrib <- read.csv('prior_params_wo_bds145.csv',header = FALSE)
vector_input_indices = c(63,68,72,99,71,65);
df_prior <- data.frame(param = c('th9', 'mn10', 'th10', 'th13', 'thr10', 'hg10'),
                       mean = prior_distrib[vector_input_indices,1],
                       sd =   prior_distrib[vector_input_indices,2])

################################################
###    LOAD THE ROBUST CALIB RESULTS
processed_data_directory <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/processed/"
setwd(processed_data_directory)

## Number of simulations of the metamodel for choosing robust calibration
Nomega <- 500
n_lhs <- 4500
beta <- 0.005
idx_of_conditional_minimizer <- 174 
minimizers <- npyLoad(paste("LHS_robust_minimizers_LogN02_beta005_condmin",
                            idx_of_conditional_minimizer,
                            "_Nomega",Nomega,"_nlhs",n_lhs,
                            "_dim6_sample_final.npy", sep=""))

minimizer_names_ordered <- c('$\\hat{x}^{*(r_i)}_{s\\_PCE}$',
                             '$\\bar{x}^{*}$',
                             '${x}^*_{\\hat{f}_s<0.01}$',
                             '${x}^*_{\\hat{f}_s<0.02}$',                             
                             '$\\mu_{prior}^{*}$')

names_params <- c("$\\theta_{s.inter}$", "$mn_{deep}$",
                  "$\\theta_{s.deep}$","$\\theta_{s.surf}$", "$\\theta_{r.deep}$")
df_minimizers <- data.frame(minimizers[c(1,2,3,4,5),vector_input_indices[1:5]])
colnames(df_minimizers) <- names_params
df_minimizers$robust.level <- minimizer_names_ordered
head(df_minimizers)

str_beta <- '005' # beta regularization as a character.
## load data
setwd(paste(working_directory, "/scripts_PCE/metamodels/", sep=""))
A_multidx <- read.csv(paste('union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb',str_beta,'.csv',sep=""), header = TRUE)
prior_distrib <- read.csv('prior_params_wo_bds145.csv',header = FALSE)
load(file = paste('c_alpha_KDE_new_beta', str_beta, '.RData', sep = "")) # c_alpha_KDE_new


# A unique variable to be compatible with optimr package
surrogate_cost_function_analytical_fixedrain_1input <- function(x)
{eval_surrogate_OLS_in_point_x(c(x, -3.50), df_prior, A_multidx, coeff_PCE_fixed_rain)}

# Number of parameters to be optimized
n_params <- 1
param_idxs_out_of_5 <- 1 # 1:5

# Isoprobabilistic transform
IsoProb_z_to_x <- function(z, df_prior){
  z*df_prior$sd[param_idxs_out_of_5] + df_prior$mean[param_idxs_out_of_5]
}

## load KDE
Nomega <- 50
c_alpha_KDE_new_temp <- c_alpha_KDE_new
c_alpha_KDE_new <- c_alpha_KDE_new_temp[1:Nomega,]
c_alpha_KDE_new_multidx_r <- t(c_alpha_KDE_new)

## Initialize sets
n.params <- 5
n_lhs_plot <- 60

# generate LHS in the parameter space
z_LHS_quantiles <- maximinLHS(n_lhs_plot, n.params, dup = 1)
unif_vec_idxs <- c(1:5); gauss_vec_idxs <-c(1:5) # dirty fix distrib. types
z_LHS <- matrix(NaN,n_lhs_plot,n.params)
# Transform to N(0,1)
z_LHS[,gauss_vec_idxs] <- qnorm(z_LHS_quantiles[,gauss_vec_idxs], mean = 0, sd = 1)#

## choose the fixed point for all other parameters
fixed_point <- df_minimizers[df_minimizers$robust.level=='${x}^*_{\\hat{f}_s<0.01}$',1:5]
fixed_point
fixed_point_name <- 'exset01'

################################################################
######   EVALUATE STOCH METAMODEL IN THE FIXED POINT, 
##               vary params one by one
#################################################################

dfdf <- c()
dfdf_pesh <- c()
list_of_plots_excset <- list(); z <- 1

for (param_idxs_out_of_5 in c(1:n.params)){
  print(paste(param_idxs_out_of_5, 'out of ',n.params))
  df <- c()
  # Loop on ALL trajectories
  for (rain_idx in seq(1:length(c_alpha_KDE_new_multidx_r[1,]))){
    print(paste(rain_idx, 'out of ',length(c_alpha_KDE_new_multidx_r[1,])))
    
    # Read pce coefficients corresponding to this trajectory
    c_alpha <- c_alpha_KDE_new_multidx_r[,rain_idx]
    
    # Loop on the space-filling design in the parameter space
    for (idx_theta in seq(1, n_lhs_plot)) {
      
      # Fix one point from the unitary LHS  # Take its value in the parameter space
      zk <- z_LHS[idx_theta,];
      # read the fixed point
      theta_k <- c(as.numeric(fixed_point), df_prior$mean[6])
      # change only chosen params.
      #theta_k[param_idxs_out_of_5] <- IsoProb_z_to_x(zk,df_prior = df_prior)
      theta_k[param_idxs_out_of_5] <- IsoProb_z_to_x(zk,df_prior = df_prior)[param_idxs_out_of_5]
      # Evaluate the trajectory in this point
      value_J <- eval_surrogate_OLS_in_point_x(theta_k, df_prior, A_multidx, c_alpha_KDE_new_multidx_r[,rain_idx])
      
      # Add line in dataframe to be traced
      df <- rbind(df, c(rain_idx, theta_k[1:n.params], value_J))
    }
  }
  colnames(df) <- c('rain.idx', names_params, 'value.J')
  ## Keep all values, for then running PESHM.
  dfdf_pesh <- rbind(dfdf_pesh, df)
  
  ## Keep only the varying values of the parameter
  dfdf <- rbind(dfdf, df[,c(1, (param_idxs_out_of_5 + 1), 7)])
}

if(length(dfdf[,1])/(Nomega*n_lhs_plot) == 5){print('OK')}
dfdf_results <- data.frame(dfdf)
colnames(dfdf_results) <- c('rain.idx', 'variable', 'value.J')

## add varying parameter name
dfdf_results$parameter.name <- rep(names_params, each = Nomega*n_lhs_plot)
head(dfdf_results)
str(dfdf_results)

# Convert variables to the appropriate types if necessary
dfdf_results$variable <- as.numeric(dfdf_results$variable)
dfdf_results$value.J <- as.numeric(dfdf_results$value.J)
dfdf_results$parameter.name <- factor(dfdf_results$parameter.name, 
                                      levels = c("$\\theta_{s.inter}$",
                                                 "$mn_{deep}$","$\\theta_{s.deep}$","$\\theta_{s.surf}$", "$\\theta_{r.deep}$"))# levels = ASI_order)
dfdf_results$rain.idx <- as.factor(as.numeric(dfdf_results$rain.idx))
str(dfdf_results)

dfdf_minimizers <- melt(df_minimizers, id.vars = 'robust.level')
colnames(dfdf_minimizers) <- c('robust.level', 'parameter.name', 'value.calib')
dfdf_minimizers$robust.level <- as.factor(dfdf_minimizers$robust.level)

head(dfdf_results)
head(dfdf_minimizers)
dfdf_results <- merge(x = dfdf_results, y = dfdf_minimizers,
                      by.x = 'parameter.name', by.y = 'parameter.name')

######################################
####    GG PLOT
colors_minimizers <- c('${x}^*_{\\hat{f}_s<0.01}$'= "#386cb0",'${x}^*_{\\hat{f}_s<0.02}$' = "#f0027f")
minimizers_to_plot <-  c('${x}^*_{\\hat{f}_s<0.01}$','${x}^*_{\\hat{f}_s<0.02}$')

df_plot <- dfdf_results[dfdf_results$robust.level %in% minimizers_to_plot,]

library(dplyr)
df_plot <- group_by(df_plot, parameter.name, rain.idx)
df_plot_none <- filter(df_plot, rep(all(value.J >= 0.02),n()))
df_plot_mid <- filter(df_plot, rep(all(value.J >= 0.01) & any(value.J < 0.02),n()))
df_plot_bottom <- filter(df_plot, rep(any(value.J < 0.01),n()))

df_plot_none$color.seuil <- '$\\hat{f}_s>0.02$'
df_plot_mid$color.seuil <- '$\\hat{f}_s<0.02$'
df_plot_bottom$color.seuil <- '$\\hat{f}_s<0.01$'

colors_minimizers
#colors_seuil <- c(df_plot_none$color.seuil = 'gray25'
#                  df_plot_mid$color.seuil = '#f0027f'
#                  df_plot_bottom$color.seuil = '#386cb0')

df_plot <- rbind(df_plot_none, df_plot_mid, df_plot_bottom)

colors_minimizers <- c('${x}^*_{\\hat{f}_s<0.01}$'= "#386cb0",'${x}^*_{\\hat{f}_s<0.02}$' = "#f0027f",
                       '$\\hat{f}_s>0.02$' = 'gray25', '$\\hat{f}_s<0.02$' = '#fe7abf', '$\\hat{f}_s<0.01$' = '#709ad1')


gg <- ggplot(data = df_plot, aes(x = variable, y = value.J, group = rain.idx, color = color.seuil)) + 
  geom_line(alpha = 1) + theme_bw() +
  geom_vline(aes(color = robust.level, xintercept = value.calib),linetype="dashed",linewidth = 1.0)+
  facet_wrap(vars(parameter.name), scales = "free_x", ncol = 3) + xlab('') + ylab('') +
  scale_color_manual('',values=colors_minimizers) + 
  geom_hline(yintercept = 0.01, color="#386cb0",linetype="solid",linewidth = 1.0)+
  geom_hline(yintercept = 0.02, color="#f0027f",linetype="solid",linewidth = 1.0)+
  theme(legend.position = "bottom")   
gg

# louis
gg <- ggplot() + 
  geom_line(data = df_plot, aes(x = variable, y = value.J, group = rain.idx, color = color.seuil), alpha = 0.7) + theme_bw() +
  geom_vline(data = df_plot, aes(color = robust.level, xintercept = value.calib),linetype="dashed", linewidth = 1.0)+
  facet_wrap(vars(parameter.name), scales = "free_x", ncol = 3) + xlab('') + ylab('') +
  scale_color_manual('',values=colors_minimizers) + 
  geom_hline(yintercept = 0.01, color="#386cb0",linetype="solid",linewidth = 1.0)+
  geom_hline(yintercept = 0.02, color="#f0027f",linetype="solid",linewidth = 1.0)+
  theme(legend.position = "bottom")   
gg

######################################
####    SAVE PLOT to tikZ
#setwd("~/Documents/REDACTION/manuscrit_radisic/5_CALIB_MOIST_PROF/figures/")
#tikz(paste('MM_stoch_beta',str_beta,'_moistprof_fix_',fixed_point_name,'.tex', sep = ''), standAlone = TRUE, width=5.5, height=4.0)
#gg
#dev.off()
#tools::texi2dvi(paste('MM_stoch_beta',str_beta,'_moistprof_fix_',fixed_point_name,'.tex', sep = ''),pdf=T)

###################################################################
###    SAVE EXPDESIGN RESULTS FOR NEW PESHMELBA SIMULATIONS
###################################################################
head(dfdf_pesh)
length(dfdf_pesh[,1]) # 50*60*5 (n_omega x n_lhs x n.params)
length(dfdf_pesh[1,]) 

## Arrange minimizers for new PESHMELBA simulations
r.new.rains <- 50 # number of rains for same parameter values
n.dim <- 5 * 60 # number of parameter values per rain
if (length(dfdf_pesh[,1]) != r.new.rains * n.dim){print('error dimensions')}

## read previous LHS (for other parameter values)
raw_data_directory <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/raw/"
setwd(raw_data_directory)

## load an old LHS for the format
old_LHS <- npyLoad("LHS_rep5_50_R500_dim6_sample_final.npy")
## Take only the dimension needed for varying the parameter values, 
# this will then berepeated for simulations under varying rains in peshmelba
single_LHS <- old_LHS[1:n.dim, 1:145]

## Write the new values for the exp design for new peshmelba simulations
# fix only one rain (no matter which one, they all have the same exp design)
single_LHS[,vector_input_indices[1]] <- dfdf_pesh[dfdf_pesh[,1] == 1, 2]
single_LHS[,vector_input_indices[2]] <- dfdf_pesh[dfdf_pesh[,1] == 1, 3]
single_LHS[,vector_input_indices[3]] <- dfdf_pesh[dfdf_pesh[,1] == 1, 4]
single_LHS[,vector_input_indices[4]] <- dfdf_pesh[dfdf_pesh[,1] == 1, 5]
single_LHS[,vector_input_indices[5]] <- dfdf_pesh[dfdf_pesh[,1] == 1, 6]
single_LHS[,vector_input_indices[6]] <- df_prior$mean[6] # fix the 6th parameter

## Check that the parameter values are in the good places 
head(single_LHS[,vector_input_indices])

## Repeat the same exp design as many times as there are new rains
training_set_r <- c()
for (idx_r in 1:r.new.rains)
{
  training_set_r <- rbind(training_set_r, single_LHS)
}

## Save new LHS
#processed_data_directory <- '~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/processed/'
#setwd(processed_data_directory)
#npySave(paste("LHS_MMtrajectories_LogN02_beta005_condmin",idx_of_conditional_minimizer,
#              "_Nomega",Nomega,"_nlhs",n_lhs,
#              "_dim6_sample_final.npy", sep=""), training_set_r)

#################################################################
## read the PESHMELBA simulations cost function in the test experimental design of the parameters
setwd(paste(working_directory, "/data/raw/", sep=""))
Y_test  <- npyLoad('Ymoisture_profile_LHS25_100_R198_dim6_YzeronS06_hourly.npy')
true_rain_idx <- 53
Y_true <- Y_test[true_rain_idx*100+1,]

Y_new_traj <- npyLoad("Ymoisture_profile_MMs_traj_dim6_YzeronS06_LogN02.npy")

setwd(paste(working_directory, "/data/processed/", sep=""))
x_new_traj <- npyLoad('LHS_MMtrajectories_LogN02_beta005_condmin174_Nomega50_nlhs4500_dim6_sample_final.npy')[1:6000,] # plutot 50

###############################
## Calculate cost function
if (length(Y_new_traj[1,]) != length(Y_true))
{stop("Error : dimensions not matching")}
delta_Y <- Y_new_traj - t(matrix(Y_true, length(Y_true),length(Y_new_traj[,1])))
delta_Y2 <- delta_Y^2
#  weights = c(0.5,1,2,3,4,5,6,10,15,20,25,30,35,40,45,50,55,65,75,100,150,200,250,300,400)
weights <- rep(1,length(Y_new_traj[1,]))
J_train <- delta_Y2%*%weights
x_J_train_rains <- data.frame(th9   = x_new_traj[,63], 
                              mn10  = x_new_traj[,68], 
                              th10  = x_new_traj[,72], 
                              th13  = x_new_traj[,99], 
                              thr10 = x_new_traj[,71], 
                              hg10  = x_new_traj[,65], 
                              J0 = J_train) # ML order
# Calculate the Jb term,regularization for mn and thetar
Jb.train <- beta*(((x_J_train_rains$mn10 - df_prior$mean[2])/df_prior$sd[2])^2 +
                    ((x_J_train_rains$thr10 - df_prior$mean[5])/df_prior$sd[5])^2)
x_J_train_rains$J <- x_J_train_rains$J0 + Jb.train

head(x_J_train_rains)
head(dfdf_pesh)


