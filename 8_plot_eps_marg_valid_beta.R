## plot eps margi wrt beta
rm(list=ls())
#source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_functions.R")
#source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")
#metamodel_folder <- '~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/scripts_PCE/metamodels/'
#figures_plot <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/figures/"
#excursion_set_results <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/processed/excursion_sets/"

metamodel_folder <- '~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/scripts_PCE/metamodels/'
working_directory <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/"
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_functions.R")
setwd(paste(working_directory, "/scripts_PCE/metamodels/", sep=""))

## fix beta
beta <- 0.000
str_beta <- '000'

prior_distrib <- read.csv('prior_params_wo_bds145.csv',header = FALSE)
vector_input_indices = c(63,68,72,99,71,65);
df_prior <- data.frame(param = c('th9', 'mn10', 'th10', 'th13', 'thr10', 'hg10'),
                       mean = prior_distrib[vector_input_indices,1],
                       sd =   prior_distrib[vector_input_indices,2])

#################################################################
## read the PESHMELBA simulations cost function in the test experimental design of the parameters
setwd(paste(working_directory, "/data/raw/", sep=""))

## read the observation
Y_test  <- npyLoad('Ymoisture_profile_LHS25_100_R198_dim6_YzeronS06_hourly.npy')
true_rain_idx <- 53
Y_true <- Y_test[true_rain_idx*100+1,]

## read validation set
nR_validation <- 299 # number of trajectories used for validation
Y_validation <- npyLoad('Ymoisture_profile_LHS25_100_R499_dim6_YzeronS06_LogN02.npy')[(1+200*100):(499*100),]
x_validation <- npyLoad('LHS25_100_R500_dim6_sample_final.npy')[1:(100*nR_validation),]

###############################
###    VALIDATION CONST FUNCTION
###############################
if (length(Y_validation[1,]) != length(Y_true))
{stop("Error : dimensions not matching")}
delta_Y <- Y_validation - t(matrix(Y_true, length(Y_true),length(Y_validation[,1])))
delta_Y2 <- delta_Y^2
#  weights = c(0.5,1,2,3,4,5,6,10,15,20,25,30,35,40,45,50,55,65,75,100,150,200,250,300,400)
weights <- rep(1,length(Y_validation[1,]))
J_validation <- delta_Y2%*%weights
x_J_validation_rains <- data.frame(th9 = x_validation[,63], 
                                   mn10  = x_validation[,68], 
                                   th10  = x_validation[,72], 
                                   th13  = x_validation[,99], 
                                   thr10 = x_validation[,71], 
                                   hg10  = x_validation[,65], 
                                   J0 = J_validation) # ML order
# Calculate the Jb term,regularization for mn and thetar
Jb.validation <- beta*(((x_J_validation_rains$mn10 - df_prior$mean[2])/df_prior$sd[2])^2 +
                         ((x_J_validation_rains$thr10 - df_prior$mean[5])/df_prior$sd[5])^2)
x_J_validation_rains$J <- x_J_validation_rains$J0 + Jb.validation


###############################################
#####     BOOTSTRAPPING
###############################################
## choose bootstrap number
N_boot <- 100

NX_original_sample <- 1:100  # number of LHS points in parameter space
Ntraj_original_sample <- 1:299 # number of trajectories
NX_num_samples_boot <- 100  # fix the number of (bootstrap) to sample with repetition
Ntraj_num_samples_boot <- 299  # fix the number of (bootstrap) to sample with repetition

vect_params <- 1:5
eps_marg_valid_bootstrap <- c()

## read all trajectories corresponding to bootstrap resample NX_num_samples_boot
for (n_boot in seq(1:NX_num_samples_boot)){
  
  ## bootstrap the validation set in the parameter space
  NX_boot_idxs <- sample(NX_original_sample, size = NX_num_samples_boot, replace = TRUE)
  
  ##  bootstrapping the trajectories
  boot_trajectories1 <- sample(Ntraj_original_sample, size = Ntraj_num_samples_boot, replace = TRUE)
  boot_trajectories2 <- sample(Ntraj_original_sample, size = Ntraj_num_samples_boot, replace = TRUE)
  
  ## keep data for two histograms in each LHS point
  data1_all_histograms <- data2_all_histograms  <- c()
  
  ## loop over the LHS points (as they appear in the boot resample)
  for (idx_test_x_param in NX_boot_idxs){
    
    ## fix the parameter value to a point from the test set rain dataset from 1-nb_x_param (1 is the "unknown" true parameter)
    x_param_fixed <- x_J_validation_rains[idx_test_x_param, vect_params]
    
    # the test histogram in a fixed test point x
    data1 <- x_J_validation_rains[(boot_trajectories1-1)*100 + idx_test_x_param,]
    data2 <- x_J_validation_rains[(boot_trajectories2-1)*100 + idx_test_x_param,]
    
    ## verify that they are identical
    
    ## keep track of histogram in each of the parameter points
    data1_all_histograms <- cbind(data1_all_histograms, data1$J)
    data2_all_histograms  <- cbind(data2_all_histograms, data2$J)   
    
    ## Compare the two sets of data (two histograms) with a KS test.
    kstest <- ks.test(data1$J, data2$J) # p-value small : probability of coming from the same distribution is small
    #print(paste('idx_x = ', idx_test_x_param, 'p.val', kstest$p.value))
    
    ## Plot and save THE TWO HISTOGRAMS
    df <- data.frame(value_J = c(data1$J, data2$J), 
                     type = c(rep('data1', length(data1$J)), rep('data2', length(data2$J))))
    
    gg_hist <- ggplot(df, aes(x = value_J, y =..density.., fill = type, color = type)) + 
      geom_histogram(alpha = 0.2, position="identity", binwidth = 0.001) + theme_bw()  + xlim(c(0,0.15)) +  
      ggtitle(paste('Comparison of histograms in test_point_idx = ',idx_test_x_param,'; KS p-val = ', trunc(100*kstest$p.value)/100,sep=''))
    gg_hist
    
    ## save figures
    #setwd(paste(working_directory, "/output_figures", sep=""))
    #png(paste('histogram_bootstrap_test_rains_test_parameter_LogN02_idx', idx_test_x_param, '.png', sep=""))
    #print(gg_hist)
    #dev.off()
  }
  
  ###     MAKE DATA FRAME WITH MARGINAL HISTOGRAMS IN VALIDATION POINTS
  df_marginals <- data.frame(rbind(data1_all_histograms, data2_all_histograms))
  ## each column is a parameter point (100 columns), the rows correspond to trajectories
  colnames(df_marginals) <- NX_boot_idxs
  df_marginals$type <- as.factor(c(rep('data1',length(data1_all_histograms[,1])),
                                   rep('data2',length(data2_all_histograms[,1]))))
  head(df_marginals)
  
  ##################################
  ###    CALCULATE EPS MARGINAL (for one bootstrap realization of LHS params)
  ##################################
  # Loop over all validation points
  vector_dWS_valid_points <- c() 
  vector_stddev_X_valid_points <- c()
  N_X <- length(NX_boot_idxs)
  
  for (x_idx in NX_boot_idxs){
    # Calculate the histograms
    hist_data1 <- hist(df_marginals[df_marginals$type=='data1', x_idx], breaks = 30, plot = FALSE)
    hist_data2 <- hist(df_marginals[df_marginals$type=='data2', x_idx], breaks = 30, plot = FALSE)
    stddev_X <- sd(df_marginals[df_marginals$type=='data2', x_idx])
    
    # Extract the histogram counts and break points
    counts1 <- hist_data1$counts
    breaks1 <- hist_data1$breaks
    
    counts2 <- hist_data2$counts
    breaks2 <- hist_data2$breaks
    
    # Compute the cumulative sum of the counts
    cumulative_counts1 <- cumsum(counts1)
    cumulative_counts2 <- cumsum(counts2)
    
    # Normalize the cumulative sum to get the CDFs
    cdf_values1 <- cumulative_counts1 / sum(counts1)
    cdf_values2 <- cumulative_counts2 / sum(counts2)
    
    # Create data frames for the CDFs
    cdf_data1 <- data.frame(
      x = breaks1[-length(breaks1)],  # Use the left edge of each bin
      cdf = cdf_values1
    )
    
    cdf_data2 <- data.frame(
      x = breaks2[-length(breaks2)],  # Use the left edge of each bin
      cdf = cdf_values2
    )
    wasserstein_distance <- wasserstein1d(cdf_data1$cdf, cdf_data2$cdf, p = 2)
    #print(wasserstein_distance)
    # save all dWS in all valid points
    vector_dWS_valid_points <- cbind(vector_dWS_valid_points, wasserstein_distance)
    # save all stddev in all valid points
    vector_stddev_X_valid_points <- cbind(vector_stddev_X_valid_points, stddev_X)
  }
  eps_marg <- NaN
  if (length(vector_dWS_valid_points) == N_X & length(vector_stddev_X_valid_points) == N_X){
    eps_marg <- mean(vector_dWS_valid_points/vector_stddev_X_valid_points)
  }
  
  # keep eps marg of each bootstrap realization in a vector
  eps_marg_valid_bootstrap <- c(eps_marg_valid_bootstrap, eps_marg)
}

# save the eps margi of bootstrapped validation
setwd(metamodel_folder)
save(eps_marg_valid_bootstrap, file = paste('eps_marg_valid_bootstrap_beta',str_beta,'.RData', sep =''))

quantile_95 <- quantile(eps_marg_valid_bootstrap, probs = 0.95)
quantile_05 <- quantile(eps_marg_valid_bootstrap, probs = 0.05)

quantile_25 <- quantile(eps_marg_valid_bootstrap, probs = 0.25)
quantile_75 <- quantile(eps_marg_valid_bootstrap, probs = 0.75)

#################################
######       PLOTTING       #####
#################################
eps_marg_beta<- c()

# read data
setwd(metamodel_folder)
for (beta in c(0,0.002,0.005,0.01)){
  load(file = paste("eps_marg_beta", beta, ".RData"))
  eps_marg_beta <- rbind(eps_marg_beta, c(beta, eps_marg)) 
}

# transform to data frame
df_eps_marg_beta<- data.frame(eps_marg_beta)
colnames(df_eps_marg_beta) <- c('beta', 'eps_marg')
head(df_eps_marg_beta)

# plot
ggplot(df_eps_marg_beta, aes(x = beta, y = eps_marg)) + geom_point() + theme_bw() + 
  geom_rect(aes(ymin=min(quantile_05), ymax=max(quantile_95), xmin=-Inf, xmax=Inf), alpha=0.1) + 
  geom_rect(aes(ymin=min(quantile_25), ymax=max(quantile_75), xmin=-Inf, xmax=Inf), alpha=0.2) + 
  geom_hline(yintercept = mean(eps_marg_valid_bootstrap), color = 'red') # choose the one corresponding to beta
