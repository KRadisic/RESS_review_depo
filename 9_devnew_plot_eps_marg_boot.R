##########################################################
###     READ C_ALPHA FROM TRAINING SET
###     INFER DISTRIBUTION OF C_ALPHA
###     PRODUCE NEW SAMPLES C_ALPHA_NEW
###     READ THE VALIDATION C_ALPHA
###     COMPARE : LATENT SPACES, PCE COEFFS AND MARGINAL DISTIBS
###     TO DO : CALCULATE GLOBAL QUALITY OF MARGINAL DISTRIBUTION
###     TO DO : COMPARE QUALITY OF MARGINAL DISTIBUTIONS FOR DIFFERENT ZONES.

rm(list=ls())
library('ggplot2')
library('GGally')
library('ks') # kde
library('EQL') # hermite
library('RcppCNPy') # open .npy data
library('transport')  # For Wasserstein distance calculation
library('boot')  # For  bootstrapping

metamodel_folder <- '~/Documents/RESS_review_depo/metamodels/'
working_directory <- '~/Documents/RESS_review_depo/'

source("~/Documents/RESS_review_depo/0_functions.R")

###########################################
## FUNCTION reproject the PCA scores back to the original basis
pca.recon <- function(pca, PC_scores, n_pc){
  mu <- matrix(rep(pca$center, nrow(PC_scores)), nrow = nrow(PC_scores), byrow = T)
  recon <- PC_scores[,1:n_pc] %*% t(pca$rotation[,1:n_pc]) + mu
  return(recon)
}

## Open the c_alpha train and validation from the OLS
setwd(metamodel_folder)
beta <- 0.005;
str_beta <- '005'; # 00, 002, 005, 01
ml_coeff_multidx_r_TRAIN <- read.csv(paste("ml_multicoeff_c_alpha_OLS_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb",str_beta,".csv", sep=""), header = FALSE)
ml_coeff_multidx_r_TEST <- read.csv(paste("ml_multicoeff_c_alpha_OLS_TEST_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb",str_beta,".csv", sep=""), header = FALSE)

coeff_multidx_r_TRAIN <- data.frame(t(ml_coeff_multidx_r_TRAIN)) # 200 realizations of c_alpha
coeff_multidx_r_TEST <- data.frame(t(ml_coeff_multidx_r_TEST)) # 500 realizations of c_alpha

## Open the basis A_multidx from ML
A_multidx <- read.csv(paste('union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb',str_beta,'.csv',sep=""), header = TRUE)

## Get the degrees of the parameters in the basis
param_order_A_basis <- paste(data.frame(t(A_multidx)))
param_order_A_basis # to be used for c_alpha column names

# Define test data indices 
test_rain_idx <- c(201:499)
calpha_test <- coeff_multidx_r_TEST[test_rain_idx,]
head(calpha_test)
#colnames(calpha_test) <- param_order_A_basis

# open priors from Emilie marginals
setwd(working_directory)
prior_distrib <- read.csv('prior_params_wo_bds145.csv',header = FALSE)
vector_input_indices = c(63,68,72,99,71,65);
df_prior <- data.frame(param = c('th9', 'mn10', 'th10', 'th13', 'thr10', 'hg10'),
                       mean = prior_distrib[vector_input_indices,1],
                       sd =   prior_distrib[vector_input_indices,2])

#################################################################
## read the PESHMELBA simulations cost function in the test experimental design of the parameters

## read the observation
Y_test  <- npyLoad('Ymoisture_profile_LHS25_100_R198_dim6_YzeronS06_hourly.npy')
true_rain_idx <- 53
Y_true <- Y_test[true_rain_idx*100+1,]

## read train and validation set
nR_train <- 200 # number of trajectories used for training the c_alpha
nR_validation <- 299 # number of trajectories used for validation

Y_train <- npyLoad('Ymoisture_profile_LHS_rep5_50_R500_dim6_YzeronS06_LogN02.npy')[1:(50*nR_train),]
Y_validation <- npyLoad('Ymoisture_profile_LHS25_100_R499_dim6_YzeronS06_LogN02.npy')[(1+200*100):(499*100),]

x_train <- npyLoad('LHS_rep5_50_R500_dim6_sample_final.npy')[1:(50*nR_train),] # plutot 50
x_validation <- npyLoad('LHS25_100_R500_dim6_sample_final.npy')[1:(100*nR_validation),]

###############################
###    VALIDATION COST FUNCTION
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

#######################################################
##      choose bootstrap number (for trajectories)
#######################################################
Nboot <- 8 #max(ml_coeff_multidx_r_TEST_boot[33,]) # maybe only 8 ... #50
#Ntraj_original_sample <- 1:100 # 1:200 # number of trajectories
# no need for this as the bootstrapping is done in ML
Nx_num_samples_boot <- c(50, 100) # c(50,100,150,200)

vect_params <- 1:5
df_eps_marg_train_bootstrap <- c()

#n_pc <- 2
#Nx_valid <- 100
#n_boot <- 1

# Loop over the inference type (KDE indep., joint, # of PCs)
for (n_pc in c(2,3,4))
{
  # Loop over Nx = 25, 50, 75
  for (Nx_valid in Nx_num_samples_boot){
    eps_marg_train_traj_bootstrap <- c()
    print(paste('Nx_valid = ', Nx_valid))
    
    setwd(working_directory)
    
    ml_coeff_multidx_r_TEST_boot <- read.csv(paste("ml_boot_multicoeff_c_alpha_OLS_TESTADD_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb",
                                                   str_beta,"Nx_valid",Nx_valid,".csv", sep=""), header = FALSE)
    calpha_train_temp <- data.frame(t(ml_coeff_multidx_r_TEST_boot)) 
    
    ## read all trajectories corresponding to bootstrap resample NX_num_samples_boot
    for (n_boot in seq(1:Nboot)){
      print(paste('n_boot = ', n_boot, 'out of Nboot = ', Nboot))
      
      calpha_train <- calpha_train_temp[calpha_train_temp$X33 == n_boot,1:32]
      
      ##  bootstrapping the indices of training trajectories
      #boot_trajectories1 <- sample(Ntraj_original_sample, size = Nx_valid, replace = TRUE)
      train_rain_idx <- c(1:81,83:200) #boot_trajectories1 #c(1:81,83:200)  # Define training data indices
      #calpha_train <- coeff_multidx_r_TRAIN[train_rain_idx,]
      #head(calpha_train)
      ##colnames(calpha_train) <- param_order_A_basis
      

      ###   PCATAKE A LOWER NUMBER OF N_PC
      pca <- prcomp(calpha_train)
      df_pca <- as.data.frame(pca$x)
      c_alpha_hat_training_reconstructed <- pca.recon(pca, pca$x, n_pc) 
      
      ###     KERNEL DENSITY ESTIMATOR
      kde_fitted <- kde(df_pca[,1:n_pc]) #, h= 0.1
      #kde_fitted1 <- kde(df_pca[,1]) #, h= 0.1
      #kde_fitted2 <- kde(df_pca[,2]) #, h= 0.1
      
      # Simulate new points with KDE
      new_pc_points <- rkde(200, kde_fitted)
      #new_pc_points1 <- rkde(200, kde_fitted1)
      #new_pc_points2 <- rkde(200, kde_fitted2)
      df_new_points <- data.frame(new_pc_points)
      #df_new_points <- data.frame(new_pc_points1,new_pc_points2)
      head(df_new_points)
      
      #################################################
      ###  COMPARE NEW KDE POINTS REPROJECTED ON ORIGINAL BASIS
      ## these are the new c_alpha from the metamodel
      c_alpha_KDE_new <- pca.recon(pca, as.matrix(df_new_points[,1:n_pc]), n_pc)
      c_alpha_KDE_new <- data.frame(c_alpha_KDE_new)
      #colnames(c_alpha_KDE_new) <- colnames(calpha_train)

      #################################################################
      data1_all_test_histograms <- c() # here we will keep the results of the hypothesis testing
      data2_all_metamodel_histograms <- c() # here we will keep the results of the hypothesis testing
      
      #################################################################
      ## check nb_x_param
      nb_x_param <- 100 # number of points in the parameter space in which histograms are observed
      dim_x_param <- 6 # dimension of the parameter space
      vect_params <- 1:dim_x_param
      idx_test_x_param <-6
      
      ###############################
      ###    VALIDATION COST FUNCTION
      idx_test_x_param <- 1
      # loop over the 100 different points in the parameter space, one at a time
      for (idx_test_x_param in seq(1,nb_x_param)){
        ## fix the parameter value to a point from the test set rain dataset from 1-nb_x_param (1 is the "unknown" true parameter)
        x_param_fixed <- x_J_validation_rains[idx_test_x_param, vect_params]
        # read the corresponding cost function realizations
        df_test_J_histogram_in_fixed_point <- x_J_validation_rains[idx_test_x_param + nb_x_param*seq(0,nR_validation) , ] # choose from 1-100 (1 is the true parameter though)
        # verify it is the same point in the parameter space repeated as many time as R_train + R_test
        if (table(table(df_test_J_histogram_in_fixed_point[vect_params]))[[1]] != 1)
        {print("error : the test set is not evaluated in a fix point of the parameter space")}
        # the test histogram in a fixed test point x # keep only the realizations of the test set cost function
        data1 <- test_J_histogram_in_fixed_point <- df_test_J_histogram_in_fixed_point$J[train_rain_idx]
        # the metamodel histogram in this same point # in the fixed point x_param_fixed evaluate all PCE as with all c_alpha
        data2 <- new_J_histogram_in_fixed_point <- eval_surrogate_OLS_in_point_x_variable_rain(x_param_fixed, df_prior, A_multidx, 
                                                                                               t(c_alpha_KDE_new)) 
        # about 0.01 sec per simulation, 12 sec per 1000 simu, 1200 sec = 60min = 1 hour  per 1000*100 simu
        
        data1_all_test_histograms <- cbind(data1_all_test_histograms, data1)
        data2_all_metamodel_histograms  <- cbind(data2_all_metamodel_histograms, data2)
        
        ## Plot and save THE TWO HISTOGRAMS and the TWO SAMPLE CUMULATIVE DENSITY FUNCTIONS FOR COMPARISON
        df <- data.frame(value_J = c(test_J_histogram_in_fixed_point, new_J_histogram_in_fixed_point), 
                         type = c(rep('real model with test rains', length(data1)), rep('stochastic MM', length(data2))))
        
        gg_hist <- ggplot(df, aes(x = value_J, y =..density.., fill = type, color = type)) + 
          geom_histogram(alpha = 0.2, position="identity", binwidth = 0.001) + theme_bw()  + xlim(c(0,0.15)) +  
          ggtitle(paste('Comparison of histograms in test_point_idx = ',idx_test_x_param,sep=''))
        
        ## save figures
        setwd(paste(working_directory, "/output_figures_bootstrap", sep=""))
        png(paste('histogram_bootstrap_test_rains_test_parameter_LogN02_idx', idx_test_x_param, 
                  'Nboot', n_boot, 'Nx_valid', Nx_valid, '.png', sep=""))
        print(gg_hist)
        dev.off()
      }

      ######################################################################
      ###     MAKE DATA FRAME WITH MARGINAL HISTOGRAMS IN VALIDATION POINTS
      ######################################################################
      df_marginals <- data.frame(rbind(data1_all_test_histograms, data2_all_metamodel_histograms))
      colnames(df_marginals) <- paste("$f\\_s({x}\\_{", seq(1,length(df_marginals[1,])),"},\\cdot)$",sep = "")
      # ", d{WS} = ", trunc(1000*vector_dWS_valid_points[idxs_chosen_valid_points])/1000,
      df_marginals$type <- as.factor(c(rep('test',length(data1_all_test_histograms[,1])),
                                       rep('KDE',length(data2_all_metamodel_histograms[,1]))))
      
      ##################################
      ###    CALCULATE EPS MARGINAL
      ##################################
      # Loop over all validation points
      vector_dWS_valid_points <- c() 
      vector_stddev_X_valid_points <- c()
      N_X <- 100
      for (x_idx in seq(1,N_X)){
        # Calculate the histograms
        hist_data1 <- hist(df_marginals[df_marginals$type=='KDE', x_idx], breaks = 30, plot = FALSE)
        hist_data2 <- hist(df_marginals[df_marginals$type=='test', x_idx], breaks = 30, plot = FALSE)
        stddev_X <- sd(df_marginals[df_marginals$type=='test', x_idx])
        
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
        print(wasserstein_distance)
        # save all dWS in all valid points
        vector_dWS_valid_points <- cbind(vector_dWS_valid_points, wasserstein_distance)
        # save all stddev in all valid points
        vector_stddev_X_valid_points <- cbind(vector_stddev_X_valid_points, stddev_X)
      }
      eps_marg <- NaN
      if (length(vector_dWS_valid_points) == N_X & length(vector_stddev_X_valid_points) == N_X){
        eps_marg <- mean(vector_dWS_valid_points/vector_stddev_X_valid_points)
      }
      
      print(paste('eps_marg = ', eps_marg))
      eps_marg_train_traj_bootstrap <- c(eps_marg_train_traj_bootstrap, eps_marg)
      print(paste('beta = ', beta))
      
    } # close loop for calculating eps marg bootstraps on a fixed number of trajectories
    df_temp <- data.frame(eps_marg = eps_marg_train_traj_bootstrap, 
                          num_traj = rep(Nx_valid, Nboot))
    df_eps_marg_train_bootstrap <- rbind(df_eps_marg_train_bootstrap, df_temp)
  }# close loop over all trajectories numbers
  
  setwd(metamodel_folder)
  save(df_eps_marg_train_bootstrap, file = paste('df_eps_marg_train_bootNx_jointPC',n_pc,'_marg_bootstrap_beta',str_beta,'.RData', sep=""))
}

#################################
######       PLOTTING       #####
#################################
# load the eps margi of bootstrapped validation
setwd(metamodel_folder)
load(paste('eps_marg_valid_bootstrap_beta',str_beta,'.RData', sep =''))

quantile_95 <- quantile(eps_marg_valid_bootstrap, probs = 0.95)
quantile_05 <- quantile(eps_marg_valid_bootstrap, probs = 0.05)

quantile_25 <- quantile(eps_marg_valid_bootstrap, probs = 0.25)
quantile_75 <- quantile(eps_marg_valid_bootstrap, probs = 0.75)

print(c(str_beta, quantile_95, quantile_75, quantile_25, quantile_05))

# load the eps margi of bootstrapped validation
load(paste('df_eps_marg_train_PC1_bootstrap_beta',str_beta,'.RData', sep =""))
temp1 <- df_eps_marg_train_bootstrap
load(paste('df_eps_marg_train_PC2_marg_bootstrap_beta',str_beta,'.RData', sep =""))
temp2.marg <- df_eps_marg_train_bootstrap
load(paste('df_eps_marg_train_PC2_bootstrap_beta',str_beta,'.RData', sep =""))
temp2 <- df_eps_marg_train_bootstrap
load(paste('df_eps_marg_train_PC3_bootstrap_beta',str_beta,'.RData', sep =""))
temp3 <- df_eps_marg_train_bootstrap
load(paste('df_eps_marg_train_PC4_bootstrap_beta',str_beta,'.RData', sep =""))
dfdf_eps_marg_train_bootstrap <- rbind(temp1, temp2.marg, temp2, temp3, df_eps_marg_train_bootstrap)
dfdf_eps_marg_train_bootstrap$nPCs <- c(rep('nPC=1 KDE joint', length(temp1$eps_marg)), 
                                        rep('nPC=2 KDE ind.', length(temp2.marg$eps_marg)),
                                        rep('nPC=2 KDE joint', length(temp2$eps_marg)),
                                        rep('nPC=3 KDE joint', length(temp3$eps_marg)),
                                        rep('nPC=4 KDE joint', length(df_eps_marg_train_bootstrap$eps_marg)))
head(dfdf_eps_marg_train_bootstrap)

## PLOT
gg_eps_marg_train_bootstrap <- ggplot(dfdf_eps_marg_train_bootstrap, aes(x=num_traj, y=eps_marg, color = nPCs,  group = interaction(num_traj,nPCs))) + 
  geom_boxplot(outlier.colour=NA)  + theme_bw() + # ggtitle(paste('eps marg beta = ',beta, sep ='')) +
  scale_color_manual(name = "", values = c('nPC=1 KDE joint' = "#253494", 'nPC=2 KDE ind.' = "#225ea8", 'nPC=2 KDE joint' = "#2c7fb8",
                                           'nPC=3 KDE joint' = "#41b6c4", 'nPC=4 KDE joint' = "#a1dab4"))+
  #geom_rect(aes(ymin=min(quantile_25), ymax=max(quantile_75), xmin=-Inf, xmax=Inf), color = 'gray') + 
  annotate("rect", ymin=min(quantile_25), ymax=max(quantile_75), xmin=-Inf, xmax=Inf, fill = 'black', alpha=0.5) + 
  annotate("rect", ymin=min(quantile_05), ymax=max(quantile_95), xmin=-Inf, xmax=Inf, fill = 'black', alpha=0.4) + 
  xlab('$N_{train}$') + ylab('$\\epsilon_{marg}$') +
  geom_hline(yintercept = mean(eps_marg_valid_bootstrap), color = 'black', linetype = 'dashed', linewidth = 1) 
gg_eps_marg_train_bootstrap

## to tikZ
setwd("~/Documents/REDACTION/manuscrit_radisic/6_CALIB_MOIST_PROF/figures/")
#tikz(paste('gg_eps_marg_train_bootstrap_beta',str_beta,'.tex', sep =''), standAlone = TRUE, width=4.0, height=2.0)
tikz(paste('gg_eps_marg_jointvsindep_train_bootstrap_beta',str_beta,'.tex', sep =''), standAlone = TRUE, width=4.0, height=2.0)
gg_eps_marg_train_bootstrap
dev.off()
#tools::texi2dvi(paste('gg_eps_marg_train_bootstrap_beta',str_beta,'.tex', sep =''),pdf=T)
tools::texi2dvi(paste('gg_eps_marg_jointvsindep_train_bootstrap_beta',str_beta,'.tex', sep =''),pdf=T)


