##########################################################
###     READ C_ALPHA FROM TRAINING SET
###     INFER DISTRIBUTION OF C_ALPHA
###     PRODUCE NEW SAMPLES C_ALPHA_NEW
###     READ THE VALIDATION C_ALPHA
###     COMPARE : LATENT SPACES, PCE COEFFS AND MARGINAL DISTIBS
###     TO DO : CALCULATE GLOBAL QUALITY OF MARGINAL DISTRIBUTION
###     TO DO : COMPARE QUALITY OF MARGINAL DISTIBUTIONS FOR DIFFERENT ZONES.

rm(list=ls())
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_functions.R")
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")
metamodel_folder <- '~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/scripts_PCE/metamodels/'
figures_plot <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/figures/"
excursion_set_results <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/processed/excursion_sets/"

## Open the c_alpha train and validation from the OLS
setwd(metamodel_folder)
beta <- 0.005;
str_beta <- '005'; # 00, 002, 005, 01
ml_coeff_multidx_r_TRAIN <- read.csv(paste("ml_multicoeff_c_alpha_OLS_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb",str_beta,".csv", sep=""), header = FALSE)
ml_coeff_multidx_r_TEST <- read.csv(paste("ml_multicoeff_c_alpha_OLS_TEST_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb",str_beta,".csv", sep=""), header = FALSE)
coeff_multidx_r_TRAIN <- data.frame(t(ml_coeff_multidx_r_TRAIN)) # 200 realizations of c_alpha
coeff_multidx_r_TEST <- data.frame(t(ml_coeff_multidx_r_TEST)) # 500 realizations of c_alpha
# 10 14 23 43 73 80 82 86 87 123 130 131 151 155 164 173 175 178 188

## Check that the coefficient c_alpha TRAIN VS TEST have the same behavior 
coeff_raw <- rbind(coeff_multidx_r_TRAIN, coeff_multidx_r_TEST) 
coeff_raw$type <- c(rep('TRAIN', length(coeff_multidx_r_TRAIN[,1])),rep('TEST', length(coeff_multidx_r_TEST[,1])))
# ggpairs(coeff_raw, aes(color=type))
# if they are different, it shows a sensitivity to the exp design of parameters
# thus, the Ntrain should probably be augmented

## Open the basis A_multidx from ML
A_multidx <- read.csv(paste('union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb',str_beta,'.csv',sep=""), header = TRUE)

## Get the degrees of the parameters in the basis
param_order_A_basis <- paste(data.frame(t(A_multidx)))
param_order_A_basis # to be used for c_alpha column names

###########################################
###     PRINCIPAL COMPONENT ANALYSIS
###########################################
# Define training data indices 
train_rain_idx <- c(1:81,83:200)
calpha_train <- coeff_multidx_r_TRAIN[train_rain_idx,]
head(calpha_train)
colnames(calpha_train) <- param_order_A_basis

# Define test data indices 
test_rain_idx <- c(201:499)
#test_rain_idx <- c(1:81,83:200)
calpha_test <- coeff_multidx_r_TEST[test_rain_idx,]
head(calpha_test)
colnames(calpha_test) <- param_order_A_basis

# PCA
pca <- prcomp(calpha_train)
#save('pca', file="pca.RData")

## PLOT Scree plot and PC scores
fviz_screeplot(pca, ncp=10)
n_pc_percentage <- trunc(1000*cumsum((pca$sdev)^2/sum((pca$sdev)^2)))/10 # 0.8731214 0.9807023 0.9890268
n_pc_percentage
ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC2)) + geom_point() + theme_bw()
ggplot(as.data.frame(pca$x), aes(x = PC1, y = PC3)) + geom_point() + theme_bw()

## choose the number of PCs to be saved (n_pc) and save in a csv
n_pc <- length(A_multidx[,1]) # choose number of principal components
df_pc <- as.data.frame(pca$x)[,seq(1,n_pc)]
row.names(df_pc) <- mapply(paste,'rain_',train_rain_idx, sep = "")
#write.csv(df_pc, "PC_train_set.csv", row.names = TRUE)
df_pca <- as.data.frame(pca$x)

gg_pc_train <- ggplot(df_pca, aes(x = PC1, y = PC2)) + geom_point(shape=1) +theme_bw()
gg_pc_train
## still round U shape..

###########################################
# Project the test set onto the PCA basis
PC_scores_test <- predict(pca, newdata = calpha_test)

###########################################
#####     SEE WHAT WAS LOST IN THE PCA
## reproject the PCA scores back to the original basis
pca.recon <- function(pca, PC_scores, n_pc){
  mu <- matrix(rep(pca$center, nrow(PC_scores)), nrow = nrow(PC_scores), byrow = T)
  recon <- PC_scores[,1:n_pc] %*% t(pca$rotation[,1:n_pc]) + mu
  return(recon)
}
# COMPARE WHEN ALL NPC ARE KEPT : for a high number of n_pc this should be very similar to coeff_multidx_r_TRAIN
#c_alpha_hat_training_reconstructed <- pca.recon(pca, pca$x, n_pc) 

## PLOT c_alpha_hat_training_reconstructed vs x_train, to see what was "lost" in PCA step
#df_c_alpha_trainvs_pcareprojected <- rbind(c_alpha_hat_training_reconstructed,calpha_train, calpha_test)
#df_c_alpha_trainvs_pcareprojected$type <- c(rep('recon.train', length(c_alpha_hat_training_reconstructed[,1])), 
#                                            rep('train', length(calpha_train[,1])), 
#                                            rep('test', length(calpha_test[,1])))
#mean_magnitudes <- colMeans(abs(df_c_alpha_trainvs_pcareprojected[1:length(A_multidx[,1])]), na.rm = TRUE)
#top_6_columns <- names(sort(mean_magnitudes, decreasing = TRUE)[1:6])
#ggpairs(df_c_alpha_trainvs_pcareprojected[, c("type", top_6_columns)], columns = 2:7, aes(color = type, alpha = 0.5)) +theme_bw()+
#  ggtitle(paste('PCE coefficients space for nPC = ', n_pc, ' beta = 0.',str_beta, sep = "")) +
#  theme(legend.position = 'bottom') +
#  guides(alpha = "none") 

###########################################
###     TAKE A LOWER NUMBER OF N_PC
###########################################
n_pc <- 4 ### 3 or 4 PCA components
c_alpha_hat_training_reconstructed <- pca.recon(pca, pca$x, n_pc) 

## PLOT c_alpha_hat_training_reconstructed vs x_train, to see what was "lost" in PCA step
df_c_alpha_trainvs_pcareprojected <- rbind(c_alpha_hat_training_reconstructed,
                                           calpha_train,
                                           calpha_test)
df_c_alpha_trainvs_pcareprojected$type <- c(rep('recon.train', length(c_alpha_hat_training_reconstructed[,1])), 
                                            rep('train', length(calpha_train[,1])),
                                            rep('test', length(calpha_test[,1])))
# top 6 columns are chosen
mean_magnitudes <- colMeans(abs(df_c_alpha_trainvs_pcareprojected[1:length(A_multidx[,1])]), na.rm = TRUE)
top_6_columns <- names(sort(mean_magnitudes, decreasing = TRUE)[1:9])
ggpairs(df_c_alpha_trainvs_pcareprojected[, c("type", top_6_columns)], columns = 2:10,
        aes(pch = type,color = type, alpha = 0.5),
        # upper = list(continuous = wrap("density", alpha = 0.6), combo = "box"),
        lower = list(continuous = wrap("points", alpha = 0.5, size=1, pch = 1), 
                     combo = wrap("dot", alpha = 0.4, size=0.2) )) +theme_bw() +
  ggtitle(paste('PCE coefficients space for nPC = ', n_pc, ' beta =', beta, sep = "")) +
  theme(legend.position = 'bottom') 

###########################################
###     KERNEL DENSITY ESTIMATOR
###########################################
kde_fitted <- kde(df_pca[,1:n_pc]) #, h= 0.1

# Simulate new points with KDE
new_simu_KDE <- 1000
new_pc_points <- rkde(new_simu_KDE, kde_fitted)
df_new_points <- data.frame(new_pc_points)
head(df_new_points)

## make plot calpha train vs calpha new
df_latent_trainvsKDE <- rbind(df_new_points,df_pca[,1:n_pc], PC_scores_test[,1:n_pc])
colnames(df_latent_trainvsKDE) <- paste('PC', seq(1,4), ' (',n_pc_percentage[1:4],'\\%)',sep='')
df_latent_trainvsKDE$type <- c(rep('KDE', length(df_new_points[,1])), 
                               rep('train', length(df_pca[,1])), 
                               rep('test', length(PC_scores_test[,1])))

###########################################
###     Figure 2 : latent space KDE vs TEST
###########################################
# Ensure 'type' is a factor
df_latent_trainvsKDE$type <- as.factor(df_latent_trainvsKDE$type)

# Define a mapping of types to pch values and colors
pch_mapping <- c("KDE" = 4, "train" = 1, "test" = 1)  
color_mapping <- c("KDE" = "#d95f02", "train" = "#7570b3", "test" = "#1b9e77") 

# Create the ggpairs plot with different pch types and colors
gg_figure_latent_test <- ggpairs(
  df_latent_trainvsKDE, columns = 1:n_pc, aes(color = type, alpha = 0.5),
  upper = list(continuous = wrap("density", alpha = 0.6), combo = "box"),
  lower = list(continuous = function(data, mapping) {
    ggplot(data, mapping) +
      geom_point(aes(pch = type), alpha = 0.5, size = 1) +
      scale_shape_manual(values = pch_mapping)+
      scale_color_manual(values = color_mapping)}, 
  combo = wrap("dot", alpha = 0.1, size = 0.2)), legend = 1) + theme_bw()+
  scale_color_manual('',values = color_mapping) +
  scale_fill_manual('',values = color_mapping) +
  ggtitle(paste('Latent space for nPC = ', n_pc, ' beta = ', beta, sep = "")) +
  theme(legend.position = 'bottom') +
  guides(alpha = "none")  # Add this line to remove alpha from the legend

# Print the plot
print(gg_figure_latent_test)

## to tikZ
#setwd("~/Documents/REDACTION/manuscrit_radisic/5_CALIB_MOIST_PROF/figures/")
#tikz(paste('gg_figure_latent_test_beta',str_beta,'.tex', sep = ''), standAlone = TRUE, width=7.3, height=5.0)
#print(gg_figure_latent_test) # does not add a page alone
#dev.off()
#tools::texi2dvi(paste('gg_figure_latent_test_beta',str_beta,'.tex', sep = ''),pdf=T)

#################################################
###  COMPARE NEW KDE POINTS REPROJECTED ON ORIGINAL BASIS
## these are the new c_alpha from the metamodel
c_alpha_KDE_new <- pca.recon(pca, as.matrix(df_new_points[,1:n_pc]), n_pc)
c_alpha_KDE_new <- data.frame(c_alpha_KDE_new)
colnames(c_alpha_KDE_new) <- colnames(calpha_train)

## Save the new c_alpha KDE
#save(c_alpha_KDE_new, file = paste('c_alpha_KDE_new_beta',str_beta,'.RData', sep = ""))

## make plot calpha train vs calpha new
df_c_alpha_trainvsKDE <- rbind(c_alpha_KDE_new,calpha_train, calpha_test)
df_c_alpha_trainvsKDE$type <- c(rep('$c_{\\alpha}^{new}$', length(c_alpha_KDE_new[,1])), 
                                rep('$c_{\\alpha}^{(train)}$', length(calpha_train[,1])), 
                                rep('$c_{\\alpha}^{(test)}$', length(calpha_test[,1])))

###########################################
###     Figure 3 : c_alpha KDE vs TEST
###########################################
# Step 1: Calculate the mean magnitude for each column
#mean_magnitudes <- colMeans(abs(df_c_alpha_trainvsKDE[1:length(A_multidx[,1])]), na.rm = TRUE)
var_magnitudes <- sapply(abs(df_c_alpha_trainvsKDE[1:length(A_multidx[,1])]), var)
top_columns <- names(sort(var_magnitudes, decreasing = TRUE)[1:9])
df_c_alpha_trainvsKDE$type <- as.factor(df_c_alpha_trainvsKDE$type)

# pretty labels
color_mapping_calpha <- c("$c_{\\alpha}^{new}$" = "#d95f02", "$c_{\\alpha}^{(train)}$" = "#7570b3", "$c_{\\alpha}^{(test)}$" = "#1b9e77") 
pch_mapping_calpha <- c("$c_{\\alpha}^{new}$" = 4, "$c_{\\alpha}^{(train)}$" = 1, "$c_{\\alpha}^{(test)}$" = 1) 


# Create the ggpairs plot with different pch types and colors
gg_figure_calpha_test <- ggpairs(
  df_c_alpha_trainvsKDE[, c("type", top_columns)],
  columns = c(2,3,4,5,6),
  aes(color = type, alpha = 0.5),
  upper = list(continuous = wrap("density", alpha = 0.6), combo = "box"),
  lower = list(continuous = function(data, mapping) {
    ggplot(data, mapping) +
      geom_point(aes(pch = type), alpha = 0.5, size = 1) +
      scale_shape_manual(values = pch_mapping_calpha) +
     scale_color_manual(values = color_mapping_calpha) # scatterplot colors
  }, combo = wrap("dot", alpha = 0.4, size = 0.2)), legend = 1) + theme_bw()+
  scale_color_manual('', values = color_mapping_calpha) + # density plots
  scale_fill_manual('', values = color_mapping_calpha) + # diagonal fill
  ggtitle(paste('Comparison of PCE coefficients for nPC = ', n_pc,' beta = ', beta, sep = "")) +
  theme(legend.position = 'bottom') +
  guides(alpha = "none")  # Add this line to remove alpha from the legend
 
# Print the plot
print(gg_figure_calpha_test)

## to tikZ
#setwd("~/Documents/REDACTION/manuscrit_radisic/5_CALIB_MOIST_PROF/figures/")
#tikz(paste('gg_figure_calpha_test_beta',str_beta,'.tex', sep=""), standAlone = TRUE, width=9.5, height=6.5)
#print(gg_figure_calpha_test) # does not add a page alone
#dev.off()
#tools::texi2dvi(paste('gg_figure_calpha_test_beta',str_beta,'.tex', sep=""),pdf=T)

## to png
setwd("~/Documents/REDACTION/manuscrit_radisic/5_CALIB_MOIST_PROF/figures/")
png(paste('gg_figure_calpha_test_beta',str_beta,'.png', sep=""))
print(gg_figure_calpha_test)
dev.off()

## save the new coefficients C_ALPHA
#head(c_alpha_KDE_new)
#setwd(metamodel_folder)
#write.csv(c_alpha_KDE_new, file="c_alpha_hat_new_LogN02_truerain53_Jb01.csv", row.names = FALSE)

##########################################################
###     COMPARE HISTOGRAMS OF STOCH MM AND TRAIN SET
##########################################################
###     read newly sampled c_alpha values (from the metamodel)
###     then project these PCE coefficients on the physical space
###     compare the histograms generated by the metamodel with 
###     a test set of histograms generated by PESHMELBA
###     the sets set is created with new rains and parameter points
###     in multiple fixed points of the parameter space
##########################################################

#rm(list=ls())
## change to your working directory
working_directory <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/"
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_functions.R")
setwd(paste(working_directory, "/scripts_PCE/metamodels/", sep=""))

## Read the new c_alpha obtained as projections of the GMM to the PCA basis
#c_alpha_KDE_new <- read.csv('c_alpha_hat_new_LogN02_truerain53_Jb01.csv', header = TRUE) # coeffs Fabio (KDE 2 dim)
#head(c_alpha_KDE_new)

# Read the c_alpha training the basis, coming from the OLS previously fitted
#ml_coeff_multidx_r_TRAIN <- read.csv("ml_multicoeff_c_alpha_OLS_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb01.csv", header = FALSE)
#coeff_multidx_r_TRAIN <- data.frame(t(ml_coeff_multidx_r_TRAIN))
#train_rain_idx <- 1:200
#calpha_train <- coeff_multidx_r_TRAIN[train_rain_idx,]
#head(calpha_train)

# open A_multidx from ML
#A_multidx <- read.csv('union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb01.csv',header = TRUE)

# open priors from Emilie marginals
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

## read train and validation set
nR_train <- 200 # number of trajectories used for training the c_alpha
nR_validation <- 299 # number of trajectories used for validation

Y_train <- npyLoad('Ymoisture_profile_LHS_rep5_50_R500_dim6_YzeronS06_LogN02.npy')[1:(50*nR_train),]
Y_validation <- npyLoad('Ymoisture_profile_LHS25_100_R499_dim6_YzeronS06_LogN02.npy')[(1+200*100):(499*100),]

x_train <- npyLoad('LHS_rep5_50_R500_dim6_sample_final.npy')[1:(50*nR_train),] # plutot 50
x_validation <- npyLoad('LHS25_100_R500_dim6_sample_final.npy')[1:(100*nR_validation),]

###############################
###    TRAIN COST FUNCTION
if (length(Y_train[1,]) != length(Y_true))
{stop("Error : dimensions not matching")}
delta_Y <- Y_train - t(matrix(Y_true, length(Y_true),length(Y_train[,1])))
delta_Y2 <- delta_Y^2
#  weights = c(0.5,1,2,3,4,5,6,10,15,20,25,30,35,40,45,50,55,65,75,100,150,200,250,300,400)
weights <- rep(1,length(Y_train[1,]))
J_train <- delta_Y2%*%weights
x_J_train_rains <- data.frame(th9   = x_train[,63], 
                              mn10  = x_train[,68], 
                              th10  = x_train[,72], 
                              th13  = x_train[,99], 
                              thr10 = x_train[,71], 
                              hg10  = x_train[,65], 
                              J0 = J_train) # ML order
# Calculate the Jb term,regularization for mn and thetar
Jb.train <- beta*(((x_J_train_rains$mn10 - df_prior$mean[2])/df_prior$sd[2])^2 +
                    ((x_J_train_rains$thr10 - df_prior$mean[5])/df_prior$sd[5])^2)
x_J_train_rains$J <- x_J_train_rains$J0 + Jb.train

#################################################################
KS_pvalues <- c() # here we will keep the results of the hypothesis testing
data1_all_test_histograms <- c() # here we will keep the results of the hypothesis testing
data2_all_metamodel_histograms <- c() # here we will keep the results of the hypothesis testing

## check nb_x_param
nb_x_param <- 50 # number of points in the parameter space in which histograms are observed
dim_x_param <- 6 # dimension of the parameter space
vect_params <- 1:dim_x_param
idx_test_x_param <-6

###############################
###    TRAIN COST FUNCTION
# loop over the 100 different points in the parameter space, one at a time
for (idx_test_x_param in seq(1,nb_x_param)){
  ## fix the parameter value to a point from the test set rain dataset from 1-nb_x_param (1 is the "unknown" true parameter)
  x_param_fixed <- x_J_train_rains[idx_test_x_param, vect_params]
  # read the corresponding cost function realizations
  df_test_J_histogram_in_fixed_point <- x_J_train_rains[idx_test_x_param + nb_x_param*seq(0,nR_train) , ] # choose from 1-100 (1 is the true parameter though)
  # verify it is the same point in the parameter space repeated as many time as R_train + R_test
  if (table(table(df_test_J_histogram_in_fixed_point[vect_params]))[[1]] != 1)
  {print("error : the test set is not evaluated in a fix point of the parameter space")}
  # keep only the realizations of the test set cost function
  test_J_histogram_in_fixed_point <- df_test_J_histogram_in_fixed_point$J[train_rain_idx]
  # in the fixed point x_param_fixed evaluate all PCE as with all c_alpha
  new_J_histogram_in_fixed_point <- eval_surrogate_OLS_in_point_x_variable_rain(x_param_fixed, df_prior, A_multidx, 
                                                                                t(c_alpha_KDE_new)) 
  # about 0.01 sec per simulation, 12 sec per 1000 simu, 1200 sec = 60min = 1 hour  per 1000*100 simu
  ## Compare the two sets of data (two histograms) with a KS test.
  data1 <- test_J_histogram_in_fixed_point # the test histogram in a fixed test point x
  data2 <- new_J_histogram_in_fixed_point # the metamodel histogram in this same point
  kstest <- ks.test(data1, data2) # p-value small : probability of coming from the same distribution is small
  
  data1_all_test_histograms <- cbind(data1_all_test_histograms, data1)
  data2_all_metamodel_histograms  <- cbind(data2_all_metamodel_histograms, data2)
  
  KS_pvalues <- c(KS_pvalues, kstest$p.value)
  print(kstest$p.value)
  
  ## Plot and save THE TWO HISTOGRAMS and the TWO SAMPLE CUMULATIVE DENSITY FUNCTIONS FOR COMPARISON
  df <- data.frame(value_J = c(test_J_histogram_in_fixed_point, new_J_histogram_in_fixed_point), 
                   type = c(rep('real model with test rains', length(data1)), rep('stochastic MM', length(data2))))
  
  gg_hist <- ggplot(df, aes(x = value_J, y =..density.., fill = type, color = type)) + 
    geom_histogram(alpha = 0.2, position="identity", binwidth = 0.001) + theme_bw()  + #xlim(c(0,0.35)) +  
    ggtitle(paste('Comparison of histograms in test_point_idx = ',idx_test_x_param,'; KS p-val = ', trunc(100*kstest$p.value)/100,sep=''))
  gg_hist
  
  #gg_cdf <- ggplot(df, aes(x = value_J, fill = type, color = type)) + 
  #  stat_ecdf(geom = "step") + xlim(c(0,0.35)) +  theme_bw() + 
  #  ggtitle(paste('Comparison of CDF in test_point_idx = ',idx_test_x_param,
  #                '; KS p-val = ', trunc(100*kstest$p.value)/100,sep=''))
  ## save figures
  setwd(paste(working_directory, "/output_figures", sep=""))
  png(paste('histogram_test_rains_test_parameter_LogN02_idx', idx_test_x_param, '.png', sep=""))
  print(gg_hist)
  dev.off()
  #png(paste('cdf_test_rains_test_parameter_idx', idx_test_x_param, '.png', sep=""))
  #print(gg_cdf)
  #dev.off()
}
###############################
###    TRAIN COST FUNCTION
#  null hypothesis : x and y were drawn from the same distribution is performed. 
sum(KS_pvalues>0.05) # YES these histograms come from the same distributions 
sum(KS_pvalues<0.05) # NO these histograms do not come from the same distributions


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

#################################################################
KS_pvalues <- c() # here we will keep the results of the hypothesis testing
data1_all_test_histograms <- c() # here we will keep the results of the hypothesis testing
data2_all_metamodel_histograms <- c() # here we will keep the results of the hypothesis testing

## check nb_x_param
nb_x_param <- 100 # number of points in the parameter space in which histograms are observed
dim_x_param <- 6 # dimension of the parameter space
vect_params <- 1:dim_x_param
idx_test_x_param <-6

###############################
###    VALIDATION CONST FUNCTION
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
  
  ## Compare the two sets of data (two histograms) with a KS test.
  kstest <- ks.test(data1, data2) # p-value small : probability of coming from the same distribution is small
  
  data1_all_test_histograms <- cbind(data1_all_test_histograms, data1)
  data2_all_metamodel_histograms  <- cbind(data2_all_metamodel_histograms, data2)
  
  KS_pvalues <- c(KS_pvalues, kstest$p.value)
  print(paste('idx_x = ', idx_test_x_param, 'p.val', kstest$p.value))
  
  ## Plot and save THE TWO HISTOGRAMS and the TWO SAMPLE CUMULATIVE DENSITY FUNCTIONS FOR COMPARISON
  df <- data.frame(value_J = c(test_J_histogram_in_fixed_point, new_J_histogram_in_fixed_point), 
                   type = c(rep('real model with test rains', length(data1)), rep('stochastic MM', length(data2))))
  
  gg_hist <- ggplot(df, aes(x = value_J, y =..density.., fill = type, color = type)) + 
    geom_histogram(alpha = 0.2, position="identity", binwidth = 0.001) + theme_bw()  + xlim(c(0,0.15)) +  
    ggtitle(paste('Comparison of histograms in test_point_idx = ',idx_test_x_param,'; KS p-val = ', trunc(100*kstest$p.value)/100,sep=''))
  gg_hist
  
  #gg_cdf <- ggplot(df, aes(x = value_J, fill = type, color = type)) + 
  #  stat_ecdf(geom = "step") + xlim(c(0,0.35)) +  theme_bw() + 
  #  ggtitle(paste('Comparison of CDF in test_point_idx = ',idx_test_x_param,
  #                '; KS p-val = ', trunc(100*kstest$p.value)/100,sep=''))
  ## save figures
  setwd(paste(working_directory, "/output_figures", sep=""))
  png(paste('histogram_test_rains_test_parameter_LogN02_idx', idx_test_x_param,'_beta',str_beta,'.png', sep=""))
  print(gg_hist)
  dev.off()
  #png(paste('cdf_test_rains_test_parameter_idx', idx_test_x_param, '.png', sep=""))
  #print(gg_cdf)
  #dev.off()
}
#  null hypothesis : x and y were drawn from the same distribution is performed. 
sum(KS_pvalues>0.05) # YES these histograms come from the same distributions 
sum(KS_pvalues<0.05) # NO these histograms do not come from the same distributions

######################################################################
###     MAKE DATA FRAME WITH MARGINAL HISTOGRAMS IN VALIDATION POINTS
######################################################################
df_marginals <- data.frame(rbind(data1_all_test_histograms, data2_all_metamodel_histograms))
colnames(df_marginals) <- paste("$i=", seq(1,length(df_marginals[1,])),"$",sep = "")
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

eps_marg
beta
# > eps_marg 20.51922 beta 0
# > eps_marg 25.21369 beta 0.002
# > eps_marg 20.52382 beta 0.005
# > eps_marg 25.07291 beta 0.01
# save eps_marg wrt beta
setwd(metamodel_folder)
save(eps_marg, file = paste("eps_marg_beta", beta, ".RData"))

######################################################################
###     get eps marg for different domains...?
######################################################################

###     get eps marg for bootstrap on the validation set.


##########################################################################
###     Figure 4 : Marginal distributions in chosen validation points
##########################################################################
# pretty labels
df_marginals$type.pretty <- NaN
df_marginals[df_marginals$type == 'test',]$type.pretty <- '$f_s(\\mathbf{x}_{test}^{(i)}, \\cdot)$'
df_marginals[df_marginals$type == 'KDE',]$type.pretty <- '$\\hat{f}_s(\\mathbf{x}_{test}^{(i)}, \\cdot)$'

pch_mapping_Y <- c('$\\hat{f}_s(\\mathbf{x}_{test}^{(i)}, \\cdot)$' = 4, 
                   '$f_s(\\mathbf{x}_{test}^{(i)}, \\cdot)$' = 1)
color_mapping_Y <- c('$\\hat{f}_s(\\mathbf{x}_{test}^{(i)}, \\cdot)$' = '#d95f02', 
                   '$f_s(\\mathbf{x}_{test}^{(i)}, \\cdot)$' = '#1b9e77')

# Create the ggpairs plot with different pch types and colors
#idxs_chosen_valid_points <- c(38,39,40,41,42,43,44,45) #
idxs_chosen_valid_points <- c(10,20,30,40,50) #
gg_figure_marginals <- ggpairs(
  df_marginals, columns = idxs_chosen_valid_points,
  aes(color = type.pretty, alpha = 0.3),
  upper = list(continuous = wrap("density", alpha = 0.5), combo = "box"),
  lower = list(continuous = function(data, mapping) {
    ggplot(data, mapping) +
      geom_point(aes(pch = type.pretty), alpha = 0.3, size = 1) +
      scale_shape_manual(values = pch_mapping_Y) +
    scale_color_manual('',values = color_mapping_Y)
  }, combo = wrap("dot", alpha = 0.6, size = 0.2)), legend = 1) + theme_bw()+
  scale_color_manual('',values = color_mapping_Y)+
  scale_fill_manual('',values = color_mapping_Y) +
  ggtitle(paste('Comparison of marginals for nPC = ', n_pc, ', beta = ',beta,sep = "")) +
  theme(legend.position = 'bottom') +
  guides(alpha = "none")  # Add this line to remove alpha from the legend

# Print the plot
print(gg_figure_marginals)

## to tikZ
setwd("~/Documents/REDACTION/manuscrit_radisic/5_CALIB_MOIST_PROF/figures/")
tikz('gg_figure_marginals.tex', standAlone = TRUE, width=7.3, height=5.0)
print(gg_figure_marginals) # does not add a page alone
dev.off()
tools::texi2dvi('gg_figure_marginals.tex',pdf=T)
### Marginal distributions WRT trajectories number # in : 10_plot_eps_marg_train_boot_PC1

###################################################################
####   EVALUATE EMULATOR IN NEW POINTS TO GET ROBUST MINIMIZERS
###################################################################
# Isoprobabilistic transform
IsoProb_z_to_x <- function(z, df_prior){
  z*df_prior$sd[1:5] + df_prior$mean[1:5]
}
# A unique variable to be compatible with optimr package
surrogate_cost_function_analytical_fixedrain_1input <- function(x)
{eval_surrogate_OLS_in_point_x(c(x, -3.50), df_prior, A_multidx, coeff_PCE_fixed_rain)}

## fix rain 
n_params <- 5
c_alpha_KDE_new_multidx_r <- t(c_alpha_KDE_new)

###################################################################
### CONDITIONAL MINIMIZERS for each new trajectory (if regred-based robust parameters.)
###################################################################

all_results_minimize_fixed_rain <- c()
for (rain_idx in seq(1,length(c_alpha_KDE_new_multidx_r[1,])))
{
  print(rain_idx)
  coeff_PCE_fixed_rain <- c_alpha_KDE_new_multidx_r[,rain_idx] # to get the coeff_PCE_fixed_rain
  results_minimize_fixed_rain <- optim(df_prior$mean[1:5], 
                                       surrogate_cost_function_analytical_fixedrain_1input, gr = NULL,
                                       method ="L-BFGS-B",
                                       lower = mapply(qnorm, rep(0.001, n_params), mean = df_prior$mean[1:5], sd = df_prior$sd[1:5]), 
                                       upper = mapply(qnorm, rep(0.999, n_params), mean = df_prior$mean[1:5], sd = df_prior$sd[1:5]),
                                       control = list(), hessian = FALSE)
  all_results_minimize_fixed_rain <- rbind(all_results_minimize_fixed_rain, c(results_minimize_fixed_rain$par, results_minimize_fixed_rain$value))
}
df_all_results_minimize_fixed_rain <- data.frame(all_results_minimize_fixed_rain)
names(df_all_results_minimize_fixed_rain) <- c(df_prior$param[1:5], 'J_min')
head(df_all_results_minimize_fixed_rain)

## if negative, put to zero
df_all_results_minimize_fixed_rain$J_min_positive <- ifelse(df_all_results_minimize_fixed_rain$J_min < 0,
                                                            0, df_all_results_minimize_fixed_rain$J_min)

#setwd('~/code/4_PCE_GMM_relative_regret_min/outputs/')
#df <- read.csv(paste('cond_min_OLS_J_GMM1000_Ntrain50_dim5_NEWTON_without_hg_moistureprofile_truerain1_trueparam1_YzeronS06.csv',sep=""), row.names = 1)
## add a column for relative regret

## Choose relative regred quantities
relreg_qty1 <- 0.1
relreg_qty2 <- 0.2
relreg_qty3 <- 0.5

df_all_results_minimize_fixed_rain$J_relative1 <- df_all_results_minimize_fixed_rain$J_min_positive * (1 + relreg_qty1)
df_all_results_minimize_fixed_rain$J_relative2 <- df_all_results_minimize_fixed_rain$J_min_positive * (1 + relreg_qty2)
df_all_results_minimize_fixed_rain$J_relative3 <- df_all_results_minimize_fixed_rain$J_min_positive * (1 + relreg_qty3)

setwd(metamodel_folder)
#save(df_all_results_minimize_fixed_rain, file = "df_all_results_minimize_fixed_rain_beta005_PC4KDE.RData")
load("df_all_results_minimize_fixed_rain_beta005_PC4KDE.RData")

## add another relative regret level
relreg_qty4 <- 2.0
df_all_results_minimize_fixed_rain$J_relative4 <- df_all_results_minimize_fixed_rain$J_min_positive * (1 + relreg_qty4)

head(df_all_results_minimize_fixed_rain)


##############################################
####           EVALUATE EMULATOR          ####
##############################################
excset_val1 <- 0.01
excset_val2 <- 0.02
excset_val3 <- 0.03

## Initialize sets
admissible_points_excset1 <- c(); admissible_points_excset2 <- c(); admissible_points_excset3 <- c()
admissible_points_relreg1 <- c(); admissible_points_relreg2 <- c(); admissible_points_relreg3 <- c()
admissible_points_relreg4 <- c()

## Generate a space-filling design in the parameter space (the MM is evaluated)
n_lhs <- 1000 #10000 # number of sample points
## Create LHS quantiles
z_LHS_quantiles <- maximinLHS(n_lhs, n_params, dup = 1)
# Transform to N(0,1)
z_LHS <- qnorm(z_LHS_quantiles, mean = 0, sd = 1)#
hist(z_LHS[,1])
#z_LHS1 <- cbind(rnorm(n_lhs), rnorm(n_lhs), rnorm(n_lhs), rnorm(n_lhs), rnorm(n_lhs))

## Get different Robust estimators with accept/reject algorithm
# Loop on new trajectories
for (rain_idx in seq(1:length(c_alpha_KDE_new_multidx_r[1,])))
{
  print(paste('omega idx = ', rain_idx))
  # Read pce coefficients corresponding to this trajectory
  c_alpha <- c_alpha_KDE_new_multidx_r[,rain_idx]
  
  # Loop on the space-filling design in the parameter space
  for (idx_theta in seq(n_lhs)){
    # Fix one point from the unitary LHS
    zk <- z_LHS[idx_theta,]
    # Take its value in the parameter space
    theta_k <- c(IsoProb_z_to_x(zk,df_prior = df_prior),-3.50)
    # Evaluate the trajectory in this point
    value_J <- eval_surrogate_OLS_in_point_x(theta_k, df_prior, A_multidx, c_alpha_KDE_new_multidx_r[,rain_idx])
    
    # Accept the point if it fulfills the criteria: 
    if (value_J < df_all_results_minimize_fixed_rain$J_relative1[rain_idx]){
      admissible_points_relreg1 <- rbind(admissible_points_relreg1, c(theta_k, rain_idx, idx_theta))}
    if (value_J < df_all_results_minimize_fixed_rain$J_relative2[rain_idx]){
      admissible_points_relreg2 <- rbind(admissible_points_relreg2, c(theta_k, rain_idx, idx_theta))}
    if (value_J < df_all_results_minimize_fixed_rain$J_relative3[rain_idx]){
      admissible_points_relreg3 <- rbind(admissible_points_relreg3, c(theta_k, rain_idx, idx_theta))}    
    if (value_J < df_all_results_minimize_fixed_rain$J_relative3[rain_idx]){
      admissible_points_relreg4 <- rbind(admissible_points_relreg4, c(theta_k, rain_idx, idx_theta))}
    
    if (value_J < excset_val1){
      admissible_points_excset1 <- rbind(admissible_points_excset1, c(theta_k, rain_idx, idx_theta))}
    if (value_J < excset_val2){
      admissible_points_excset2 <- rbind(admissible_points_excset2, c(theta_k, rain_idx, idx_theta))}
    if (value_J < excset_val3){
      admissible_points_excset3 <- rbind(admissible_points_excset3, c(theta_k, rain_idx, idx_theta))}
  }
}
# prevu :2h00 pour 1000x1000 simulations; commencé à 14h17; fini a 16h17 (2h); 11h55-14h00
# prevu 15h pour 10000x1000 simulations; début 11h36; fin 2AM

################################
###      SAVE RESULTS
################################

# format columns names
#df_admissible_points_excset1 <- data.frame(admissible_points_excset1)
#df_admissible_points_excset2 <- data.frame(admissible_points_excset2)
#df_admissible_points_excset3 <- data.frame(admissible_points_excset3)

#colnames(df_admissible_points_excset1) <- c(df_prior$param, 'rain_idx', 'idx_theta')
#colnames(df_admissible_points_excset2) <- c(df_prior$param, 'rain_idx', 'idx_theta')
#colnames(df_admissible_points_excset3) <- c(df_prior$param, 'rain_idx', 'idx_theta')

#head(df_admissible_points_excset1)
#head(df_admissible_points_excset2)
#head(df_admissible_points_excset3)

# format columns names
#df_admissible_points_relreg1 <- data.frame(admissible_points_relreg1)
#df_admissible_points_relreg2 <- data.frame(admissible_points_relreg2)
#df_admissible_points_relreg3 <- data.frame(admissible_points_relreg3)
#df_admissible_points_relreg4 <- data.frame(admissible_points_relreg4)

#colnames(df_admissible_points_relreg1) <- c(df_prior$param, 'rain_idx', 'idx_theta')
#colnames(df_admissible_points_relreg2) <- c(df_prior$param, 'rain_idx', 'idx_theta')
#colnames(df_admissible_points_relreg3) <- c(df_prior$param, 'rain_idx', 'idx_theta')
#colnames(df_admissible_points_relreg4) <- c(df_prior$param, 'rain_idx', 'idx_theta')

#head(df_admissible_points_relreg1)
#head(df_admissible_points_relreg2)
#head(df_admissible_points_relreg3)
#head(df_admissible_points_relreg4)


# save results
setwd(excursion_set_results)
#save(df_admissible_points_excset1, file = paste("df_admissible_points_excset",excset_val1,"_LHS.RData", sep = ""))
#save(df_admissible_points_excset2, file = paste("df_admissible_points_excset",excset_val2,"_LHS.RData", sep = ""))
#save(df_admissible_points_excset3, file = paste("df_admissible_points_excset",excset_val3,"_LHS.RData", sep = ""))

#save(df_admissible_points_relreg1, file = paste("df_admissible_points_relreg",relreg_qty1,"_LHS.RData", sep = ""))
#save(df_admissible_points_relreg2, file = paste("df_admissible_points_relreg",relreg_qty2,"_LHS.RData", sep = ""))
#save(df_admissible_points_relreg3, file = paste("df_admissible_points_relreg",relreg_qty3,"_LHS.RData", sep = ""))
#save(df_admissible_points_relreg4, file = paste("df_admissible_points_relreg",relreg_qty4,"_LHS.RData", sep = ""))

###########################################
## Load admissible sets for each trajectory realization
excursion_set_results <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/processed/excursion_sets/"
setwd(excursion_set_results)
load(paste("df_admissible_points_excset",excset_val1,"_LHS.RData", sep = ""))
load(paste("df_admissible_points_excset",excset_val2,"_LHS.RData", sep = ""))

load(paste("df_admissible_points_relreg",relreg_qty1,"_LHS.RData", sep = ""))
load(paste("df_admissible_points_relreg",relreg_qty2,"_LHS.RData", sep = ""))
load(paste("df_admissible_points_relreg",relreg_qty3,"_LHS.RData", sep = ""))
load(paste("df_admissible_points_relreg",relreg_qty4,"_LHS.RData", sep = ""))
df_admissible_points_excset0.02Nomega500_nLHS_eachLHS5000

########################
###   EXCURSION SETS
########################
## Count how many times each theta_idx appears in the excursion sets
value_counts1 <- table(df_admissible_points_excset1$idx_theta)
value_counts2 <- table(df_admissible_points_excset2$idx_theta)

## Which theta value is the most frequent
most_frequent_value1 <- names(value_counts1)[which.max(value_counts1)]
most_frequent_value2 <- names(value_counts2)[which.max(value_counts2)]

theta_k1 <- df_admissible_points_excset1[which(df_admissible_points_excset1$idx_theta == most_frequent_value1)[1],]
theta_k2 <- df_admissible_points_excset2[which(df_admissible_points_excset2$idx_theta == most_frequent_value2)[1],]

# Print the most frequent value
theta_k1
theta_k2

## Compare the top 5 most frequent values
top5_most_frequent_value1 <- names(sort(value_counts1, decreasing = TRUE)[1:5])
top5_most_frequent_value2 <- names(sort(value_counts2, decreasing = TRUE)[1:5])

## see how distant these values are one from another
df_admissible_points_excset1[which(df_admissible_points_excset1$idx_theta == top5_most_frequent_value1[1])[1],1:6]
df_admissible_points_excset1[which(df_admissible_points_excset1$idx_theta == top5_most_frequent_value1[2])[1],1:6]
df_admissible_points_excset1[which(df_admissible_points_excset1$idx_theta == top5_most_frequent_value1[3])[1],1:6]

df_admissible_points_excset2[which(df_admissible_points_excset2$idx_theta == top5_most_frequent_value2[1])[1],1:6]
df_admissible_points_excset2[which(df_admissible_points_excset2$idx_theta == top5_most_frequent_value2[2])[1],1:6]
df_admissible_points_excset2[which(df_admissible_points_excset2$idx_theta == top5_most_frequent_value2[3])[1],1:6]

################################
###    RELATIVE REGRET ESTIMATOR
################################
value_counts_relreg1 <- table(df_admissible_points_relreg1$idx_theta)
value_counts_relreg2 <- table(df_admissible_points_relreg2$idx_theta)
value_counts_relreg3 <- table(df_admissible_points_relreg3$idx_theta)
value_counts_relreg4 <- table(df_admissible_points_relreg4$idx_theta)

most_frequent_value_relreg1 <- names(value_counts_relreg1)[which.max(value_counts_relreg1)]
most_frequent_value_relreg2 <- names(value_counts_relreg2)[which.max(value_counts_relreg2)]
most_frequent_value_relreg3 <- names(value_counts_relreg3)[which.max(value_counts_relreg3)]
most_frequent_value_relreg4 <- names(value_counts_relreg4)[which.max(value_counts_relreg4)]

theta_k1_relreg <- df_admissible_points_relreg1[which(df_admissible_points_relreg1$idx_theta == most_frequent_value_relreg1)[1],1:6]
theta_k2_relreg <- df_admissible_points_relreg2[which(df_admissible_points_relreg2$idx_theta == most_frequent_value_relreg2)[1],1:6]
theta_k3_relreg <- df_admissible_points_relreg3[which(df_admissible_points_relreg3$idx_theta == most_frequent_value_relreg3)[1],1:6]
theta_k4_relreg <- df_admissible_points_relreg4[which(df_admissible_points_relreg4$idx_theta == most_frequent_value_relreg4)[1],1:6]

## Compare the top 5 most frequent values
top5_most_frequent_value_relreg1 <- names(sort(value_counts_relreg1, decreasing = TRUE)[1:5])
top5_most_frequent_value_relreg2 <- names(sort(value_counts_relreg2, decreasing = TRUE)[1:5])
top5_most_frequent_value_relreg3 <- names(sort(value_counts_relreg3, decreasing = TRUE)[1:5])
top5_most_frequent_value_relreg4 <- names(sort(value_counts_relreg4, decreasing = TRUE)[1:5])

# Print the most frequent value
rbind(theta_k1_relreg,theta_k2_relreg,theta_k3_relreg,theta_k4_relreg)
rbind(top5_most_frequent_value_relreg1, top5_most_frequent_value_relreg2, 
      top5_most_frequent_value_relreg3, top5_most_frequent_value_relreg4)

## see how distant these values are one from another
df_admissible_points_relreg1[which(df_admissible_points_relreg1$idx_theta == top5_most_frequent_value1[1])[1],1:6]
df_admissible_points_relreg1[which(df_admissible_points_relreg1$idx_theta == top5_most_frequent_value1[2])[1],1:6]
df_admissible_points_relreg1[which(df_admissible_points_relreg1$idx_theta == top5_most_frequent_value1[3])[1],1:6]


##########################################################
#####    PLOT CONDITIONAL MINIMIZERS 
setwd("~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/scripts_PCE/metamodels/")
df_min_LARS <- read.csv('conditional_minimizers_LARS_errorLogN02_Jb_beta005LHS_rep5_50_R500_Ntrain_hourly_BDDquantiles50.csv', header = TRUE)[1:200,]
#df_min_LARS <- read.csv('conditional_minimizers_OLS_errorLogN02_Jb_beta0.005_LHS_rep5_50_R500_Ntrain_hourly_BDDquantiles50.csv', header = TRUE) 

# results of the conditional one shot minimization
head(df_min_LARS)
colnames(df_min_LARS) <- c("$\\theta_{s.inter.}$", "$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep.}$", " $hg_{.deep}$", "f")

## Format conditional minimizers in dataframe for ggplot
R <- 200; N_train <- 50; n_params <- 6; true_test_idx <- 1;

df_cond_minimizers_LARS <- data.frame(cond.min = as.vector(as.matrix(df_min_LARS)[seq(n_params,R),]),
                                      parameter = rep(colnames(df_min_LARS), each = length(seq(n_params,R))))

df_cond_minimizers_LARS$type <- '$\\hat{x}_{PCE}^{*(r)}$'
df_cond_min <- df_cond_minimizers_LARS

######  ROBUST AND TRUE RESULTS AND FORMATTING FOR PLOT : VERTICAL LINES

## Load previous results from other robust calibrations ...
## Load true
setwd('~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/raw/')
x_test <- npyLoad("LHS25_100_R198_dim6_sample_final.npy")
x_true <- x_test[1,c(63,68,72,99,71,65)] # matlab order
#x_pce_sum <- c(0.317852,	0.120000,	0.3416294,	0.2812774,	0.03000000,	-2.000000) # matlab order

## Load prior add true in dataframe for plotting
n_params_plot <- 5
df_prior_true_mean <- df_prior[1:n_params_plot,]
df_prior_true_mean$param <- c("$\\theta_{s.inter.}$", "$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep.}$")
df_prior_true_mean$true_parameter <- x_true[1:n_params_plot]
#df_prior_true_mean$true_rain <- c(as.matrix(df_min_LARS)[1,1:n_params_plot])
df_prior_true_mean$exc_set1 <- as.numeric(theta_k1[1:n_params_plot])
df_prior_true_mean$exc_set2 <- as.numeric(theta_k2[1:n_params_plot])
df_prior_true_mean$relreg1 <- as.numeric(theta_k1_relreg[1:n_params_plot])
df_prior_true_mean$relreg2 <- as.numeric(theta_k2_relreg[1:n_params_plot])
df_prior_true_mean$relreg3 <- as.numeric(theta_k3_relreg[1:n_params_plot])
df_prior_true_mean$relreg4 <- as.numeric(theta_k4_relreg[1:n_params_plot])
#df_prior_true_mean$sum <- x_pce_sum[1:n_params_plot]
head(df_prior_true_mean)

###################################################################
###   PLOT PARAMETER VALUES FROM MINIMIZATION RESULTS
###################################################################

list_of_plots_excset <- list(); z <- 1
for (parameter in df_prior_true_mean$param){
  df_cond_min_param <- df_cond_min[df_cond_min$parameter==parameter,]
  prior_distrib_param <- cbind(df_prior_true_mean[df_prior_true_mean$param == parameter,])
  
  conditional_minimizers <- ggplot(df_cond_min_param, aes(x=cond.min, fill=type)) +
    geom_histogram(position = "identity", aes(y = ..density..),alpha = 0.50, bins = 50, color="black")+
    xlim(qnorm(0.0001, mean = prior_distrib_param$mean, sd = prior_distrib_param$sd),
         qnorm(0.9999, mean = prior_distrib_param$mean, sd = prior_distrib_param$sd)) +
    scale_fill_manual(name = "Conditional minimizers",values=c("skyblue")) + xlab('') + ylab('') +
    facet_grid(. ~ parameter, scales = "free_x" ) + theme_bw() + theme(legend.position="none") +
    stat_function(fun = dnorm, args = list(mean = prior_distrib_param$mean, sd = prior_distrib_param$sd), color ="blue") +
    geom_vline(data = prior_distrib_param,aes(xintercept = true_parameter, color="$x_{true}$"),linetype="solid",linewidth = 1.5)+
    #geom_vline(data = prior_distrib_param,aes(xintercept = true_rain, color="x.pce.z.true"),linetype="dashed",linewidth = 1.5)+
    geom_vline(data = prior_distrib_param,aes(xintercept = exc_set1, color=paste("exc.set f \\textless ", excset_val1 , sep = "")),linetype="dashed",linewidth = 1.0)+
    geom_vline(data = prior_distrib_param,aes(xintercept = exc_set2, color=paste("exc.set f \\textless ", excset_val2 , sep = "")),linetype="dashed",linewidth = 1.0)#+
  #scale_color_manual(name = "Robust estimator", values = c(Conditional = "#c4c0b3", x.true = "red",
  #                                                         exc.set1 = "yellow", exc.set2 = "orange", exc.set2 = "", exc.set3 = "orange",
  #                                                         relreg2 = "violet", relreg2 = "forestgreen"))
  list_of_plots_excset[[z]] <- conditional_minimizers
  z <- z+1
}

## add legend as sixth plot # Using the cowplot package
conditional_minimizers <- conditional_minimizers + theme(legend.position="left")
legend <- get_legend(conditional_minimizers)
list_of_plots_excset[[z]] <- legend
conditional_minimizers_over_prior_true <- do.call(grid.arrange, c(list_of_plots_excset, ncol = 3))
conditional_minimizers_over_prior_true

## to tikZ
setwd("~/Documents/REDACTION/manuscrit_radisic/5_CALIB_MOIST_PROF/figures/")
tikz(paste('gg_arranged_cond_min_robust_estim_excset_PCE_LogN02_Jb',str_beta,'_YzeronS06.tex', sep=''), standAlone = TRUE, width=5.5, height=4.0)
print(grid::grid.draw(conditional_minimizers_over_prior_true)) # does not add a page alone
dev.off()
tools::texi2dvi(paste('gg_arranged_cond_min_robust_estim_excset_PCE_LogN02_Jb',str_beta,'_YzeronS06.tex', sep=''),pdf=T)

###################################################################
###    SAVE MINIMIZATION RESULTS FOR NEW PESHMELBA SIMULATIONS
###################################################################
## Arrange minimizers for new PESHMELBA simulations
r.new.rains <- 300 # number of rains for same parameter values
n.dim <- 6 # number of parameter values per rain

## read previous LHS (for other parameter values)
raw_data_directory <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/raw/"
setwd(raw_data_directory)
old_LHS <- npyLoad("LHS_rep5_50_R500_dim6_sample_final.npy")
single_LHS <- old_LHS[1:n.dim, 1:145]

# random number from 1 to 200 for the conditional minimizer
idx_of_conditional_minimizer <- 174 #c(15,82,113,122) 
x_bar_PCE <- colMeans(df_min_LARS)

## Write the new values for the exp design for new peshmelba simulations
single_LHS[1, vector_input_indices] <- as.numeric(df_min_LARS[idx_of_conditional_minimizer, 1:6])
single_LHS[2, vector_input_indices] <- as.numeric(c(x_bar_PCE[1:5], -3.50))
single_LHS[3, vector_input_indices] <- as.numeric(c(df_prior_true_mean$relreg1, -3.50))
single_LHS[4, vector_input_indices] <- as.numeric(c(df_prior_true_mean$relreg3, -3.50))
single_LHS[5, vector_input_indices] <- as.numeric(c(df_prior_true_mean$exc_set1, -3.50))
single_LHS[6, vector_input_indices] <- as.numeric(c(df_prior_true_mean$exc_set2, -3.50))

head(single_LHS[,vector_input_indices])
head(single_LHS[1:6,vector_input_indices]) # good parameters in good places

## Repeat the same exp design as many times as there are new rains
training_set_r <- c()
for (idx_r in 1:r.new.rains)
{
  training_set_r <- rbind(training_set_r, single_LHS)
}

training_set_r[,72]

## Save new LHS
processed_data_directory <- '~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/processed/'
setwd(processed_data_directory)
npySave(paste("LHS_robust_minimizers_LogN02_beta005_condmin",idx_of_conditional_minimizer,"_dim6_sample_final.npy", sep=""), training_set_r)



