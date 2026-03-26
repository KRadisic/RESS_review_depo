##########################################################
###     READ C_ALPHA FROM TRAINING SET
###     INFER DISTRIBUTION OF C_ALPHA
###     PRODUCE NEW SAMPLES C_ALPHA_NEW
###     READ THE VALIDATION C_ALPHA
###     COMPARE : LATENT SPACES, PCE COEFFS AND MARGINAL DISTIBS

rm(list=ls())

library('ggplot2')
library('GGally')
library('ks')

working_dir <- '~/Documents/RESS_review_depo/'

## Open the raw data c_alpha train and test from the OLS
setwd(working_dir)
beta <- 0.005;
str_beta <- '005'; # 00, 002, 005, 01
ml_coeff_multidx_r_TRAIN <- read.csv(paste("ml_multicoeff_c_alpha_OLS_TRAIN_of_R500_Ntrain50_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb",str_beta,".csv", sep=""), header = FALSE)
ml_coeff_multidx_r_TEST <- read.csv(paste("ml_multicoeff_c_alpha_OLS_TEST_of_R500_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb",str_beta,".csv", sep=""), header = FALSE)
coeff_multidx_r_TRAIN <- data.frame(t(ml_coeff_multidx_r_TRAIN))
coeff_multidx_r_TEST <- data.frame(t(ml_coeff_multidx_r_TEST))

## Check that the coefficient c_alpha TRAIN VS TEST have the same behavior 
coeff_raw <- rbind(coeff_multidx_r_TRAIN, coeff_multidx_r_TEST) 
coeff_raw$type <- c(rep('TRAIN', length(coeff_multidx_r_TRAIN[,1])),rep('TEST', length(coeff_multidx_r_TEST[,1])))
# ggpairs(coeff_raw, aes(color=type))
# if they are different, it shows a sensitivity to the exp design of parameters

## Open the basis A_multidx from ML
A_multidx <- read.csv(paste('old_union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb',str_beta,'.csv',sep=""), header = TRUE)

## Get the degrees of the parameters in the basis
param_order_A_basis <- paste(data.frame(t(A_multidx)))
param_order_A_basis # to be used for c_alpha column names

## Set training trajectories' indices 
train_rain_idx <- c(1:81,83:200)
calpha_train <- coeff_multidx_r_TRAIN[train_rain_idx,]
head(calpha_train)
colnames(calpha_train) <- param_order_A_basis

## Set test trajectories' indices 
test_rain_idx <- c(201:499)
calpha_test <- coeff_multidx_r_TEST[test_rain_idx,]
head(calpha_test)
colnames(calpha_test) <- param_order_A_basis

###########################################
###     PRINCIPAL COMPONENT ANALYSIS
###########################################
## Deduce PCA from train c alpha
pca_train <- prcomp(calpha_train)
n_pc_percentage <- trunc(1000*cumsum((pca_train$sdev)^2/sum((pca_train$sdev)^2)))/10 #  87.2  98.0  98.9  99.5  99.8
n_pc_percentage

## Get the maximal number of PCs
n_pc <- length(A_multidx[,1]) # choose number of principal components
df_pca <- as.data.frame(pca_train$x)

## Quick ploc of PC1, PC2 : round U shape
gg_pc_train <- ggplot(df_pca, aes(x = PC1, y = PC2)) + geom_point(shape=1) +theme_bw()
gg_pc_train

###########################################
##   Project the TEST set onto the PCA basis obtained from the TRAIN set
###########################################
PC_scores_test <- predict(pca_train, newdata = calpha_test)

## Function reprojects PCA scores (PC_scores) to the physical space, through the specified basis (pca)
pca.recon <- function(pca, PC_scores, n_pc){
  mu <- matrix(rep(pca$center, nrow(PC_scores)), nrow = nrow(PC_scores), byrow = T)
  recon <- PC_scores[,1:n_pc] %*% t(pca$rotation[,1:n_pc]) + mu
  return(recon)
}

###########################################
##    PLOTTING to see what is lost in the PCA reduction
print(n_pc)

# COMPARE WHEN ALL NPC ARE KEPT : for a high number of n_pc this should be very similar to coeff_multidx_r_TRAIN
c_alpha_hat_training_reconstructed <- pca.recon(pca_train, pca_train$x, n_pc) 
c_alpha_hat_test_reconstructed <- pca.recon(pca_train, PC_scores_test, n_pc) 

## PLOT c_alpha_hat_training_reconstructed vs x_train, to see what was "lost" in PCA step
df_c_alpha_trainvs_pcareprojected <- rbind(c_alpha_hat_training_reconstructed, 
                                           c_alpha_hat_test_reconstructed, 
                                           calpha_train, calpha_test)
df_c_alpha_trainvs_pcareprojected$type <- c(rep('recon.train', length(c_alpha_hat_training_reconstructed[,1])), 
                                            rep('recon.test', length(c_alpha_hat_test_reconstructed[,1])), 
                                            rep('train', length(calpha_train[,1])), 
                                            rep('test', length(calpha_test[,1])))

# Define a mapping of types to pch values
pch_mapping <- c("recon.train" = 4,"recon.test" = 4, "train" = 1, "test" = 1)  

var_magnitudes <- sapply(abs(df_c_alpha_trainvs_pcareprojected[1:length(A_multidx[,1])]), var)
top_columns <- names(sort(var_magnitudes, decreasing = TRUE)[1:9])
#mean_magnitudes <- colMeans(abs(df_c_alpha_trainvs_pcareprojected[1:length(A_multidx[,1])]), na.rm = TRUE)
#top_columns <- names(sort(mean_magnitudes, decreasing = TRUE)[1:6])

ggpairs(
  df_c_alpha_trainvs_pcareprojected[, c("type", top_columns)], columns = 2:6, aes(color = type, alpha = 0.5),
  upper = list(continuous = wrap("density", alpha = 0.6), combo = "box"),
  lower = list(continuous = function(data, mapping) {
    ggplot(data, mapping) +
      geom_point(aes(pch = type), alpha = 0.8, size = 3) +
      scale_shape_manual(values = pch_mapping)}, 
    combo = wrap("dot", alpha = 0.1, size = 0.2)), legend = 1) +theme_bw()+
  ggtitle(paste('PCE coefficients space for nPC = ', n_pc, ' beta = 0.',str_beta, sep = "")) +
  theme(legend.position = 'bottom') +
  guides(alpha = "none") 
## perfect reconstruction when using all the PCA dimensions

###########################################
###     TAKE A LOWER NUMBER OF N_PC
###########################################
n_pc <- 4 ### 3 or 4 PCA components

# Transform the TRAIN projection coefficients to the original space of PCE coefficients
c_alpha_hat_training_reconstructed <- pca.recon(pca_train, pca_train$x, n_pc) 

# Transform the TEST projection coefficients back to the original space of PCE coefficients
# the PCA basis found on the TRAIN set is used ()pca_train
c_alpha_hat_test_reconstructed <- pca.recon(pca_train, PC_scores_test, n_pc) 

## PLOT c_alpha_hat_training_reconstructed vs x_train, to see what was "lost" in PCA step
df_c_alpha_trainvs_pcareprojected <- rbind(c_alpha_hat_training_reconstructed, 
                                           c_alpha_hat_test_reconstructed, 
                                           calpha_train, calpha_test)
df_c_alpha_trainvs_pcareprojected$type <- c(rep('recon.train', length(c_alpha_hat_training_reconstructed[,1])), 
                                            rep('recon.test', length(c_alpha_hat_test_reconstructed[,1])), 
                                            rep('train', length(calpha_train[,1])), 
                                            rep('test', length(calpha_test[,1])))
ggpairs(
  df_c_alpha_trainvs_pcareprojected[, c("type", top_columns)], columns = 2:6,
        aes(color = type, alpha = 0.5),
  upper = list(continuous = wrap("density", alpha = 0.6), combo = "box"),
  lower = list(continuous = function(data, mapping) {
    ggplot(data, mapping) +
      geom_point(aes(pch = type), alpha = 0.8, size = 3) +
      scale_shape_manual(values = pch_mapping)}, 
    combo = wrap("dot", alpha = 0.1, size = 0.2)), legend = 1) + theme_bw() +
  ggtitle(paste('PCE coefficients space for nPC = ', n_pc, ' beta = 0.', str_beta, sep = "")) +
  theme(legend.position = 'bottom') +
  guides(alpha = "none") 

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
###     Figure 7 : latent space KDE vs TEST
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
c_alpha_KDE_new <- pca.recon(pca_train, as.matrix(df_new_points[,1:n_pc]), n_pc)
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
###     Figure 8 : c_alpha KDE vs TEST
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
#setwd("~/Documents/REDACTION/manuscrit_radisic/5_CALIB_MOIST_PROF/figures/")
#png(paste('gg_figure_calpha_test_beta',str_beta,'.png', sep=""))
#print(gg_figure_calpha_test)
#dev.off()
