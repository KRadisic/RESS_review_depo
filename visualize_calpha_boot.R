##########################################################
rm(list=ls())

library('ggplot2')
library('GGally')

working_dir <- '~/Documents/RESS_review_depo/'

## Open the raw data c_alpha from the OLS
setwd(working_dir)
beta <- 0.005;
str_beta <- '005'; 
Nx_boot <- 100
ml_coeff_multidx_r_TEST_boot <- read.csv(paste("ml_boot_multicoeff_c_alpha_OLS_TESTADD_of_R500_dim6_pesh_profmoist_Jpce_errorLogN02_truerain53_Jb",str_beta,"Nx_", Nx_boot,".csv", sep=""), header = FALSE)
coeff_multidx_r_TEST_boot <- data.frame(t(ml_coeff_multidx_r_TEST_boot)) 

## Get the degrees of the parameters in the basis
A_multidx <- read.csv(paste('union_reduced_A_basis_Jpce_errorLogN02_truerain53_Jb',str_beta,'_boot.csv',sep=""), header = TRUE)
param_order_A_basis <- paste(data.frame(t(A_multidx)))
param_order_A_basis # to be used for c_alpha column names

## Set column names with basis multi indices
colnames(coeff_multidx_r_TEST_boot) <- c(param_order_A_basis, "boot_idx")
coeff_multidx_r_TEST_boot$boot_idx <- as.factor(coeff_multidx_r_TEST_boot$boot_idx)

## Columns from Figure8 of the paper
cols_fig8 <- c("c(0, 0, 0, 0, 0, 0)", "c(0, 1, 0, 0, 0, 0)", 
               "c(0, 0, 0, 0, 1, 0)", "c(0, 0, 1, 0, 0, 0)", "c(0, 0, 2, 0, 0, 0)")
##    PLOT
var_magnitudes <- sapply(abs(coeff_multidx_r_TEST_boot[1:(length(coeff_multidx_r_TEST_boot) - 1)]), var)
top_columns <- names(sort(var_magnitudes, decreasing = TRUE)[1:9])

## no legend
ggpairs(
  coeff_multidx_r_TEST_boot[, c("boot_idx", cols_fig8)], columns = 2:6,
  aes(color = boot_idx, alpha = 0.5),
  upper = list(continuous = wrap("density", alpha = 0.6), combo = "box"),
  lower = list(continuous = function(data, mapping) {
    ggplot(data, mapping) +
      geom_point(pch = 1, alpha = 0.5, size = 1)}, 
    combo = wrap("dot", alpha = 0.1, size = 0.2))) + theme_bw()

## with legend
#ggpairs(
#  coeff_multidx_r_TEST_boot[, c("boot_idx", cols_fig8)], columns = 2:6,
#  aes(color = boot_idx, alpha = 0.5),
#  upper = list(continuous = wrap("density", alpha = 0.6), combo = "box"),
#  lower = list(continuous = function(data, mapping) {
#    ggplot(data, mapping) +
#      geom_point(pch = 1, alpha = 0.5, size = 1)}, 
#    combo = wrap("dot", alpha = 0.1, size = 0.2)), legend = 1) + theme_bw()+
#  theme(legend.position = 'bottom') +
#  guides(alpha = "none") 
