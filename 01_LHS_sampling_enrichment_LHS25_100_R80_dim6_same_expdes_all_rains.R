########################
###   LHS sampling   ###
## add new 100 lhs points to the test set.
########################
rm(list=ls()) 
library('lhs')
library("RcppCNPy") # needed to read .npy
set.seed(111222)

n.dim <- 6
test_set_r80 <- c()

## Read old LHS to enrich
setwd("/Users/icj_ecl/Documents/RESS_review_depo/")
old_LHS_test <- as.matrix(npyLoad("LHS25_100_R500_dim6_sample_final.npy"))

## Prepare columns whose values should not be changed
factors  <- c(63,65,68,71,72,99)
sorted_imp_factors <- sort(factors, decreasing = FALSE)
all_factors <- c(1:145)
unimportant_factors <- setdiff(all_factors,factors)

# open priors from Emilie marginals
prior_distrib <- read.csv('prior_params_wo_bds145.csv',header = FALSE)

## Transform old LHS to its quantiles
old_LHS_test_cdf <- c()
for (idx_param in sorted_imp_factors){
  old_LHS_test_cdf <- cbind(old_LHS_test_cdf, mapply(pnorm, old_LHS_test[1:100,idx_param], 
                                   prior_distrib[idx_param,1], prior_distrib[idx_param,2]))
}

## Augment LHS cdf
new_LHS_test_cdf <- augmentLHS(old_LHS_test_cdf, 200-100) # test set

## Transform to parameter values
new_LHS_test <- c()
for (idx in seq(1,6)){
  idx_param <- sorted_imp_factors[idx]
  new_LHS_test <- cbind(new_LHS_test, mapply(qnorm, new_LHS_test_cdf[,idx], 
                                                     prior_distrib[idx_param,1], prior_distrib[idx_param,2]))
}

## Visualize distributions
hist(new_LHS_test[,1])
hist(new_LHS_test[,2])
hist(new_LHS_test[,3])
hist(new_LHS_test[,4])
hist(new_LHS_test[,5])
hist(new_LHS_test[,6])

## SAME LHS for each rain

## Take the new Nx points, with only 10 old Nx points (to ensure that these 10 are the same as old simulations) 
Nx_new_simus <- 110
sub_new_LHS_test <- new_LHS_test[(200-Nx_new_simus+1):200,]

#########
## only with 300 rains
r.test <- 300

for (r_idx in seq(1,r.test)){
  # set of all parameters, those that are not influential are set to the values from old.
  test_set <- old_LHS_test[1:Nx_new_simus,]
  i = 1
  for (idx in sorted_imp_factors){
    test_set[,idx] <- sub_new_LHS_test[,i]
    i = i + 1
  }
  test_set_r80 <- rbind(test_set_r80,test_set)
}

write.csv(test_set_r80,"LHS110_200_R300_dim6_sample_final.csv", row.names = FALSE)
npySave("LHS110_200_R300_dim6_sample_final.npy", test_set_r80)
