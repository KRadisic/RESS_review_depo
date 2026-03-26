##########################################################
#####   PLOT NEW SIMULATIONS OF THE COST FUNCTION
#####   AFTER THE ROBUST MINIMIZATION
rm(list=ls())

library('ggplot2')
library('GGally')
library('RcppCNPy') # open .npy data
library('tikzDevice')
library("reshape2") # melt
library('ggpubr') # ggarrange

working_directory <- '~/Documents/RESS_review_depo/'

source("~/Documents/RESS_review_depo/0_functions.R")

## Number of simulations of the metamodel for choosing robust calibration
Nomega <- 500
n_lhs <- 4500
beta <- 0.005
idx_of_conditional_minimizer <- 174 

setwd(working_directory)
df_minimizers <- npyLoad(paste("LHS_robust_minimizers_LogN02_beta005_condmin",idx_of_conditional_minimizer,
                               "_Nomega",Nomega,"_nlhs",n_lhs,
                               "_dim6_sample_final.npy", sep=""))
y_new_simulations <- npyLoad(paste("Ymoisture_profile_calib_relreg1_meanprior_dim6_YzeronS06_LogN02.npy", sep=""))

vector_input_indices <- c(63, 65, 68, 71, 72, 99)
df_minimizers[c(1,2,3,4,5,6),vector_input_indices] #df_min_LARS[idx_of_conditional,x_bar_PCE,relreg1,relreg3, exc_set1,exc_set2

# Ordered names of minimizers (in the way they are ordered in 7_get_PCA_from_train_set_sigma02_compare_histograms.R)
minimizer_names_ordered <- c('$\\hat{x}^{*(r_i)}_{s\\_PCE}$',
                             '$\\bar{x}^{*}$',
                             '${x}^*_{\\hat{f}_s<0.01}$',
                             '${x}^*_{\\hat{f}_s<0.02}$',                             
                             '$\\mu_{prior}^{*}$')

colors_minimizers <- c('$\\hat{x}^{*(r_i)}_{s\\_PCE}$' = "#7fc97f",
                       '$\\bar{x}^{*}$' = "#fdc086",
                       '${x}^*_{\\hat{f}_s<0.01}$' = "#386cb0",
                       '${x}^*_{\\hat{f}_s<0.02}$' = "#f0027f",
                       '$\\mu_{prior}^{*}$' = "#beaed4")
## Load true simulation
true_test_idx <- 1
true_rain_idx <- 53
y_true <- npyLoad("Ymoisture_profile_LHS25_100_R198_dim6_YzeronS06_hourly.npy")[true_rain_idx*100 + true_test_idx,]


# open priors from Emilie marginals
prior_distrib <- read.csv('prior_params_wo_bds145.csv',header = FALSE)
vector_input_indices = c(63,68,72,99,71,65);
df_prior <- data.frame(param = c('th9', 'mn10', 'th10', 'th13', 'thr10', 'hg10'),
                       mean = prior_distrib[vector_input_indices,1],
                       sd =   prior_distrib[vector_input_indices,2])

## Calculate cost function (squared diff)
J0 <- c() 
for (idx_row in 1:2500){
  J0 <- c(J0, sum((y_new_simulations[idx_row,] - y_true)^2))
}

## Calculate regularization term
Jb <- beta*(((df_minimizers[,68] - df_prior$mean[2])/df_prior$sd[2])^2 +
              ((df_minimizers[,71] - df_prior$mean[5])/df_prior$sd[5])^2)
## Transform to data frame for plotting
J <- J0+ Jb
n.rains.test <- 500
rep.min.types <- rep(minimizer_names_ordered, n.rains.test)
df_cost_functions <- data.frame(Cost.function = J, minimizer = rep.min.types)

#############################################
#######     PLOTS OF HISTOGRAMS OF J0
#############################################
new_cost_function <- ggplot(df_cost_functions[df_cost_functions$minimizer %in% minimizer_names_ordered[c(1,3,4,5)],], 
                            aes(x = Cost.function, color = minimizer, pattern=minimizer))+ xlab('Cost function $f_s(x^{*}, \\cdot)$') + ylab('') +
  geom_histogram(position = "identity", aes(y = ..density.., fill = minimizer),alpha = 0.50, bins = 50, color="black")+
  scale_fill_manual('Calibrated values',values=colors_minimizers) + theme_bw() 
new_cost_function

####################################################
######     MAKE TABLE OF STATS OF NEW COST FUNCTION
####################################################
stats_of_new_costfct <- c()
for (idx_of_robust_min in seq(1, length(minimizer_names_ordered)))
{ stats_of_new_costfct <- 
  rbind(stats_of_new_costfct,
        c(
          round(mean(df_cost_functions[df_cost_functions$minimizer == minimizer_names_ordered[idx_of_robust_min],]$Cost.function, na.rm = TRUE),digits=4),
          round(sum(df_cost_functions[df_cost_functions$minimizer == minimizer_names_ordered[idx_of_robust_min],]$Cost.function>0.01, na.rm = TRUE)/n.rains.test,digits=2),
          round(sum(df_cost_functions[df_cost_functions$minimizer == minimizer_names_ordered[idx_of_robust_min],]$Cost.function>0.02, na.rm = TRUE)/n.rains.test,digits=2),
          round(max(df_cost_functions[df_cost_functions$minimizer == minimizer_names_ordered[idx_of_robust_min],]$Cost.function, na.rm = TRUE),digits=3),
          round(var(df_cost_functions[df_cost_functions$minimizer == minimizer_names_ordered[idx_of_robust_min],]$Cost.function, na.rm = TRUE),digits=7)))
}
df_stats_of_new_costfct <- data.frame(stats_of_new_costfct)
rownames(df_stats_of_new_costfct) <- minimizer_names_ordered
colnames(df_stats_of_new_costfct) <- c('$\\mu(f(\\cdot,\\cdot))$', 'n>0.01', 'n>0.02', 'max', 'var')
head(df_stats_of_new_costfct)

#setwd(processed_data_directory)
#write.csv(df_stats_of_new_costfct, file = 'table_results_robust_simulations_J0.csv')

#####################################################
#######     PREPARE DATAFRAME FOR PLOTTING
#####################################################
df_y_new_simulations <- data.frame(y_new_simulations)
df_y_new_simulations$minimizer <- rep.min.types

colnames(df_y_new_simulations)<-c(seq(1,25),'minimizer')
df_y_new_simulations$rain.idx <- factor(rep(seq(1,n.rains.test), each = length(minimizer_names_ordered)))

dfdf_y_new_simulations <- melt(df_y_new_simulations, id.vars = c("minimizer", "rain.idx"))
dfdf_y_new_simulations$variable <- as.numeric(dfdf_y_new_simulations$variable)

#############################################
#######     PLOTS OF MOISTURE PROFILES
#############################################
new_profiles <- ggplot(data = dfdf_y_new_simulations[dfdf_y_new_simulations$minimizer %in% minimizer_names_ordered[c(1,3,4,5)], ], 
                       aes(x = variable, y = value, group = interaction(minimizer,rain.idx))) + 
  scale_color_manual('Calib. value',values=colors_minimizers) + coord_flip() + scale_x_reverse() +
  geom_line(aes(color = minimizer), alpha = 0.50, linetype='solid') + 
  geom_point(aes(x=1,y=y_true[1])) +  geom_point(aes(x=2,y=y_true[2])) +
  geom_point(aes(x=3,y=y_true[3])) +  geom_point(aes(x=4,y=y_true[4])) +
  geom_point(aes(x=5,y=y_true[5])) +  geom_point(aes(x=6,y=y_true[6])) +
  geom_point(aes(x=7,y=y_true[7])) +  geom_point(aes(x=8,y=y_true[8])) +
  geom_point(aes(x=9,y=y_true[9])) +  geom_point(aes(x=10,y=y_true[10])) +
  geom_point(aes(x=11,y=y_true[11])) +  geom_point(aes(x=12,y=y_true[12])) +
  geom_point(aes(x=13,y=y_true[13])) +  geom_point(aes(x=14,y=y_true[14])) +
  geom_point(aes(x=15,y=y_true[15])) +  geom_point(aes(x=16,y=y_true[16])) +
  geom_point(aes(x=17,y=y_true[17])) +  geom_point(aes(x=18,y=y_true[18])) +
  geom_point(aes(x=19,y=y_true[19])) +  geom_point(aes(x=20,y=y_true[20])) +
  geom_point(aes(x=21,y=y_true[21])) +  geom_point(aes(x=22,y=y_true[22])) +
  geom_point(aes(x=23,y=y_true[23])) +  geom_point(aes(x=24,y=y_true[24])) +
  geom_point(aes(x=25,y=y_true[25])) +
  geom_vline(xintercept = 0.5, color="#FE9929",linetype="dashed",linewidth = 1.0)+
  geom_vline(xintercept = 2.5, color="#FE9929",linetype="dashed",linewidth = 1.0)+
  geom_vline(xintercept = 3, color="#D95F0E",linetype="dashed",linewidth = 1.0)+
  geom_vline(xintercept = 7.5, color="#D95F0E",linetype="dashed",linewidth = 1.0)+
  geom_vline(xintercept = 8, color="#993404",linetype="dashed",linewidth = 1.0)+
  geom_vline(xintercept = 25, color="#993404",linetype="dashed",linewidth = 1.0)+
  #ggtitle('PESHMELBA after calibration') +
  theme_bw() + xlab('Cell') + ylab('Soil moisture [cm3/cm3]')
new_profiles

ggarrange(new_cost_function, new_profiles, ncol=2, nrow=1, common.legend = TRUE, legend="right")

##########################
###    SAVE PLOTS
#Sys.setenv(PATH = paste("/usr/bin", Sys.getenv("PATH"), sep=":"))
## to tikZ
#setwd(working_directory)
#tikz('new_profiles_and_costfct.tex', standAlone = TRUE, width=7.5, height=3.0)
#ggarrange(new_cost_function, new_profiles, ncol=2, nrow=1, common.legend = TRUE, legend="right")
#dev.off()
#tools::texi2dvi('new_profiles_and_costfct.tex',pdf=T) ## if not enough memory, run with lualatex

#tikz('new_cost_function.tex', standAlone = TRUE, width=4.0, height=2.5)
#new_cost_function
#dev.off()
#tools::texi2dvi('new_cost_function.tex',pdf=T)




