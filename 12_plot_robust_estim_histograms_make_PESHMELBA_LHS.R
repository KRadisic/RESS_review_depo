##################################################################
###   OPEN ROBUSTNESS RESULTS (after convergence scripts)
###   - plot the values on the parameter space and compare with classical calibration
###   - save parameter values in format for new PESHMELBA simulations
##################################################################
rm(list=ls())
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_functions.R")

excset_val1 <- 0.01
excset_val2 <- 0.02

## Number of simulations of the metamodel for choosing robust calibration
Nomega <- 500
n_lhs <- 4500
str_beta <- '005'
beta <- 0.005

## Load admissible sets for each trajectory realization
excursion_set_results <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/processed/excursion_sets/"
setwd(excursion_set_results)
load(paste("df_admissible_points_excset",excset_val1,"Nomega",Nomega,"_nLHS_eachLHS",n_lhs,".RData", sep = ""))
load(paste("df_admissible_points_excset",excset_val2,"Nomega",Nomega,"_nLHS_eachLHS",n_lhs,".RData", sep = ""))

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

##########################################################
#####    PLOT CONDITIONAL MINIMIZERS 
#########################################################
setwd("~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/scripts_PCE/metamodels/")
df_min_LARS <- read.csv('conditional_minimizers_LARS_errorLogN02_Jb_beta005LHS_rep5_50_R500_Ntrain_hourly_BDDquantiles50.csv', header = TRUE)[1:200,]
#df_min_LARS <- read.csv('conditional_minimizers_OLS_errorLogN02_Jb_beta0.005_LHS_rep5_50_R500_Ntrain_hourly_BDDquantiles50.csv', header = TRUE) 

# results of the conditional one shot minimization
head(df_min_LARS)
colnames(df_min_LARS) <- c("$\\theta_{s.inter.}$", "$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep}$", " $hg_{.deep}$", "f")

R <- 200; n_params <- 6; true_test_idx <- 1;
df_cond_minimizers_LARS <- data.frame(cond.min = as.vector(as.matrix(df_min_LARS)[seq(n_params,R),]),
                                      parameter = rep(colnames(df_min_LARS), each = length(seq(n_params,R))))

df_cond_minimizers_LARS$type <- '$\\hat{x}_{s\\_PCE}^{*(r)}$'
df_cond_min <- df_cond_minimizers_LARS

######  ROBUST AND TRUE RESULTS AND FORMATTING FOR PLOT : VERTICAL LINES
## Load prior add true in dataframe for plotting
working_directory <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/"
setwd(paste(working_directory, "/scripts_PCE/metamodels/", sep=""))
prior_distrib <- read.csv('prior_params_wo_bds145.csv',header = FALSE)
vector_input_indices = c(63,68,72,99,71,65);
df_prior <- data.frame(param = c('th9', 'mn10', 'th10', 'th13', 'thr10', 'hg10'),
                       mean = prior_distrib[vector_input_indices,1],
                       sd =   prior_distrib[vector_input_indices,2])

## Load previous results from other robust calibrations ... Load true
setwd('~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/raw/')
x_test <- npyLoad("LHS25_100_R198_dim6_sample_final.npy")
x_true <- x_test[1,c(63,68,72,99,71,65)] # matlab order

n_params_plot <- 5
df_prior_true_mean <- df_prior[1:n_params_plot,]
df_prior_true_mean$param <- c("$\\theta_{s.inter.}$", "$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep}$")
df_prior_true_mean$true_parameter <- x_true[1:n_params_plot]
#df_prior_true_mean$true_rain <- c(as.matrix(df_min_LARS)[1,1:n_params_plot])
df_prior_true_mean$exc_set1 <- as.numeric(theta_k1[1:n_params_plot])
df_prior_true_mean$exc_set2 <- as.numeric(theta_k2[1:n_params_plot])
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
    scale_fill_manual(name = "Conditional minimizers",values=c("#7fc97f")) + xlab('') + ylab('') +
    facet_grid(. ~ parameter, scales = "free_x" ) + theme_bw() + theme(legend.position="none") +
    stat_function(fun = dnorm, args = list(mean = prior_distrib_param$mean, sd = prior_distrib_param$sd), color ="blue") +
    #geom_vline(data = prior_distrib_param,aes(xintercept = true_parameter, color="$x_{true}$"),
    #           linetype="dashed",linewidth = 1.0)+
    geom_vline(data = prior_distrib_param,aes(xintercept = mean, color="$\\mu_{prior}$"),
               linetype="dashed",linewidth = 1.0)+
    geom_vline(data = prior_distrib_param,aes(xintercept = exc_set1, color=paste("${x}^*_{\\hat{f}_s<", excset_val1,"}$", sep = "")),
               linetype="solid",linewidth = 1.0)+
    geom_vline(data = prior_distrib_param,aes(xintercept = exc_set2, color=paste("${x}^*_{\\hat{f}_s<", excset_val2,"}$", sep = "")),
               linetype="solid",linewidth = 1.0)+
  scale_color_manual(name = "", values = c("$\\mu_{prior}$" = "#beaed4",#"$x_{true}$" = "red",
                                           "${x}^*_{\\hat{f}_s<0.01}$" = "#386cb0",
                                           "${x}^*_{\\hat{f}_s<0.02}$" = "#f0027f"))
  list_of_plots_excset[[z]] <- conditional_minimizers
  z <- z+1
}
## add legend as sixth plot # Using the cowplot package
conditional_minimizers <- conditional_minimizers + theme(legend.position="left")
legend <- get_legend(conditional_minimizers)
list_of_plots_excset[[z]] <- legend
conditional_minimizers_over_prior_true <- do.call(grid.arrange, c(list_of_plots_excset, ncol = 3))
conditional_minimizers_over_prior_true

###################################################################
###   Figure 5.13
###################################################################
## to tikZ
#setwd("~/Documents/REDACTION/manuscrit_radisic/5_CALIB_MOIST_PROF/figures/")
#tikz(paste('gg_arranged_cond_min_robust_estim_excset_PCE_LogN02_Jb',str_beta,'_YzeronS06.tex', sep=''), standAlone = TRUE, width=5.5, height=4.0)
#print(grid::grid.draw(conditional_minimizers_over_prior_true)) # does not add a page alone
#dev.off()
#tools::texi2dvi(paste('gg_arranged_cond_min_robust_estim_excset_PCE_LogN02_Jb',str_beta,'_YzeronS06.tex', sep=''),pdf=T)

###################################################################
###    SAVE MINIMIZATION RESULTS FOR NEW PESHMELBA SIMULATIONS
###################################################################
## Arrange minimizers for new PESHMELBA simulations
r.new.rains <- 500 # number of rains for same parameter values
n.dim <- 5 # number of parameter values per rain

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
single_LHS[3, vector_input_indices] <- as.numeric(c(df_prior_true_mean$exc_set1, -3.50))
single_LHS[4, vector_input_indices] <- as.numeric(c(df_prior_true_mean$exc_set2, -3.50))
single_LHS[5, vector_input_indices] <- as.numeric(c(df_prior_true_mean$mean, -3.50))

#[1,] 0.3251800 0.1655000 0.3391800 0.2332000 0.07157900 -2.00
#[2,] 0.3118763 0.1741481 0.3419935 0.2739860 0.06492179 -3.50
#[3,] 0.3026490 0.1788592 0.3446818 0.2889340 0.06320647 -3.50 # <0.01
#[4,] 0.3098395 0.1771087 0.3415738 0.2991149 0.07514100 -3.50 # <0.02
#[5,] 0.3322500 0.1791000 0.3160100 0.3375000 0.06117900 -3.56

head(single_LHS[,vector_input_indices])

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
npySave(paste("LHS_robust_minimizers_LogN02_beta005_condmin",idx_of_conditional_minimizer,
              "_Nomega",Nomega,"_nlhs",n_lhs,
              "_dim6_sample_final.npy", sep=""), training_set_r)


