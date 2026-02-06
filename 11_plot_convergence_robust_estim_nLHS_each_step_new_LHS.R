###################################################################
####   EVALUATE EMULATOR IN NEW POINTS TO GET ROBUST MINIMIZERS
###################################################################
working_directory <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/"
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_functions.R")
setwd(paste(working_directory, "/scripts_PCE/metamodels/", sep=""))
prior_distrib <- read.csv('prior_params_wo_bds145.csv',header = FALSE)
vector_input_indices = c(63,68,72,99,71,65);
df_prior <- data.frame(param = c('th9', 'mn10', 'th10', 'th13', 'thr10', 'hg10'),
                       mean = prior_distrib[vector_input_indices,1],
                       sd =   prior_distrib[vector_input_indices,2])

# Isoprobabilistic transform
IsoProb_z_to_x <- function(z, df_prior){
  z*df_prior$sd[1:5] + df_prior$mean[1:5]
}

# A unique variable to be compatible with optimr package
surrogate_cost_function_analytical_fixedrain_1input <- function(x)
{eval_surrogate_OLS_in_point_x(c(x, -3.50), df_prior, A_multidx, coeff_PCE_fixed_rain)}

# Number of parameters to be optimized
n_params <- 5

## Fix values for excursion sets
excset_val1 <- 0.01; excset_val2 <- 0.02; excset_val3 <- 0.03

## Generate a space-filling design in the parameter space (the MM is evaluated)
N_LHS <- c(0,20,30,50,100,150, 200, 300, 400, 500, 750,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000) # it has to have 0 at the beginning

## load KDE
Nomega <- 500
c_alpha_KDE_new_temp <- c_alpha_KDE_new
c_alpha_KDE_new <- c_alpha_KDE_new_temp[1:Nomega,]
c_alpha_KDE_new_multidx_r <- t(c_alpha_KDE_new)

## Get different Robust estimators with accept/reject algorithm
# prevu :2h00 pour 1000x1000 simulations; commencé à 14h17; fini a 16h17 (2h); 11h55-14h00
# 50..5000 x 500 = 14 650 000 simu = 30h.

for (idx_n_lhs in seq(2,length(N_LHS))){ # starts from 2 to compensate the 0 in the beginning
  ## Initialize sets
  admissible_points_excset1 <- c(); admissible_points_excset2 <- c(); admissible_points_excset3 <- c()
  n_lhs <- N_LHS[idx_n_lhs]
  print(paste('n_lhs = ', n_lhs))
  
  # generate LHS in the parameter space
  z_LHS_quantiles <- randomLHS(n_lhs, n_params)
  z_LHS <- qnorm(z_LHS_quantiles, mean = 0, sd = 1)
  
  # Loop on ALL trajectories
  if(n_lhs %in% c(5500,6000)){
    for (rain_idx in seq(1:length(c_alpha_KDE_new_multidx_r[1,]))){
      print(paste('omega idx = ', rain_idx))
      # Read pce coefficients corresponding to this trajectory
      c_alpha <- c_alpha_KDE_new_multidx_r[,rain_idx]
      # Loop on the space-filling design in the parameter space
      #for (idx_theta in seq(N_LHS[idx_n_lhs-1] + 1, n_lhs)) 
      for (idx_theta in seq(1, n_lhs)) 
        # do not repeat the points in which the MM was already evaluated !
      {
        # Fix one point from the unitary LHS  # Take its value in the parameter space
        zk <- z_LHS[idx_theta,]; theta_k <- c(IsoProb_z_to_x(zk,df_prior = df_prior),-3.50)
        # Evaluate the trajectory in this point
        value_J <- eval_surrogate_OLS_in_point_x(theta_k, df_prior, A_multidx, c_alpha_KDE_new_multidx_r[,rain_idx])
        # Accept the point if it fulfills the criteria: 
        if (value_J < excset_val1){
          admissible_points_excset1 <- rbind(admissible_points_excset1, c(theta_k, rain_idx, idx_theta))}
        if (value_J < excset_val2){
          admissible_points_excset2 <- rbind(admissible_points_excset2, c(theta_k, rain_idx, idx_theta))}
        if (value_J < excset_val3){
          admissible_points_excset3 <- rbind(admissible_points_excset3, c(theta_k, rain_idx, idx_theta))}
      }
    }
    
    # Save results for the number of samples n_lhs.
    df_admissible_points_excset1 <- data.frame(admissible_points_excset1)
    df_admissible_points_excset2 <- data.frame(admissible_points_excset2)
    df_admissible_points_excset3 <- data.frame(admissible_points_excset3)
    colnames(df_admissible_points_excset1) <- c(df_prior$param, 'rain_idx', 'idx_theta')
    colnames(df_admissible_points_excset2) <- c(df_prior$param, 'rain_idx', 'idx_theta')
    colnames(df_admissible_points_excset3) <- c(df_prior$param, 'rain_idx', 'idx_theta')
    
    setwd(excursion_set_results)
    save(df_admissible_points_excset1, file = paste("df_admissible_points_excset",excset_val1,"Nomega",Nomega,"_nLHS_eachLHS",n_lhs,".RData", sep = ""))
    save(df_admissible_points_excset2, file = paste("df_admissible_points_excset",excset_val2,"Nomega",Nomega,"_nLHS_eachLHS",n_lhs,".RData", sep = ""))
    save(df_admissible_points_excset3, file = paste("df_admissible_points_excset",excset_val3,"Nomega",Nomega,"_nLHS_eachLHS",n_lhs,".RData", sep = ""))
  }
}

########################
###   GET ROBUST ESTIMATORS FROM THE ADMISSIBLE SETS
########################
excursion_set_results <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/processed/excursion_sets/"
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")

N_LHS <- c(0,20,30,50,100,150, 200, 300, 400, 500, 750,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500,6000) # it has to have 0 at the beginning
#for (Nomega in c(250,500)){
Nomega <- 500
excset_convergence_nLHS <- c()

for (idx_n_lhs in seq(2,length(N_LHS))){ # starts from 2 to compensate the 0 in the beginning
  n_lhs <- N_LHS[idx_n_lhs]
  print(paste('n_lhs = ', n_lhs))
  ## Load admissible sets for each trajectory realization
  setwd(excursion_set_results)
  load(paste("df_admissible_points_excset",excset_val1,"Nomega",Nomega,"_nLHS_eachLHS",n_lhs,".RData", sep = ""))
  load(paste("df_admissible_points_excset",excset_val2,"Nomega",Nomega,"_nLHS_eachLHS",n_lhs,".RData", sep = ""))
  load(paste("df_admissible_points_excset",excset_val3,"Nomega",Nomega,"_nLHS_eachLHS",n_lhs,".RData", sep = ""))
  
  ## Count how many times each theta_idx appears in the excursion sets
  value_counts1 <- table(df_admissible_points_excset1$idx_theta)
  value_counts2 <- table(df_admissible_points_excset2$idx_theta)
  value_counts3 <- table(df_admissible_points_excset2$idx_theta)
  
  ## Which theta value is the most frequent
  most_frequent_value1 <- names(value_counts1)[which.max(value_counts1)]
  most_frequent_value2 <- names(value_counts2)[which.max(value_counts2)]
  most_frequent_value3 <- names(value_counts2)[which.max(value_counts2)]
  
  theta_k1 <- df_admissible_points_excset1[which(df_admissible_points_excset1$idx_theta == most_frequent_value1)[1],]
  theta_k2 <- df_admissible_points_excset2[which(df_admissible_points_excset2$idx_theta == most_frequent_value2)[1],]
  theta_k3 <- df_admissible_points_excset3[which(df_admissible_points_excset3$idx_theta == most_frequent_value3)[1],]
  
  # Print the most frequent value
  print(paste(theta_k1[1:5]))
  print(paste(theta_k2[1:5]))
  print(paste(theta_k3[1:5]))
  
  # concatenate results
  excset_convergence_nLHS <- rbind(excset_convergence_nLHS, c(n_lhs,excset_val1,  theta_k1[1:5]))
  excset_convergence_nLHS <- rbind(excset_convergence_nLHS, c(n_lhs,excset_val2,  theta_k2[1:5]))
  excset_convergence_nLHS <- rbind(excset_convergence_nLHS, c(n_lhs,excset_val3,  theta_k3[1:5]))
}

df_excset_convergence_nLHS <- data.frame(excset_convergence_nLHS)
str(df_excset_convergence_nLHS)
colnames(df_excset_convergence_nLHS) <- c('nLHS', 'Excursion.set', 'th9','mn10','th10','th13','thr10' )
df_excset_convergence_nLHS$nLHS <- as.character(df_excset_convergence_nLHS$nLHS)
df_excset_convergence_nLHS$'Excursion.set' <- as.character(df_excset_convergence_nLHS$'Excursion.set')
df_excset_convergence_nLHS$th9 <- as.character(df_excset_convergence_nLHS$th9)
df_excset_convergence_nLHS$mn10 <- as.character(df_excset_convergence_nLHS$mn10)
df_excset_convergence_nLHS$th10 <- as.character(df_excset_convergence_nLHS$th10)
df_excset_convergence_nLHS$th13 <- as.character(df_excset_convergence_nLHS$th13)
df_excset_convergence_nLHS$thr10 <- as.character(df_excset_convergence_nLHS$thr10)

colnames(df_excset_convergence_nLHS) <- factor(c('nLHS', 'Excursion.set',"$\\theta_{s.inter.}$",
                                                 "$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep}$"),
                                               levels = c('nLHS', 'Excursion.set',"$\\theta_{s.inter.}$","$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep}$"))
df_prior$param <-factor(c("$\\theta_{s.inter.}$","$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep}$", "hg"),
                        levels = c("$\\theta_{s.inter.}$","$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep}$"))

head(df_excset_convergence_nLHS)

dfdf_excset_convergence_nLHS <- melt(df_excset_convergence_nLHS, id.vars = c('nLHS', 'Excursion.set'))
head(dfdf_excset_convergence_nLHS)
dfdf_excset_convergence_nLHS$value <- as.double(dfdf_excset_convergence_nLHS$value)
dfdf_excset_convergence_nLHS$nLHS <- as.double(dfdf_excset_convergence_nLHS$nLHS)

## Add limits to the plot # Merge the dataframes on the specified columns
dfdf_excset_convergence_nLHS <- dfdf_excset_convergence_nLHS %>%
  left_join(df_prior, by = c("variable" = "param"))

## Plot the convergence
gg_nomega <- ggplot(dfdf_excset_convergence_nLHS[dfdf_excset_convergence_nLHS$Excursion.set %in% c(0.01,0.02),],
                    aes(x=nLHS, color = Excursion.set, y = value, group = Excursion.set))+
  facet_wrap(vars(variable), scales = "free_y") + geom_point(aes(color = Excursion.set)) + 
  geom_line(aes(color = Excursion.set)) + theme_bw() + ylab('') + xlab('card($\\mathcal{X}_{candidates}$)') +
  geom_hline(aes(yintercept=mean),linetype="solid", colour="#beaed4") +
  geom_hline(aes(yintercept=mean + 3*sd),linetype="dashed", colour="gray") +
  geom_hline(aes(yintercept=mean - 3*sd),linetype="dashed", colour="gray") +
  scale_color_manual('Excursion set threshold',values=c('0.01'='#386cb0', '0.02'='#f0027f'))+
  #ggtitle(paste('Estimated parameter from MM when nKDE = ',Nomega,' and varying Nx', sep = '')) +
  theme(legend.position = "bottom")   
gg_nomega

setwd("~/Documents/REDACTION/manuscrit_radisic/5_CALIB_MOIST_PROF/figures/")
tikz(paste('gg_convergence_Nomega',Nomega,'.tex', sep = ''), standAlone = TRUE, width=6.0, height=4.0)
gg_nomega
dev.off()
tools::texi2dvi(paste('gg_convergence_Nomega',Nomega,'.tex', sep = ''),pdf=T)

#pdf(file = paste('gg_convergence_Nomega',Nomega,'.pdf', sep = ''))
#plot(gg_nomega)
#dev.off()
#}

