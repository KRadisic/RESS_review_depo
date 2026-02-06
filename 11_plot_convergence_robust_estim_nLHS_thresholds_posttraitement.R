###################################################################
####   EVALUATE EMULATOR IN NEW POINTS TO GET ROBUST MINIMIZERS
###################################################################
rm(list=ls())

working_directory <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/"
excursion_set_results <- "~/Documents/Cours_presentations_formations/4_Presentation_poster_article/2_Conferences/JMSC_strasbourg_2024/data/processed/excursion_sets/"

source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_libraries.R")
source("~/code/1_code_article_PhysicaD_lambdaOK_YzeronS06/SCRIPTS/0_functions.R")
setwd(paste(working_directory, "/scripts_PCE/metamodels/", sep=""))
prior_distrib <- read.csv('prior_params_wo_bds145.csv',header = FALSE)
vector_input_indices = c(63,68,72,99,71,65);
df_prior <- data.frame(param = c('th9', 'mn10', 'th10', 'th13', 'thr10', 'hg10'),
                       mean = prior_distrib[vector_input_indices,1],
                       sd =   prior_distrib[vector_input_indices,2])

Nomega <- 500
## LOAD results (local, already concatenated)
setwd(excursion_set_results)
#n_lhs <- 6000
#load(paste("df_MM_realiz_Nomega",Nomega,"_nLHS",n_lhs,".RData", sep = ""))

## LOAD results (cluster, concatenate)
dfdf_dataframe <- c()
for (n_lhs in c(20, 30, 50, 100, 150, 200, 300, 750, 1000, 1500, 2000, 2500,
                3000, 3500, 4000, 4500, 5000, 5500, 6000, 7500))
{
  load(paste("df_interm_MM_realiz_Nomega",Nomega,"_nLHS",n_lhs,".RData", sep = ""))
  # Concatenate with previous cardinalities of the parameter space
  dfdf_dataframe <- rbind(dfdf_dataframe, df_dataframe)
}

########################
###   GET ROBUST ESTIMATORS FROM THE ADMISSIBLE SETS
########################
head(dfdf_dataframe)

# Set the threshold value
dfdf_excset <- c()
threshold_values <- seq(0.005,0.030,0.0025)
for (c_threshold in threshold_values){
  ## Get for each X_card the best parameter.
  result <- dfdf_dataframe %>%
    filter(value_J < c_threshold) %>% # keep only those under the threshold
    group_by(n_candidates, idx_theta) %>% # cound how many times each theta appears (for a fixed n_candidated size)
    summarise(count = n()) %>% 
    group_by(n_candidates) %>% 
    slice_max(count) %>%
    ungroup()
  
  # Print the results
  print(result)
  df_excset <- data.frame(n_candidates = result$n_candidates, idx_theta = result$idx_theta, 
                          c_threshold = c_threshold, prob.acceptance = result$count)
  dfdf_excset <- rbind(dfdf_excset, df_excset) 
}
head(dfdf_excset)

dfdf_excset$th9 <- NaN; 
dfdf_excset$mn10 <- NaN; 
dfdf_excset$th10 <- NaN; 
dfdf_excset$th13 <- NaN; 
dfdf_excset$thr10 <- NaN; 

## Put the corresponding values of the parameters.
for (idx_rows in seq(1, length(dfdf_excset$n_candidates))){
  dfdf_excset$th9[idx_rows] <- unique(dfdf_dataframe[which(dfdf_dataframe$idx_theta == dfdf_excset$idx_theta[idx_rows] &
                                                             dfdf_dataframe$n_candidates == dfdf_excset$n_candidates[idx_rows]), 'th9']); 
  dfdf_excset$mn10[idx_rows] <- unique(dfdf_dataframe[which(dfdf_dataframe$idx_theta == dfdf_excset$idx_theta[idx_rows] &
                                                              dfdf_dataframe$n_candidates == dfdf_excset$n_candidates[idx_rows]), 'mn10']); 
  dfdf_excset$th10[idx_rows] <- unique(dfdf_dataframe[which(dfdf_dataframe$idx_theta == dfdf_excset$idx_theta[idx_rows] &
                                                              dfdf_dataframe$n_candidates == dfdf_excset$n_candidates[idx_rows]), 'th10']);  
  dfdf_excset$th13[idx_rows] <- unique(dfdf_dataframe[which(dfdf_dataframe$idx_theta == dfdf_excset$idx_theta[idx_rows] &
                                                              dfdf_dataframe$n_candidates == dfdf_excset$n_candidates[idx_rows]), 'th13']); 
  dfdf_excset$thr10[idx_rows] <- unique(dfdf_dataframe[which(dfdf_dataframe$idx_theta == dfdf_excset$idx_theta[idx_rows] &
                                                               dfdf_dataframe$n_candidates == dfdf_excset$n_candidates[idx_rows]), 'thr10']); 
}

head(dfdf_excset)

## Transform for plotting
df_excset_convergence_nLHS <- dfdf_excset[,c('n_candidates', 'c_threshold', 'prob.acceptance', 'th9','mn10','th10','th13','thr10')]
str(df_excset_convergence_nLHS)

colnames(df_excset_convergence_nLHS) <- c('nLHS', 'Excursion.set', 'prob.acceptance', 'th9','mn10','th10','th13','thr10' )
df_excset_convergence_nLHS$nLHS <- as.character(df_excset_convergence_nLHS$nLHS)
df_excset_convergence_nLHS$'Excursion.set' <- as.character(df_excset_convergence_nLHS$'Excursion.set')
df_excset_convergence_nLHS$'prob.acceptance' <- as.character(df_excset_convergence_nLHS$'prob.acceptance')
df_excset_convergence_nLHS$th9 <- as.character(df_excset_convergence_nLHS$th9)
df_excset_convergence_nLHS$mn10 <- as.character(df_excset_convergence_nLHS$mn10)
df_excset_convergence_nLHS$th10 <- as.character(df_excset_convergence_nLHS$th10)
df_excset_convergence_nLHS$th13 <- as.character(df_excset_convergence_nLHS$th13)
df_excset_convergence_nLHS$thr10 <- as.character(df_excset_convergence_nLHS$thr10)

colnames(df_excset_convergence_nLHS) <- factor(c('nLHS', 'Excursion.set','prob.acceptance', "$\\theta_{s.inter.}$",
                                                 "$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep}$"),
                                               levels = c('nLHS', 'Excursion.set','prob.acceptance',"$\\theta_{s.inter.}$","$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep}$"))
df_prior$param <-factor(c("$\\theta_{s.inter.}$","$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep}$", "hg"),
                        levels = c("$\\theta_{s.inter.}$","$mn_{.deep}$", "$\\theta_{s.deep}$", "$\\theta_{s.surf.}$", "$\\theta_{r.deep}$"))

head(df_excset_convergence_nLHS)

dfdf_excset_convergence_nLHS <- melt(df_excset_convergence_nLHS, id.vars = c('nLHS', 'Excursion.set','prob.acceptance'))
head(dfdf_excset_convergence_nLHS)
dfdf_excset_convergence_nLHS$value <- as.double(dfdf_excset_convergence_nLHS$value)
dfdf_excset_convergence_nLHS$nLHS <- as.double(dfdf_excset_convergence_nLHS$nLHS)

## Add limits to the plot # Merge the dataframes on the specified columns
dfdf_excset_convergence_nLHS <- dfdf_excset_convergence_nLHS %>%
  left_join(df_prior, by = c("variable" = "param"))

#######################
###    PLOT ROBUST ESTIMATE WRT cardinality of exp design
#######################
## Plot the convergence wrt N_card, for different threshold lines
plotting_thresholds <- threshold_values #c(0.01,0.02)
gg_nomega <- ggplot(dfdf_excset_convergence_nLHS[dfdf_excset_convergence_nLHS$Excursion.set %in% c(0.01,0.02),],
                    aes(x=nLHS, color = Excursion.set, y = value, group = Excursion.set))+
  facet_wrap(vars(variable), scales = "free_y") + geom_point(aes(color = Excursion.set)) + 
  geom_line(aes(color = Excursion.set)) + theme_bw() + ylab('') + xlab('card($\\mathcal{X}_{candidates}$)') +
  geom_hline(aes(yintercept=mean),linetype="solid", colour="#beaed4") +
  geom_hline(aes(yintercept=mean + 3*sd),linetype="dashed", colour="gray") +
  geom_hline(aes(yintercept=mean - 3*sd),linetype="dashed", colour="gray") +
  scale_color_manual('Threshold $c$',values=c('0.01'='#386cb0', '0.02'='#f0027f'))+
  #ggtitle(paste('Estimated parameter from MM when nKDE = ',Nomega,' and varying Nx', sep = '')) +
  theme(legend.position = "bottom")   
gg_nomega

#######################
###    SAVE PLOT
#######################
#setwd("~/Documents/REDACTION/manuscrit_radisic/5_CALIB_MOIST_PROF/figures/")
#tikz(paste('gg_convergence_Nomega',Nomega,'.tex', sep = ''), standAlone = TRUE, width=6.0, height=4.0)
#gg_nomega
#dev.off()
#tools::texi2dvi(paste('gg_convergence_Nomega',Nomega,'.tex', sep = ''),pdf=T)

#######################
###    PLOT ROBUST ESTIMATE WRT THRESHOLD
#######################
modif_dfdf_excset_convergence_nLHS <- dfdf_excset_convergence_nLHS
modif_dfdf_excset_convergence_nLHS$nLHS <- factor(dfdf_excset_convergence_nLHS$nLHS)
modif_dfdf_excset_convergence_nLHS$Excursion.set <- as.double(dfdf_excset_convergence_nLHS$Excursion.set)
modif_dfdf_excset_convergence_nLHS$prob.acceptance <- as.double(dfdf_excset_convergence_nLHS$prob.acceptance)

## Plot robust choice wrt threshold, for different N_card lines
gg_nomega <- ggplot(modif_dfdf_excset_convergence_nLHS[modif_dfdf_excset_convergence_nLHS$Excursion.set %in% plotting_thresholds &
                                                         modif_dfdf_excset_convergence_nLHS$nLHS %in% c(4500,5000,5500,6000),],
                    aes(x=Excursion.set, color = nLHS, y = value, group = nLHS))+
  facet_wrap(vars(variable), scales = "free_y") + geom_point(aes(color = nLHS)) + 
  geom_line(aes(color = nLHS)) + theme_bw() + xlab('Threshold $c$') +
  geom_hline(aes(yintercept=mean),linetype="solid", colour="#beaed4") +
  geom_hline(aes(yintercept=mean + 3*sd),linetype="dashed", colour="gray") +
  geom_hline(aes(yintercept=mean - 3*sd),linetype="dashed", colour="gray") +
  scale_color_manual('card($\\mathcal{X}_{candidates}$)',
                     values=c('4500'='#a6cee3', '5000'='#1f78b4', '5500'='#b2df8a', '6000'='#33a02c')#, '7500'='#016c59'
  )+
  theme(legend.position = "bottom")   
gg_nomega

## Plot robust choice wrt threshold, for different N_card lines
gg_nomega_proba <- ggplot(modif_dfdf_excset_convergence_nLHS[modif_dfdf_excset_convergence_nLHS$Excursion.set %in% plotting_thresholds &
                                                         modif_dfdf_excset_convergence_nLHS$nLHS %in% c(4500,5000,5500,6000),],
                    aes(x=Excursion.set, color = nLHS, y = prob.acceptance/Nomega, group = nLHS))+
  #facet_wrap(vars(variable), scales = "free_y") + 
  geom_point(aes(color = nLHS)) + 
  geom_line(aes(color = nLHS)) + theme_bw() + xlab('Threshold $c$') + ylab('$k_{i^*}/R\_{new}$') +
  scale_color_manual('card($\\mathcal{X}_{candidates}$)',
                     #values=c('4500'='#bdc9e1', '5000'='#67a9cf', '5500'='#1c9099', '6000'='#016c59', '7500'='#016f63')
                     values=c('4500'='#a6cee3', '5000'='#1f78b4', '5500'='#b2df8a', '6000'='#33a02c') # , '7500'='#016c59'
  )+
  theme(legend.position = "right")   
gg_nomega_proba

#######################
###    SAVE PLOT
#######################
setwd("~/Documents/REDACTION/manuscrit_radisic/6_CALIB_MOIST_PROF/figures/")
tikz(paste('gg_convergence_Nomega',Nomega,'_wrt_c_newcolors_cluster.tex', sep = ''), standAlone = TRUE, width=6.0, height=4.0)
gg_nomega
dev.off()
tools::texi2dvi(paste('gg_convergence_Nomega',Nomega,'_wrt_c_newcolors_cluster.tex', sep = ''),pdf=T)

setwd("~/Documents/REDACTION/manuscrit_radisic/6_CALIB_MOIST_PROF/figures/")
tikz(paste('gg_convergence_Nomega',Nomega,'_wrt_c_proba.tex', sep = ''), standAlone = TRUE, width=4.0, height=2.0)
gg_nomega_proba
dev.off()
tools::texi2dvi(paste('gg_convergence_Nomega',Nomega,'_wrt_c_proba.tex', sep = ''),pdf=T)

#pdf(file = paste('gg_convergence_Nomega',Nomega,'.pdf', sep = ''))
#plot(gg_nomega)
#dev.off()
#}

## save workspace
setwd("~/Documents/Cours_presentations_formations/4_Presentation_poster_article/MEXICO_CALIB_CHAP6/avancement_fs_6000/")
save.image()
