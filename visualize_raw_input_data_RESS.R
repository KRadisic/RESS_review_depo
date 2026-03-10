########################################################
### VISUALIZATION INPUT RAW DATA RESS
########################################################
rm(list=ls())

library("RcppCNPy") # needed to read .npy
library("ggplot2")
library("GGally")
?RcppCNPy

data_directory <- "~/Documents/RESS_review_depo/"

setwd(data_directory)

## Load input PESHMELBA data
Nx_train <- 50
Nx_test <- 100
LHS_train <- as.matrix(npyLoad("LHS_rep5_50_R500_dim6_sample_final.npy"))
LHS_test <- as.matrix(npyLoad("LHS25_100_R500_dim6_sample_final.npy"))
#LHS_enriched_test <- as.matrix(npyLoad("LHS110_200_R300_dim6_sample_final.npy"))

## Indices of parameters used in metamodeling
indices_6_parms <- c(63, 65, 68, 71, 72, 99)
names <- read.csv('entry_factor_names.csv', header = FALSE)
names$V1[indices_6_parms]

## Arrange the data into a dataframe for ggplot plotting
plotting_data_train <- data.frame(cbind(LHS_train[1:Nx_train,indices_6_parms],rep('train', Nx_train)))
plotting_data_test <- data.frame(cbind(LHS_test[1:Nx_test,indices_6_parms],rep('test', Nx_test)))
#plotting_data_enriched_test <- data.frame(cbind(LHS_enriched_test[1:110,indices_6_parms],rep('enriched_test', 110)))

df <- data.frame(rbind(plotting_data_train, plotting_data_test))#, plotting_data_enriched_test))
df <- data.frame(as.numeric(df[,1]),as.numeric(df[,2]),as.numeric(df[,3]),
                      as.numeric(df[,4]),as.numeric(df[,5]),as.numeric(df[,6]),'type' = df[,7])
colnames(df) <- c(names$V1[indices_6_parms], 'type')
df$type <- as.factor(df$type) 

## Pair plot for raw train and test set. #c('train' = "#7570b3", 'test' = "#1b9e77") 
ggpairs(df, aes(color = type, alpha = 0.5), columns = c(1:6), legend = 2,
        lower = list(continuous = wrap("points", alpha = 0.5, size=1, pch = 1), 
                     combo = wrap("dot", alpha = 0.4, size=0.2) ),
        upper = list(continuous = wrap("points", alpha = 0.5, size=1, pch = 1), 
                     combo = wrap("dot", alpha = 0.4, size=0.2) )) +theme_bw() +
        scale_color_manual(breaks = c("test", "train"), values = c("#1b9e77", "#7570b3")) + # density plots
        scale_fill_manual(breaks = c("test", "train"), values = c("#1b9e77", "#7570b3")) + # diagonal fill
        theme(legend.position = "bottom")

