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

LHS_train <- as.matrix(npyLoad("LHS_rep5_50_R500_dim6_sample_final.npy"))
LHS_test <- as.matrix(npyLoad("LHS25_100_R500_dim6_sample_final.npy"))
LHS_enriched_test <- as.matrix(npyLoad("LHS110_200_R300_dim6_sample_final.npy"))

indices_6_parms <- c(63, 65, 68, 71, 72, 99)

plotting_data_train <- data.frame(cbind(LHS_train[1:50,indices_6_parms],rep('train', 50)))
plotting_data_test <- data.frame(cbind(LHS_test[1:100,indices_6_parms],rep('test', 100)))
plotting_data_enriched_test <- data.frame(cbind(LHS_enriched_test[1:110,indices_6_parms],rep('enriched_test', 110)))


df <- data.frame(rbind(plotting_data_train, plotting_data_test, plotting_data_enriched_test))
df <- data.frame(as.numeric(df[,1]),as.numeric(df[,2]),as.numeric(df[,3]),
                      as.numeric(df[,4]),as.numeric(df[,5]),as.numeric(df[,6]),'type' = df[,7])
head(df)

## Pair plot for raw train and test set.
ggpairs(df,aes(pch = type,color = type, alpha = 0.5),
        lower = list(continuous = wrap("points", alpha = 0.5, size=1, pch = 1), 
                     combo = wrap("dot", alpha = 0.4, size=0.2) ),
        upper = list(continuous = wrap("points", alpha = 0.5, size=1, pch = 1), 
                     combo = wrap("dot", alpha = 0.4, size=0.2) )) +theme_bw() +
  theme(legend.position = 'bottom') 

## Param. 1
df_test <- data.frame(as.numeric(df[,1]),df[,7])
colnames(df_test) <- c('value', 'type')
ggplot(data = df_test, aes(x=value, y = ..density.., color = type, fill = type)) +
  geom_histogram(alpha=0.5, position="identity") + theme_bw()

## Param. 2
df_test <- data.frame(as.numeric(df[,2]),df[,7])
colnames(df_test) <- c('value', 'type')
ggplot(data = df_test, aes(x=value, y = ..density.., color = type, fill = type)) +
  geom_histogram(alpha=0.5, position="identity") + theme_bw()

## Param. 3
df_test <- data.frame(as.numeric(df[,3]),df[,7])
colnames(df_test) <- c('value', 'type')
ggplot(data = df_test, aes(x=value, y = ..density.., color = type, fill = type)) +
  geom_histogram(alpha=0.5, position="identity") + theme_bw()

## Param. 4
df_test <- data.frame(as.numeric(df[,4]),df[,7])
colnames(df_test) <- c('value', 'type')
ggplot(data = df_test, aes(x=value, y = ..density.., color = type, fill = type)) +
  geom_histogram(alpha=0.5, position="identity") + theme_bw()

## Param. 5
df_test <- data.frame(as.numeric(df[,5]),df[,7])
colnames(df_test) <- c('value', 'type')
ggplot(data = df_test, aes(x=value, y = ..density.., color = type, fill = type)) +
  geom_histogram(alpha=0.5, position="identity") + theme_bw()

## Param. 6
df_test <- data.frame(as.numeric(df[,5]),df[,7])
colnames(df_test) <- c('value', 'type')
ggplot(data = df_test, aes(x=value, y = ..density.., color = type, fill = type)) +
  geom_histogram(alpha=0.5, position="identity") + theme_bw()

