
setwd("/gpfs/home/gtayana/BMORS2015/BME")

library(BMTME)
library(dplyr)
library(Hmisc)
library(GGally)
library(ggplot2)
library(car)
library(plyr)
library(agricolae)
library(BGLR)
library(coda)
library(caret)
library(tidyverse)
library(ggplot2)
library(BMTME)
library(lme4)
library(pheatmap)
require(ggpubr)
require(tidyverse)
require(Hmisc)
require(corrplot)
library(data.table) 
library(readr)
library(plotrix)
library(tidyr)
library(ggcorrplot)
library(lavaan)
library(semPlot)
library(OpenMx)
library(knitr)
library(kableExtra)
library(MplusAutomation)

# Geno data
######################################
Marker_dat=read.table("Transpose of Geno_GWAS2015.txt", header = TRUE, sep = "", row.names=1)
Marker_dat[1:5,1:5]

##change -1,0,1 to 0,1,2
Marker_dat <- Marker_dat + 1

#GBLUP relationship calculation 
tg=tcrossprod(scale(Marker_dat))  # transposed and cross product of the centered and scaled markers 
G=tg/mean(diag(tg)) # GBLUP calculated from <M>
G = as.matrix(G)
G[1:10,1:10]# visualize first 10 rows and columns

###############################################


#Import data
Dat2015=read.csv("BLUP2015_spread.csv", header = TRUE)
head(Dat2015)
Dat2015=Dat2015[,-c(2)]
head(Dat2015)

#order data
Y_NERF2015 = Dat2015[ which(Dat2015$Env=='NERF'), ]
Y_NERF2015 <- (Y_NERF2015[order(Y_NERF2015$GID),])

Y_SERF2015 = Dat2015[ which(Dat2015$Env=='SERF'), ]
Y_SERF2015 <- (Y_SERF2015[order(Y_SERF2015$GID),])

Y_Volga2015 = Dat2015[ which(Dat2015$Env=='Volga'), ]
Y_Volga2015 <- (Y_Volga2015[order(Y_Volga2015$GID),])

Y_Winner2015 = Dat2015[ which(Dat2015$Env=='Winner'), ]
Y_Winner2015 <- (Y_Winner2015[order(Y_Winner2015$GID),])

# Rename rows
rownames(Y_NERF2015)=1:nrow(Y_NERF2015)
rownames(Y_SERF2015)=1:nrow(Y_SERF2015)
rownames(Y_Volga2015)=1:nrow(Y_Volga2015)
rownames(Y_Winner2015)=1:nrow(Y_Winner2015)


#extract only trait data
Y1=as.matrix(Y_NERF2015[,-c(1,2)])
Y2=as.matrix(Y_SERF2015[,-c(1,2)])
Y3=as.matrix(Y_Volga2015[,-c(1,2,6)])
Y4=as.matrix(Y_Winner2015[,-c(1,2)])


#NERF
LG <- cholesky(G)#Geno data common for all the environment
########################
pheno <- data.frame(GID = Y_NERF2015[, 1], Response = Y_NERF2015[, 2])
CrossV <- CV.RandomPart(pheno, NPartitions = 10, PTesting = 0.2, set_seed = 123)
ZG <- model.matrix(~0 + as.factor(Y_NERF2015$GID))
Z.G <- ZG %*% LG
fm_BME_NERF2015 <- BME(Y = Y1, Z1 = Z.G, nIter = 15000, burnIn = 10000, thin = 2, bs = 50, testingSet = CrossV)

#Saving result from each N Partitions CV 
CV_BME_NERF2015=fm_BME_NERF2015$results
write.csv(CV_BME_NERF2015, "CV_NERF_BME_2015.csv")

#Saving summery of accuracy for all the trait in each location
Accu_BME_NERF2015=summary(fm_BME_NERF2015)
write.csv(Accu_BME_NERF2015, "Accuracy_NERF_BME_2015.csv")

###############################
#SERF
pheno <- data.frame(GID = Y_SERF2015[, 1], Response = Y_SERF2015[, 2])
CrossV <- CV.RandomPart(pheno, NPartitions = 10, PTesting = 0.2, set_seed = 123)
ZG <- model.matrix(~0 + as.factor(Y_SERF2015$GID))
Z.G <- ZG %*% LG
fm_BME_SERF2015 <- BME(Y = Y2, Z1 = Z.G, nIter = 15000, burnIn = 10000, thin = 2, bs = 50, testingSet = CrossV)

#Saving result from each N Partitions CV 
CV_BME_SERF2015=fm_BME_SERF2015$results
write.csv(CV_BME_SERF2015, "CV_SERF_BME_2015.csv")

#Saving summery of accuracy for all the trait in each location
Accu_BME_SERF2015=summary(fm_BME_SERF2015)
write.csv(Accu_BME_SERF2015, "Accuracy_SERF_BME_2015.csv")

############################
#Volga
pheno <- data.frame(GID = Y_Volga2015[, 1], Response = Y_Volga2015[, 2])
CrossV <- CV.RandomPart(pheno, NPartitions = 10, PTesting = 0.2, set_seed = 123)
ZG <- model.matrix(~0 + as.factor(Y_Volga2015$GID))
Z.G <- ZG %*% LG
fm_BME_Volga2015 <- BME(Y = Y3, Z1 = Z.G, nIter = 15000, burnIn = 10000, thin = 2, bs = 50, testingSet = CrossV)
#Saving result from each N Partitions CV 
CV_BME_Volga2015=fm_BME_Volga2015$results
write.csv(CV_BME_Volga2015, "CV_Volga_BME_2015.csv")

#Saving summery of accuracy for all the trait in each location
Accu_BME_Volga2015=summary(fm_BME_Volga2015)
write.csv(Accu_BME_Volga2015, "Accuracy_Volga_BME_2015.csv")

############################
#Winner
pheno <- data.frame(GID = Y_Winner2015[, 1], Response = Y_Winner2015[, 2])
CrossV <- CV.RandomPart(pheno, NPartitions = 10, PTesting = 0.2, set_seed = 123)
ZG <- model.matrix(~0 + as.factor(Y_Winner2015$GID))
Z.G <- ZG %*% LG
fm_BME_Winner2015 <- BME(Y = Y4, Z1 = Z.G, nIter = 15000, burnIn = 10000, thin = 2, bs = 50, testingSet = CrossV)
#Saving result from each N Partitions CV 
CV_BME_Winner2015=fm_BME_Winner2015$results
write.csv(CV_BME_Winner2015, "CV_Winner_BME_2015.csv")

#Saving summery of accuracy for all the trait in each location
Accu_BME_Winner2015=summary(fm_BME_Winner2015)
write.csv(Accu_BME_Winner2015, "Accuracy_Winner_BME_2015.csv")

#######################################
