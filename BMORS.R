setwd("/gpfs/home/gtayana/BMORS2015/BMORS")
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

#Phnotype data
#####################################
BLUP2015=read.csv("BLUP2015_spread.csv", header = TRUE)
head(BLUP2015)
BLUP2015=BLUP2015[,-c(1)]
head(Dat2015)
# Data transformation in to factor or character
BLUP2015 = transform (BLUP2015, Env = factor(Env), Line = as.character(Line))
str(BLUP2015)

#Important step and never to miss it: ordering according to "lines" and "Env" first
BLUP2015<-(BLUP2015[order(BLUP2015$Env,BLUP2015$Line),])
rownames(BLUP2015)=1:nrow(BLUP2015)
head(BLUP2015)
str(BLUP2015)

# extracting and make matrix of only phenotype data
Y <- as.matrix(BLUP2015[, -c(1,2)])
########################################


#Cross validation
#Data preparation for cross validation
pheno <- data.frame(GID = BLUP2015[, 1], Env = BLUP2015[, 2], Response = BLUP2015[, 6])
#Random k partition
CrossValidation <-CV.RandomPart(pheno, NPartitions = 10, PTesting = 0.2, set_seed = 123)

#Cross validation to use K-fold instead of random K partition
#CrossValidation <- CV.KFold(pheno, DataSetID = "GID", K = 10, set_seed = 123)
#############################################


#Genotype data
Marker_dat=read.table("Transpose of Geno_GWAS2015.txt", header = TRUE, sep = "", row.names=1)
Marker_dat[1:5,1:5]

##change -1,0,1 to 0,1,2
Marker_dat <- Marker_dat + 1

#GBLUP relationship calculation 
tg=tcrossprod(scale(Marker_dat))  # transposed and cross product of the centered and scaled markers 
G=tg/mean(diag(tg)) # GBLUP calculated from <M>
G = as.matrix(G)
G[1:10,1:10]# visualize first 10 rows and columns
##############################################


#Design matrix
LG <- cholesky(G)#Geno data
ZG <- model.matrix(~0 + as.factor(BLUP2015$Line))
Z.G <- ZG %*% LG
Z.E <- model.matrix(~0 + as.factor(BLUP2015$Env))
ZEG <- model.matrix(~0 + as.factor(BLUP2015$Line):as.factor(BLUP2015$Env))
G2 <- kronecker(diag(length(unique(BLUP2015$Env))), data.matrix(G))
LG2 <- cholesky(G2)
Z.EG <- ZEG %*% LG2
#################################################



#BMORS multi-trait and multi-location uding BRR
ETA <- list(Loc = list(X = Z.E, model = "BRR"), Gen = list(X = Z.G, model = "BRR"), LocGen = list(X = Z.EG, model = "BRR"))

#BMORS model
pm_BRR2015 <- BMORS(Y, ETA = ETA, nIter = 15000, burnIn = 10000, thin = 2, progressBar = TRUE, testingSet = CrossValidation, digits = 4)

#Saving result from each N Partitions CV 
CV_BRR2015=pm_BRR2015$results
write.csv(CV_BRR2015, "CV_BRR_BMORS_2015.csv")


#Saving summery of accuracy for all the trait in each Lociron
Accu_BRR2015=summary(pm_BRR2015)
write.csv(Accu_BRR2015, "Accuracy_BRR_BMORS_2015.csv")

# Box plot of the 2015 ac curacies 
boxplot(pm_BRR2015, select ="Pearson", las = 2)
###############################################


#BMORS multi-trait and multi-location using BayesB

ETA <- list(Loc = list(X = Z.E, model = "BayesB"), Gen = list(X = Z.G, model = "BayesB"), LocGen = list(X = Z.EG, model = "BayesB"))
#BMORS model
pm_BayesB <- BMORS(Y, ETA = ETA, nIter = 15000, burnIn = 10000, thin = 2, progressBar = TRUE, testingSet = CrossValidation, digits = 4)

#Saving result from each N Partitions CV 
CV_BayesB2015=pm_BayesB$results
write.csv(CV_BayesB2015, "CV_BayesB_BMORS_2015.csv")

#Saving summery of accuracy for all the trait in each Lociron
Accu_BayesB2015=summary(pm_BayesB)
write.csv(Accu_BayesB2015, "Accuracy_BayesB_BMORS_2015.csv")

# Box plot of the 2015 ac curacies 
boxplot(pm_BayesB, select ="Pearson", las = 2)



#BMORS multi-trait and multi-location using BayesA

ETA <- list(Loc = list(X = Z.E, model = "BayesA"), Gen = list(X = Z.G, model = "BayesA"), LocGen = list(X = Z.EG, model = "BayesA"))
#BMORS model
pm_BayesA <- BMORS(Y, ETA = ETA, nIter = 15000, burnIn = 10000, thin = 2, progressBar = TRUE, testingSet = CrossValidation, digits = 4)

#Saving result from each N Partitions CV 
CV_BayesA2015=pm_BayesA$results
write.csv(CV_BayesA2015, "CV_BayesA_BMORS_2015.csv")

#Saving summery of accuracy for all the trait in each Lociron
Accu_BayesA2015=summary(pm_BayesA)
write.csv(Accu_BayesA2015, "Accuracy_BayesA_BMORS_2015.csv")

# Box plot of the 2015 ac curacies 
boxplot(pm_BayesA, select ="Pearson", las = 2)


