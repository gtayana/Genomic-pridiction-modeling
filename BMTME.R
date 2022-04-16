setwd("/gpfs/home/gtayana/BMORS2015/BMTME")
library(BMTME)
library(Hmisc)
library(GGally)
library(ggplot2)
library(car)
library(plyr)
library(agricolae)
library(BGLR)
library(caret)
library(tidyverse)
library(BMTME)
library(lme4)
library(pheatmap)
require(ggpubr)
require(Hmisc)
require(corrplot)
library(data.table) 
library(readr)
library(plotrix)
library(tidyr)
library(knitr)
library(kableExtra)
library(MplusAutomation)

#Phnotype data
#####################################
BLUP2015=read.csv("BLUP2015_spread.csv", header = TRUE)
head(BLUP2015)
BLUP2015=BLUP2015[,-c(1,7)]
#head(Dat2015)
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
pheno <- data.frame(GID = BLUP2015[, 1], Env = BLUP2015[, 2], Response = BLUP2015[, 3])
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



#BMTME: multi-trait and multi-location
pm_BMTME2015 <- BMTME(Y = Y, X = Z.E, Z1 = Z.G, Z2 = Z.EG, nIter = 5000, burnIn = 2500, thin = 2,bs = 50, testingSet = CrossValidation)

CV_BMTME2015=pm_BMTME2015$results
write.csv(CV_BMTME2015, "CV_BMTME_2015.csv")

#Saving summery of accuracy for all the trait in each Lociron
Accu_BMTME2015=summary(pm_BMTME2015)
write.csv(Accu_BMTME2015, "Acuuracy_BMTME_2015.csv")

boxplot(pm_BMTME2015, select ="Pearson", las = 2)
