library('BGLR')

setwd('//gpfs/home/gtayana/GG_SPARSE')

gen = read.table("SDSU2017_pheno_geno.txt", header = TRUE, row.names = 1)
length(gen)
phe2017=gen[ ,1:69]
phe2017[,2]

phe2<-scale(phe2017, scale=TRUE, center = TRUE)


g2<-scale(gen[,70:ncol(gen)], scale=TRUE, center = TRUE)
G<-tcrossprod(g2)/ncol(g2)
G[1:10, 1:10]

EVD <-eigen(G)
EVD$vectors

phe3=phe2[,1:4]

# rm(list=ls())



COR1 <- matrix(nrow=1, ncol=6,NA)
colnames(COR1) <- c('Pairs','Fold','Part','SingleEnv', 'AcrossEnv', 'MxE')
ENV_pairs=c(list(c(1,2),c(1,2,3), c(1,2,3,4)))
for(p in 1:length(ENV_pairs)){
      for (m in c(0.2, 0.3, 0.4, 0.5)){
        for (k in 1:2){
          
      env <-ENV_pairs[[p]]# choose any set of environments from 1:ncol(Y)
      nEnv <- length(env)
      Y <- phe3[,env]
      n <- nrow(Y)
      percTST<-m
      nTST <- round(percTST*n)
      set.seed(k)
      tst<-sample(1:n,size=nTST,replace=FALSE)
      YNA <- Y
      YNA[tst,]<-NA
      
      
      ## Single environments models #####################################
      YHatSE <- matrix(nrow=nrow(Y),ncol=ncol(Y),NA)
      ETA <- list(G=list(K=G,model='RKHS'))
      for(i in 1:nEnv){
        prefix <- paste(colnames(Y)[i],"_",sep="")
        fm <-BGLR(y=YNA[,i],ETA=ETA,nIter=15000,burnIn=10000,saveAt=prefix)
        YHatSE[,i] <- fm$yHat
      }
      ## Across environment model (ignoring GxE) #######################
      yNA <- as.vector(YNA)
      # Fixed effect (env-intercepts)
      envID <- rep(env,each=nrow(Y))
      ETA <- list(list(~factor(envID)-1,model="FIXED"))
      # Main effects of markers
      G0 <- kronecker(matrix(nrow=nEnv,ncol=nEnv,1),G)
      ETA[[2]] <- list(K=G0,model='RKHS')
      # Model Fitting
      prefix <- paste(c('Across',colnames(Y),''),collapse='_')
      fm <- BGLR(y=yNA,ETA=ETA,nIter=15000,burnIn=10000,saveAt=prefix)
      YHatAcross <- matrix(fm$yHat,ncol=nEnv)
      ## MxE Interaction Model #########################################
      # Adding interaction terms
      for(i in 1:nEnv){
        tmp <- rep(0,nEnv) ; tmp[i] <- 1; G1 <- kronecker(diag(tmp),G)
        ETA[[(i+2)]] <- list(K=G1,model='RKHS')
      }
      # Model Fitting
      prefix <- paste(c('MxE',colnames(Y),''),collapse='_')
      fm <- BGLR(y=yNA,ETA=ETA,nIter=15000,burnIn=10000,saveAt=prefix)
      YHatInt <- matrix(fm$yHat,ncol=nEnv)
      
      #COR1 <- matrix(nrow=length(env),ncol=5,NA)
      #colnames(COR1) <- c('Fold','Part','SingleEnv', 'AcrossEnv', 'MxE')
      
      # CV1 Creating a Testing Sets for CV1
      COR <- matrix(nrow=length(env),ncol=6,NA)
      colnames(COR) <- c('Pairs','Fold','Part','SingleEnv', 'AcrossEnv', 'MxE')
      rownames(COR) <- colnames(Y)
        for(i in 1:nEnv){
        tst <- which(is.na(YNA[,i]))
        COR[i,1]= as.numeric(paste(ENV_pairs[[p]],collapse=""))
        COR[i,2]=m
        COR[i,3]=k
        COR[i,4] <- cor(Y[tst,i],YHatSE[tst,i])
        COR[i,5] <- cor(Y[tst,i],YHatAcross[tst,i])
        COR[i,6] <- cor(Y[tst,i],YHatInt[tst,i])
        
      }
      COR1=rbind(COR, COR1)
    }
  }
}
write.csv(COR1,'COR_CV1.csv' )
