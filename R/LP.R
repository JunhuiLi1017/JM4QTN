LP <- function(vecPheno,PhenoData,GenoData_EST,GenoData_QTL,vecThrVal,CChr,Interval,nPB=2000,alpha=0.05){
  if(nrow(PhenoData) != nrow(GenoData_EST)-2 && nrow(PhenoData) != nrow(GenoData_QTL)-2){
    stop("Sample sieze must be equal in PhenoData, GenoData_EST and GenoData_QTL\n")
  }
  if(!all(vecPheno %in% colnames(PhenoData))){
    stop("Trait vector must be included in PhenoData\n")
  }
  if(!all(colnames(PhenoData)[1:2] == c("Indi","Popu"))){
    stop("Top four of column in PhenoData must be Indi\tPopu\tMA\tPA\n")
  }
  if(alpha <=0 || alpha >=1){
    stop("alpha should be 0~1\n")
  }
  mapORpos <- matrix(NA,ncol(GenoData_EST),3)
  mapORpos <- as.data.frame(mapORpos)
  colnames(mapORpos) <- c("marker","chr","pos")
  mapORpos[,1] <- colnames(GenoData_EST)
  mapORpos[,2:3] <- t(GenoData_EST[1:2,])
  vmapORpos <- matrix(NA,ncol(GenoData_QTL),3)
  vmapORpos <- as.data.frame(vmapORpos)
  colnames(vmapORpos) <- c("marker","chr","pos")
  vmapORpos[,1] <- colnames(GenoData_QTL)
  vmapORpos[,2:3] <- t(GenoData_QTL[1:2,])
  
  GenoData_EST <- GenoData_EST[-c(1:2),]
  GenoData_QTL <- GenoData_QTL[-c(1:2),]
  
  nT <- length(vecPheno)                         #nT: No. of Traits
  NoChr <- nlevels(as.factor(mapORpos$chr))
  nM_V <- nrow(vmapORpos)                         #nM_V: No. of Markers 
  CChMK <- subset(mapORpos,mapORpos$chr==CChr)          #Marker on Current Chromosome
  LM <- min(Interval)
  RM <- max(Interval)
  if(CChMK[nrow(CChMK),3]<RM){
    stop("Maximum of confidence interval exceeded end of current chromosome, please check out your chr and confidence intercal")
  }
  if(CChMK[1,3]>LM){
    stop("Minimum of confidence interval exceeded start of current chromosome, please check out your chr and confidence intercal")
  }
  NumLM <- which(CChMK$pos - LM <= 0)
  NumRM <- which(CChMK$pos - RM >= 0)
  MinLM <- CChMK[NumLM[length(NumLM)],"marker"]#The Max No of Marker on the left side			
  MaxRM <- CChMK[NumRM[1],"marker"]		#The Min No of Marker on the right side	
  vChMK <- subset(vmapORpos,vmapORpos$chr==CChr)
  WhichMinLM <- which(vChMK[,1] %in% MinLM)
  WhichMaxRM <- which(vChMK[,1] %in% MaxRM)
  #* read all marker information in Confidence Interval by TMData-----------------------------
  NoLM <- which(colnames(GenoData_QTL) %in% MinLM)
  NoRM <- which(colnames(GenoData_QTL) %in% MaxRM)
  GenoData_SI <- GenoData_QTL[,NoLM:NoRM]
  CI_V <- vChMK[WhichMinLM:WhichMaxRM,3]
  names(CI_V) <- vChMK[WhichMinLM:WhichMaxRM,1]
  ###Cofactors selection using MMR function
  PhenoData[,"Popu"] <- as.numeric(PhenoData[,"Popu"])
  TRMData <- cbind(PhenoData[,c(1:4,which(colnames(PhenoData) %in% vecPheno))],GenoData_EST)
  nT <- length(vecPheno)
  nObs <- nrow(TRMData)
  nRP <- nlevels(as.factor(TRMData$Popu)) 
  Popu <- data.matrix(TRMData["Popu"])
  
  MulTestTab <- matrix(NA,nT,12)
  MulTestTab <- as.data.frame(MulTestTab)
  colnames(MulTestTab) <- c("LOD_S","lod_s","pos_s","LOD_P","lod_p","pos_p","LB_P","RB_P","LOD_L","lod_l","LOD_PL","lod_pl")
  rownames(MulTestTab) <- vecPheno
  Delta <- 0.0001   #Delta: a small value to separate CF-marker with its adj-marker
  # if it equal 0, then it means CF-marker affects its adj-marker
  ### Cofactor selection for multiple trait model
  selection <- 'stepwise'
  tolerance <- 1e-7
  Trace <- "Pillai"
  select<- "SL"
  sls <- vecThrVal[1]
  sle <- sls
  Choose <- "SBC"  
  if(nRP > 1){
    include0 <- "Popu"
    remove.cf <- c(1,2)
  }else{
    include0 <- NULL
    remove.cf <- 1
  }
  notX <- c(1,3,4)
  Class0 <- include0
  tempCF_P <- stepwise(TRMData, vecPheno, notX, include0, Class0, selection, select, sle, sls, tolerance, Trace, Choose)$variate
  tempCF_P <- tempCF_P[-remove.cf]
  if(length(tempCF_P)==0){
    warning("There is no cofactor for trait\t",vecPheno,"\twith stepwise regression for significant level\t",vecThrVal[1],"\n")
  }
  #:Complete cofactors detail Infor. after getting their name_vector, tempCF_P
  CF <- matrix(nrow=0,ncol=4)
  CF <- as.data.frame(CF)
  if(length(tempCF_P) > 0){
    for(k in 1:length(tempCF_P)){
      for(i in 1:NoChr){
        ChrMarkerID <- as.character(subset(mapORpos,mapORpos$chr==i)[,1])
        ChrMarkerPos <- subset(mapORpos,mapORpos$chr==i)[,3]
        if(tempCF_P[k] %in% ChrMarkerID){
          Cchr <- i                   #Current Chr
          Mseq <- which(ChrMarkerID %in% tempCF_P[k])
          pos <- ChrMarkerPos[Mseq]     #current CF position
          if(Mseq == 1){
            Lpos <- 0          #Lpos, Rpos : left and right position (effect region)
            Rpos <- ChrMarkerPos[Mseq+1] - Delta
          }else if(Mseq == length(ChrMarkerPos)){
            Lpos <- ChrMarkerPos[Mseq-1] + Delta
            Rpos <- ChrMarkerPos[Mseq]
          }else{
            Lpos <- ChrMarkerPos[Mseq-1] + Delta
            Rpos <- ChrMarkerPos[Mseq+1] - Delta
          }
          CF <- rbind(CF,c(Cchr,pos,Lpos,Rpos))
        }
      }#i
    }#k
    rownames(CF) <- tempCF_P
    colnames(CF) <- c("chr","pos","Lpos","Rpos")
  }
  ##2.2: QTL scanning(1D scan both for pleiotropic and single trait to every marker locus)
  nM_V <- ncol(GenoData_SI)
  Ftest_P <- matrix(nrow=1,ncol=nM_V)          #F scores of each SNP;
  Ftest_P <- as.data.frame(Ftest_P)
  colnames(Ftest_P) <- colnames(GenoData_SI)
  PrFt_P <- Ftest_P                           #Pr(>F) of each SNP;
  detRSSF_P <- Ftest_P
  LOD_M <- matrix(nrow=nT+1,ncol=nM_V)	 #LOD_M:Single Effect Statistic in Pleiotropic model;
  LOD_M <- as.data.frame(LOD_M)
  colnames(LOD_M) <- colnames(GenoData_SI)
  rownames(LOD_M) <- c(vecPheno,"Joint")
  Y <- data.matrix(TRMData[vecPheno])           #Y matrix
  #q=15
  for (q in 1:nM_V) {
    SNPq <- colnames(GenoData_SI)[q]
    ##make two models for every SNP
    #*cofactors reducing
    if(nrow(CF) > 0){
      CCFvec <- tempCF_P	#CCFvec: Current CF_vector
      if (CChr %in% CF[,1]) {
        Cpos <- CI_V[q]
        CchrCF <- subset(CF,CF[,1] %in% CChr)
        CharC <- NULL
        for(ccc in 1:nrow(CchrCF)){
          # ccc=1
          if((Cpos >= CchrCF[ccc,3]) && (Cpos <= CchrCF[ccc,4])){
            CharC <- rownames(CchrCF)[ccc]
            CCFvec <- setdiff(CCFvec,CharC)
          }
        }
      }
    }else{
      CCFvec <- NULL
    }
    #Linear Model for different Situations: 1st. Multiple Crossed Lines(MCL) and 2nd. Single Crossed Lines(SCL).
    #1st MCL
    if(length(CCFvec) == 0){
      Xf <- GenoData_SI[,SNPq]
    }else{
      Xr <- TRMData[c(CCFvec)]
      Xf <- cbind(Xr,GenoData_SI[,SNPq])
    }
    Xr <- data.matrix(Xr)
    Xf <- data.matrix(Xf)
    if(nRP>1){
      if(length(CCFvec) ==0 ){
        lmr <- lm(Y~Popu)
        lmf <- lm(Y~Popu+Xf:Popu)
      }else{
        lmr <- lm(Y~Popu+Xr:Popu)
        lmf <- lm(Y~Popu+Xf:Popu)
      }
      #2nd SCL
    }else{
      if(length(CCFvec) ==0 ){
        lmr <- lm(Y~1)
        lmf <- lm(Y~1+Xf)
      }else{
        lmr <- lm(Y~1+Xr)
        lmf <- lm(Y~1+Xf)
      }
    }
    # get RSSp and RSSr and then pleiotropic model for Joint and single effect		
    resF_P <- resid(lmf)
    resR_P <- resid(lmr)
    RSSF_P <- t(resF_P) %*% resF_P
    RSSR_P <- t(resR_P) %*% resR_P
    detRSSF <- det(RSSF_P)
    detRSSR <- det(RSSR_P)
    RSSRs <- diag(RSSR_P)
    RSSFs <- diag(RSSF_P)
    #residual df and Fvalue+Pvalue+JEs+SEs
    resdff <- df.residual(lmf)
    resdfr <- df.residual(lmr)
    rDFd <- resdfr - resdff              #difference between rDF2 and rDF2
    #Fvalue <- (detRSSR-detRSSF)/(resdfr-resdff)/(detRSSF/resdff)
    #PrFt_P[q] <- 1-pf(Fvalue,resdfr-resdff,resdff)
    #anvar <- anova(lmr,lmf,test="Wilks")
    anvar <- anova(lmr,lmf)
    Ftest_P[SNPq] <- anvar[2,"approx F"]
    PrFt_P[SNPq] <- anvar[2,"Pr(>F)"]
    detRSSF_P[SNPq] <- detRSSF
    
    LOD_M[nT+1,SNPq] <- 0.5*nObs*log10(detRSSR/detRSSF)
    #Wilks=detRSSF/detRSSR
    #Effect for each trait
    #for(nt in 1:nT){
    #nt=1
    LOD_M[1:nT,SNPq] <- 0.5*nObs*log10(RSSRs/RSSFs)
    #}
  }#q
  ##2.3: F-test for single trait with pleiotropic QTL on alpha=0.1 levels
  CriLOD_M <- vecThrVal[-1]		#PermuTest, alpha=0.1 level
  MulTestTab[1:nT,"LOD_S"] <- as.numeric(CriLOD_M[1:nT])
  MulTestTab[1:nT,"LOD_P"] <- CriLOD_M[nT+1]
  MulTestTab[1:nT,"LOD_L"] <- CriLOD_M[nT+2]
  for(nt in 1:nT){
    MulTestTab[nt,"lod_s"] <- max(LOD_M[nt,])
    MulTestTab[nt,"pos_s"] <- CI_V[which(LOD_M[nt,] %in% max(LOD_M[nt,]))]
  }
  if(!all(MulTestTab[,1] < MulTestTab[,2])){
    return(MulTestTab)
    stop("Single trait QTLs are not all significant")
  }
  MulTestTab[1,"lod_p"] <- max(LOD_M[nT+1,])
  if(MulTestTab[1,"lod_p"] >= MulTestTab[1,"LOD_P"]){
    CPName <- colnames(LOD_M)[LOD_M[nT+1,] %in% max(LOD_M[nT+1,])]
    CID <- 1 #Conf.Inter.Degrade
    #CriQRD <- 5.0	 #Criterion of QTL redundant distance(in cM)
    tempQTL_P <- NULL
    QnameVec_P <- NULL
    
    CNr <- which(LOD_M[nT+1,] %in% max(LOD_M[nT+1,]))
    Cpos <- as.vector(CI_V[CNr])
    Lflag <- 0	 #left flag of	max(LOD_M[nT+1,])
    Rflag <- 0	 #right flag of	max(LOD_M[nT+1,])
    #check left side, whether it is <	max(LOD_M[nT+1,])
    if(CNr > 1){
      if( (is.na(LOD_M[nT+1,CNr-1])) | (!(is.na(LOD_M[nT+1,CNr-1])) && (LOD_M[nT+1,CNr-1] < max(LOD_M[nT+1,]))) ){
        Lflag <- 1
      }
    }else{
      Lflag <- 1
    }
    #check right side, whether it is <	max(LOD_M[nT+1,])
    if(CNr < nM_V){
      if( (is.na(LOD_M[nT+1,CNr+1])) | (!(is.na(LOD_M[nT+1,CNr+1])) && (LOD_M[nT+1,CNr+1] < max(LOD_M[nT+1,]))) ){
        Rflag <- 1
      }
    }else{
      Rflag <- 1
    }
    if((Lflag == 1) && (Rflag == 1)){
      #*search left-bound-of-QTL-CI
      LNr <- CNr
      if(CNr > 1){
        while((LNr-1 > 0) && ((is.na(LOD_M[nT+1,LNr-1]))|(!(is.na(LOD_M[nT+1,LNr-1]))&&(LOD_M[nT+1,LNr-1] >	max(LOD_M[nT+1,])-CID)))){
          LNr <- LNr-1
        }
      }
      if(LNr == 1){
        LBpos <- round(as.vector(CI_V[LNr]),digits=3)
      }else{
        LBpos <- round(mean(as.vector(CI_V[c(LNr-1,LNr)])),digits=3)
      }
      #*search right-bound-of-QTL-CI
      RNr <- CNr
      if(CNr < nM_V){
        while((RNr+1 < nM_V) && ((is.na(LOD_M[nT+1,RNr+1]))|(!(is.na(LOD_M[nT+1,RNr+1]))&&(LOD_M[nT+1,RNr+1] > max(LOD_M[nT+1,])-CID)))){
          RNr <- RNr+1
        }
      }
      if(RNr == nM_V){
        RBpos <- round(as.vector(CI_V[RNr]),digits=3)
      }else{
        RBpos <- round(mean(as.vector(CI_V[c(RNr,RNr+1)])),digits=3)
      }
    }
    tempQTL_P <- c( max(LOD_M[nT+1,]),CChr,Cpos,LBpos,RBpos) # max(LOD_M[nT+1,]) is subset of JEs_P,data.frame
    tempQTL_P <- as.data.frame(t(tempQTL_P))
    QnameVec_P <- CPName
    
    if(length(QnameVec_P) > 0){
      rownames(tempQTL_P) <- QnameVec_P
      colnames(tempQTL_P) <- c("LOD_Joint","chr","pos","LBpos","RBpos")
    }
    ##Single-trait model test against alpha=0.1 level thresholds
  }else{
    tempQTL_P <- NULL
    stop("Pleiotropic Model is not significant!")
  }
  
  if(length(tempQTL_P)!=0){
    MulTestTab[1,"lod_p"] <- tempQTL_P[1,"LOD_Joint"]
    MulTestTab[1,"pos_p"] <- tempQTL_P[1,"pos"]
    MulTestTab[1:2,"LB_P"] <- tempQTL_P[1:2,"LBpos"]
    MulTestTab[1:2,"RB_P"] <- tempQTL_P[1:2,"RBpos"]
  }
  
  ##2.4: construct linkage model
  if(length(tempQTL_P)>0){
    hQTL_S <- NULL
    
    for(nt in 1:nT){
      if(max(LOD_M[nt,])>CriLOD_M[nt]){
        hQTL_S <- append(hQTL_S,names(which.max(LOD_M[nt,])))
      }
    }#nt 
    
    if(length(hQTL_S)>0){
      #pleiotropic model
      PorL <- NULL
      #CCFvec2 <- tempCF_P
      resF_L <- matrix(NA,nObs,nT)
      resR_L <- matrix(NA,nObs,nT)
      for(nt in 1:nT){
        #nt=1
        CCFvec2 <- tempCF_P
        CF_P2 <- CF[CCFvec2,]
        #reducing cofactors
        SNPs <- hQTL_S[nt]
        Cpos <- CI_V[SNPs]
        CchrCF <- subset(CF_P2,CF_P2[,1]==CChr)
        ChrC <- NULL
        if(nrow(CchrCF)>0){
          for(nCf in 1:nrow(CchrCF)){
            if(Cpos >= CchrCF[nCf,"Lpos"] && Cpos <= CchrCF[nCf,"Rpos"]){
              ChrC <- rownames(CchrCF)[nCf]
              CCFvec2 <- setdiff(CCFvec2,ChrC)
            }
          }#ncf
        }
        #}#nt
        #for(nt in 1:nT){
        #nt=1
        #run 2 models, reduced model and full model
        SNPs <- hQTL_S[nt]
        if(nRP>1){
          if(length(CCFvec2)==0){
            Xf1 <- data.matrix(GenoData_SI[c(SNPs)])
            lmr_L <- lm(Y[,nt]~Popu)
            lmf_L <- lm(Y[,nt] ~ Popu+Xf1:Popu)
            resR_L[,nt] <- data.matrix(resid(lmr_L))
            resF_L[,nt] <- data.matrix(resid(lmf_L))
          }else{
            Xr <- TRMData[c(CCFvec2)]
            Xf1 <- cbind(Xr,GenoData_SI[c(SNPs)])
            Xr <- data.matrix(Xr)
            Xf1 <- data.matrix(Xf1)
            lmr_L <- lm(Y[,nt]~Popu+Xr:Popu)
            lmf_L <- lm(Y[,nt] ~ Popu+Xf1:Popu)
            resR_L[,nt] <- data.matrix(resid(lmr_L))
            resF_L[,nt] <- data.matrix(resid(lmf_L))
          }
        }else{
          if(length(CCFvec2)==0){
            Xf1 <- data.matrix(GenoData_SI[,SNPs])
            lmr_L <- lm(Y[,nt] ~ 1)
            lmf_L <- lm(Y[,nt] ~ 1+Xf1)
            resR_L[,nt] <- data.matrix(resid(lmr_L))
            resF_L[,nt] <- data.matrix(resid(lmf_L))
          }else{
            Xr <- TRMData[c(CCFvec2)]
            Xf1 <- cbind(Xr,GenoData_SI[,SNPs])
            Xr <- data.matrix(Xr)
            Xf1 <- data.matrix(Xf1)
            lmr_L <- lm(Y[,nt] ~ 1+Xr)
            lmf_L <- lm(Y[,nt] ~ 1+Xf1)
            resR_L[,nt] <- data.matrix(resid(lmr_L))
            resF_L[,nt] <- data.matrix(resid(lmf_L))
          } 
        }
      }#nt
      #reduced model
      #resR_L <- resid(lmr_L)
      RSSR_L <- t(resR_L) %*% resR_L
      detRSSR_L <- det(RSSR_L)
      #full model	
      RSSF_L <- t(resF_L) %*% resF_L
      detRSSF_L <- det(RSSF_L)
      #Linkage model vs Reduced model
      #1. Approximate F test 
      
      #2. Likelihood ratio test
      dfR <- df.residual(lmr_L)
      dfF <- df.residual(lmf_L)
      #MulTestTab[1,"lod_l"] <- (df.residual(lmr)-1 - 0.5*(nT-(dfR-dfF)))*log(detRSSR_L/detRSSF_L)
      MulTestTab[1,"lod_l"] <- 0.5*nObs*log10(detRSSR_L/detRSSF_L)
      if(MulTestTab[1,"lod_l"] < MulTestTab[1,"LOD_L"]){
        return(MulTestTab)
        stop("Linakge Model is not significant!")
      }
      #likelihood ratio test(Pleiotropy VS. Linkage)
      if(length(unique(hQTL_S)) == 1){
        PorL <- 0.5*nObs*log10(as.numeric(detRSSF_P[unique(hQTL_S)])/detRSSF_L)
      }else{
        PorL <- 0.5*nObs*log10(as.numeric(detRSSF_P[QnameVec_P])/detRSSF_L)
      }
    }
    ###Parametric Bootstrap-----------------------------------------------------------
    ## information about pleiotropic QTL (CF and Pleiotropic QTL)
    #*reducing cofactors
    CCFvec3 <- tempCF_P
    SNPp <- QnameVec_P
    #reducing cofactors
    Cpos <- CI_V[SNPp]
    CchrCF <- subset(CF,CF[,1]==CChr)
    ChrC <- NULL
    if(nrow(CchrCF)>0){
      for(nCf in 1:nrow(CchrCF)){
        #nCf=1
        if(Cpos >= CchrCF[nCf,"Lpos"] && Cpos <= CchrCF[nCf,"Rpos"]){
          ChrC <- rownames(CchrCF)[nCf]
          CCFvec3	<- CCFvec3[!CCFvec3 %in% ChrC]
        }
      }#nCf
    }
    ###Estimate effect of Pleiotropic QTL
    Xcp <- cbind(TRMData[c(CCFvec3)],GenoData_SI[,QnameVec_P])
    Xcp <- data.matrix(Xcp)
    if(nRP>1)	lmp <- lm(Y~Popu+Xcp:Popu) else lmp <- lm(Y~1+Xcp)
    #resid(lmp)
    sigma.sqrt <- sqrt(deviance(lmp)/df.residual(lmp))
    Y.fit <- fitted(lmp)
    PorL_PB <- rep(NA,nPB)
    for(npb in 1:nPB){
      #simulate phenotype value
      Y.hat <- matrix(NA,nObs,nT)
      Sigma <- matrix(NA,nObs,nT)
      for(n in 1:nObs){
        Sigma[n,] <- rnorm(nT,c(0,0),sigma.sqrt)
      }
      Y.hat <- Y.fit+Sigma
      rownames(Y.hat) <- paste("Indi_",sprintf("%03d",c(1:nObs)),sep="")
      Ftest_PB <- matrix(nrow=1,ncol=nM_V)	 #F scores of each SNP;
      Ftest_PB <- as.data.frame(Ftest_PB)
      colnames(Ftest_PB) <- names(CI_V)
      PrFt_PB <- Ftest_PB									 #Pr(>F) of each SNP;
      LOD_PB <- matrix(nrow=nT+1,ncol=nM_V)	 #LOD_M:Single Effect Statistic in Pleiotropic model;
      LOD_PB <- as.data.frame(LOD_PB)
      colnames(LOD_PB) <- names(CI_V)	
      rownames(LOD_PB) <- c(vecPheno,"Joint")
      for (q in 1:nM_V) {
        SNPq <- names(CI_V)[q]
        #make two models for every SNP
        #cofactors reducing
        if(nrow(CF) > 0){
          CCFvec4 <- tempCF_P	#CCFvec: Current CF_vector
          if (CChr %in% CF[,1]) {
            Cpos <- CI_V[q]
            CchrCF <- subset(CF,CF[,1] %in% CChr)
            CharC <- NULL
            for(ccc in 1:nrow(CchrCF)){
              if((Cpos >= CchrCF[ccc,3]) && (Cpos <= CchrCF[ccc,4])){
                CharC <- rownames(CchrCF)[ccc]
                CCFvec4 <- setdiff(CCFvec4,CharC)
              }
            }
          }
        }
        #Linear Model for different Situations:
        #1st. Multiple-line Cross and 2nd. Single-line Cross.
        if(length(CCFvec4) ==0 ){
          Xf2 <- GenoData_SI[,SNPq]
          Xf2 <- data.matrix(Xf2)
          if(nRP>1){
            lmr <- lm(Y.hat~Popu)
            lmf <- lm(Y.hat~Popu+Xf2:Popu)
          }else{
            lmr <- lm(Y.hat~1)
            lmf <- lm(Y.hat~1+Xf2)
          }
        }else{
          Xr2 <- TRMData[,CCFvec4]
          Xf2 <- cbind(Xr2,GenoData_SI[,SNPq])
          Xr2 <- data.matrix(Xr2)
          Xf2 <- data.matrix(Xf2)
          if(nRP>1){
            lmr <- lm(Y.hat~Popu+Xr2:Popu)
            lmf <- lm(Y.hat~Popu+Xf2:Popu)
          }else{
            lmr <- lm(Y.hat~1+data.matrix(Xr2))
            lmf <- lm(Y.hat~1+data.matrix(Xf2))
          }
        }			
        #get RSSp and RSSr and then pleiotropic model for Joint and single effect		
        resF_P <- resid(lmf)
        resR_P <- resid(lmr)
        RSSF_P <- t(resF_P) %*% resF_P
        RSSR_P <- t(resR_P) %*% resR_P
        detRSSF <- det(RSSF_P)
        detRSSR <- det(RSSR_P)
        RSSRs <- diag(RSSR_P)
        RSSFs <- diag(RSSF_P) 
        #residual df and Fvalue+Pvalue+JEs+SEs
        resdff <- df.residual(lmf)
        resdfr <- df.residual(lmr)
        rDFd <- resdfr - resdff						 #difference between rDF2 and rDF2
        #Fvalue <- (detRSSR-detRSSF)/(resdfr-resdff)/(detRSSF/resdff)
        #PrFt_P[q] <- 1-pf(Fvalue,resdfr-resdff,resdff)
        anvar <- anova(lmr,lmf)
        Ftest_PB[SNPq] <- anvar[2,"approx F"]
        PrFt_PB[SNPq] <- anvar[2,"Pr(>F)"]	
        LOD_PB[nT+1,SNPq] <- 0.5*nObs*log10(detRSSR/detRSSF)
        #Single Effect for each trait
        for(nt in 1:nT){
          LOD_PB[nt,SNPq] <- 0.5*nObs*log10(RSSRs[nt]/RSSFs[nt])
        }			
      }#q
      
      rese_S <- matrix(NA,nObs,nT)
      colnames(rese_S) <- vecPheno
      rownames(rese_S) <- rownames(Y.hat)
      
      for(nt in 1:nT){
        Xf_S <- data.frame(TRMData[CCFvec3],GenoData_SI[,names(which.max(LOD_PB[nt,]))])
        Xf_S <- data.matrix(Xf_S)
        if(nRP == 1) lms <- lm(Y.hat[,nt]~Xf_S) else lms <- lm(Y.hat[,nt]~Popu+Xf_S:Popu)
        rese_S[,nt] <- resid(lms)
      }
      Xf_P <- cbind(TRMData[CCFvec3],GenoData_SI[,names(which.max(LOD_PB[nT+1,]))])
      Xf_P <- data.matrix(Xf_P)
      if(nRP == 1) lme <- lm(Y.hat~Xf_P) else lme <- lm(Y.hat~Popu+Xf_P:Popu)
      RSSe_S <- t(rese_S)%*%rese_S
      RSSe_P <- t(resid(lme))%*%resid(lme)
      Vec_PB <- 0.5*nObs*log10(det(RSSe_P)/det(RSSe_S))
      PorL_PB[npb] <- Vec_PB
      if(npb %% (nPB/10)==0){
        cat(npb/nPB*100,"% ParametricBootstrap","have done\t\t @", as.character(Sys.time()),"\n")
      }
    }
    CriPL <- sort(PorL_PB)[nPB*(1-alpha)]
    if(PorL > CriPL) print("The testing QTLs are linked QTL") else print("The testing QTLs are pleiotropic QTL")
    MulTestTab[1,c("LOD_PL","lod_pl")] <- c(CriPL,PorL)
    
    return(MulTestTab) 
  }else{
    print("No significance for pleiotropic QTL model")
    return(MulTestTab)
  }
  #* --------------------------
  # create and set working dir
  #*---------------------------
  mainDir=getwd()
  subDir="OUTPUT_PL"
  ifelse(!file.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
  write.table(round(MulTestTab,3),file=paste(file.path(mainDir, subDir),"/MultiTrait_PL.txt",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t")	
  write.table(PorL_PB,file=paste(file.path(mainDir, subDir),"/MultiTrait_pbPL.txt",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t")	
}
