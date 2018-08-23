JM <- function(vecPheno,vecH2,PhenoData,method,vecThrVal,GenoData_EST,GenoData_QTL=NULL){
  if(length(vecPheno) != length(vecThrVal)/2 && length(vecPheno) != length(vecH2)){
    stop("No. of trait in vecPheno must be equal to length of vector vecThrVal or vector of heredibility\n")
  }
  if(nrow(PhenoData) != nrow(GenoData_EST)-2){
    stop("Sample size must be equal in PhenoData and GenoData_EST\n")
  }
  if(!all(vecPheno %in% colnames(PhenoData))){
    stop("Trait vector must be included in PhenoData\n")
  }
  if(all(colnames(PhenoData)[1:2] != c("Indi","Popu"))){
    stop("Top Four of column in PhenoData must be Indi\tPopu\tMA\tPA\n")
  }
  if(method == "AM"){
    GenoData_QTL <- GenoData_EST
  }else if(method == "LM"){
    if(is.null(GenoData_QTL) || ncol(GenoData_QTL) < ncol(GenoData_EST)){
      stop("GenoData_QTL should be prepared by function calGenoProb")
    }
  }else{
    stop("method must be AM or LM")
  }
  nT <- length(vecPheno)
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
  NoChr <- nlevels(as.factor(mapORpos$chr))
  nM_V <- ncol(GenoData_QTL)
  MarkerID <- colnames(GenoData_QTL)
  nObs <- nrow(PhenoData)
  nRP <- nlevels(as.factor(PhenoData$Popu))
  vecPop <- levels(as.factor(PhenoData$Popu))
  Delta <- 1e-04
  Summary_D <- matrix(nrow=nT,ncol=8)
  colnames(Summary_D) <- c("nCF","naQ","R2_aQ","R2j_aQ","pG_aR2.h2","pG_aR2j.h2","neQ","R2_eQ")
  rownames(Summary_D) <- vecPheno
  vecThrVal <- as.numeric(vecThrVal)
  CriPrF <- vecThrVal[1:nT]
  Ftest_P <- matrix(nrow=nT,ncol=nM_V)
  Ftest_P <- as.data.frame(Ftest_P)
  colnames(Ftest_P) <- MarkerID
  rownames(Ftest_P) <- vecPheno
  PrFt_P <- Ftest_P
  detRSSF_P <- Ftest_P
  LOD_P <- Ftest_P
  CF_P <- NULL
  for(j in 1:nT){
    TRMdata <- cbind(PhenoData[c(1:4)],PhenoData[vecPheno[j]],GenoData_EST)
    TRMdata[,"Popu"] <- as.numeric(TRMdata[,"Popu"])
    selection='stepwise'
    tolerance = 1e-07
    Trace = "Pillai"
    select="SL"
    sle=CriPrF[j]
    sls=CriPrF[j]
    Choose="SBC"
    if(nRP > 1){
      include0 <- "Popu"
      remove.cf <- c(1,2)
    }else{
      include0 <- NULL
      remove.cf <- 1
    }
    Popu <- data.matrix(TRMdata[,"Popu"])
    y <- vecPheno[j]
    notX <- c(1,3,4)
    Class0 <- include0
    tempCF_P <- stepwise(TRMdata, y, notX, include0, Class0, selection, select, sle, sls, tolerance, Trace, Choose)$variate    
    tempCF_P <- tempCF_P[-remove.cf]
    if(length(tempCF_P)==0){
      warning("There is no cofactor for trait\t",vecPheno[j],"\twith stepwise regression for significant level\t",vecThrVal[j],"\n")
    }
    CF <- matrix(nrow=0,ncol=4)
    CF <- as.data.frame(CF)
    if(length(tempCF_P) > 0){
      for(k in 1:length(tempCF_P)){
        for(i in 1:NoChr){
          ChrMarkerID <- as.character(subset(mapORpos,mapORpos$chr==i)[,1])
          ChrMarkerPos <- subset(mapORpos,mapORpos$chr==i)[,3]
          if(tempCF_P[k] %in% ChrMarkerID){
            Cchr <- i
            Mseq <- which(ChrMarkerID %in% tempCF_P[k])
            pos <- ChrMarkerPos[Mseq]
            if(Mseq == 1){
              Lpos <- 0
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
        }
      }
      rownames(CF) <- tempCF_P
      colnames(CF) <- c("chr","pos","Lpos","Rpos")
    }
    Summary_D[j,1] <- nrow(CF)
    CF_P[j] <- list(CF)
    Y <- data.matrix(PhenoData[,y])
    TRVMData <- cbind(PhenoData[c(1:4)],PhenoData[vecPheno[j]],GenoData_QTL)
    for(i in 1:NoChr) {
      qtl <- subset(vmapORpos,vmapORpos$chr==i)
      qEnd <- nrow(qtl)
      for(q in 1:qEnd) {
        SNPq <- qtl[q,1]
        #make two models for every SNP
        #cofactors reducing
        if (nrow(CF) > 0) {
          CCFvec <- rownames(CF)
          if (i %in% CF[, 1]) {
            Cpos <- qtl[q,3]
            CchrCF <- subset(CF, CF[, 1] %in% i)
            CharC <- NULL
            for (ccc in 1:nrow(CchrCF)) {
              if ((Cpos >= CchrCF[ccc, 3]) && (Cpos <= CchrCF[ccc, 4])) {
                CharC <- rownames(CchrCF)[ccc]
                CCFvec <- CCFvec[which(!CCFvec %in% CharC)]
              }
            }
          }
          if (length(CCFvec) == 0) {
            Xf <- data.matrix(TRVMData[SNPq])
          } else {
            Xr <- data.matrix(TRVMData[c(CCFvec)])
            Xf <- data.matrix(TRVMData[c(CCFvec, SNPq)])
          }
        }else{
          Xf <- data.matrix(TRVMData[SNPq])
        }
        
        if(nRP>1){
          if ((nrow(CF) > 0 && length(CCFvec) == 0) || nrow(CF) == 0) {
            lmr <- lm(Y ~ Popu)
            lmf <- lm(Y ~ Popu + Xf:Popu)
          }
          if(nrow(CF) > 0 && length(CCFvec) != 0){
            lmr <- lm(Y ~ Popu + Xr:Popu)
            lmf <- lm(Y ~ Popu + Xf:Popu)
          }
        }else {
          if ((nrow(CF) > 0 && length(CCFvec) == 0) || nrow(CF) == 0) {
            lmr <- lm(Y ~ 1)
            lmf <- lm(Y ~ 1 + Xf)
          }
          if(nrow(CF) > 0 && length(CCFvec) != 0){
            lmr <- lm(Y ~ 1 + Xr)
            lmf <- lm(Y ~ 1 + Xf)
          }
        }
        rDFf <- df.residual(lmf)
        rDFr <- df.residual(lmr)
        rDFd <- rDFr - rDFf
        anvar <- anova(lmr,lmf)
        Ftest_P[j,SNPq] <- anvar[2,"F"]
        PrFt_P[j,SNPq] <- anvar[2,"Pr(>F)"]
        LOD_P[j,SNPq] <- 0.5*nObs*log10(1+(rDFd/rDFf)*Ftest_P[j,SNPq])
      }
    }
  }
  if(method=="LM"){
    CriLOD_D <- vecThrVal[nT+1:nT]   #PermuTest, alpha=0.10 level
    CID <- 1.0   #Conf.Inter.Degrade
    CriQRD <- 5.0   #Criterion of QTL redundant distance(in cM)
  }else if(method=="AM"){
    PrBonfferoni <- vecThrVal[nT+1:nT]/ncol(GenoData_QTL) ##Bonfferoni correct
  }
  if (method == "LM") {
    CriLOD_D <- vecThrVal[nT + 1:nT]
    CID <- 1
    CriQRD <- 5
  }else if (method == "AM") {
    PrBonfferoni <- vecThrVal[nT + 1:nT]/ncol(GenoData_QTL)
  }
  ltaQTN_D <- NULL
  
  for (j in 1:nT) {
    TRVMData <- cbind(PhenoData[c(1:4)],PhenoData[vecPheno[j]],GenoData_QTL)
    tempQTN <- NULL
    QnameVec <- NULL
    EstQTN <- NULL
    if(method=="LM"){
      for (i in 1:NoChr){
        qtl <- subset(vmapORpos,vmapORpos$chr==i)
        CCend <- nrow(qtl)
        CCLOD <- LOD_P[j,which(colnames(LOD_P) %in% qtl[,1])]
        PreTestVec <- CCLOD[which(CCLOD >= CriLOD_D[j])]
        if(length(PreTestVec) > 0){
          for(k in 1:length(PreTestVec)){
            CPeak <- PreTestVec[k]
            CNr <- which( qtl[,1] %in% names(PreTestVec)[k])
            Cpos <- as.numeric(qtl[CNr,3])
            Lflag <- 0
            Rflag <- 0
            if(CNr > 1){
              if( (is.na(CCLOD[CNr-1])) | (!(is.na(CCLOD[CNr-1])) && (CCLOD[CNr-1] < CPeak)) ){
                Lflag <- 1
              }
            }else{
              Lflag <- 1
            }
            if(CNr < CCend){
              if( (is.na(CCLOD[CNr+1])) | (!(is.na(CCLOD[CNr+1])) && (CCLOD[CNr+1] < CPeak)) ){
                Rflag <- 1
              }
            }else{
              Rflag <- 1
            }
            if((Lflag == 1) && (Rflag == 1)){
              LNr <- CNr
              if(CNr > 1){
                while((LNr-1 > 0) && ((is.na(CCLOD[LNr-1]))|(!(is.na(CCLOD[LNr-1]))&&(CCLOD[LNr-1] > CPeak-CID)))){
                  LNr <- LNr-1
                }
              }
              if(LNr == 1){
                LBpos <- round(as.numeric(qtl[LNr,3]),digits=3)
              }else{
                LBpos <- round(mean(as.numeric(qtl[c(LNr-1,LNr),3])),digits=3)
              }
              RNr <- CNr
              if(CNr < CCend){
                while((RNr+1 < CCend) && ((is.na(CCLOD[RNr+1]))|(!(is.na(CCLOD[RNr+1]))&&(CCLOD[RNr+1] > CPeak-CID)))){
                  RNr <- RNr+1
                }
              }
              if(RNr == CCend){
                RBpos <- round(as.numeric(qtl[RNr,3]),digits=3)
              }else{
                RBpos <- round(mean(as.numeric(qtl[c(RNr,RNr+1),3])),digits=3)
              }
              if(length(QnameVec)==0){
                tempQTN <- c(CPeak[1,1],i,Cpos,LBpos,RBpos)
                tempQTN <- as.data.frame(t(tempQTN))
                QnameVec <- append(QnameVec,names(PreTestVec)[k])
              }else{
                NrLR <- nrow(tempQTN)
                if((tempQTN[NrLR,2] == i) && ((Cpos - tempQTN[NrLR,3]) <= CriQRD)){
                  if(tempQTN[NrLR,1] <= CPeak[1,1]){
                    tempQTN[NrLR,] <- c(CPeak[1,1],i,Cpos,LBpos,RBpos)
                    tempQTN <- as.data.frame(tempQTN)
                    QnameVec[NrLR] <- names(PreTestVec)[k]
                  }
                }else{
                  tempQTN <- rbind(tempQTN,c(CPeak[1,1],i,Cpos,LBpos,RBpos))
                  tempQTN <- as.data.frame(tempQTN)
                  QnameVec <- append(QnameVec,names(PreTestVec)[k])
                }
              }
            }
          }
        }
      }
      if(length(QnameVec) >= 1){
        rownames(tempQTN) <- QnameVec
        colnames(tempQTN) <- c("LOD","chr","pos","LBpos","RBpos")
      }
    }else if(method=="AM"){
      QnameVec <- mapORpos[PrFt_P[j,] < PrBonfferoni[j],1]
      if(length(QnameVec) != 0){
        tempQTN <- cbind(t(PrFt_P[j,mapORpos[,1] %in% QnameVec]),mapORpos[mapORpos[,1] %in% QnameVec,2:3])
        if(length(QnameVec) >= 1){
          rownames(tempQTN) <- QnameVec
          colnames(tempQTN) <- c("Pr","chr","pos")
        }
      }
    }
    if(length(QnameVec) == 0){
      warning("There is no QTN detected for trait\t",vecPheno[j],"\twith significant level\t",vecThrVal[j],"\tand lod\t",vecThrVal[nT+j],"\n")
    }
    if (length(QnameVec)!=0) {
      tempQTN_Ori <- rownames(tempQTN)
      if (length(QnameVec) > 0) {
        tempQTN <- tempQTN[order(tempQTN[,1],decreasing = TRUE),]
      }
    }
    if(length(QnameVec) == 0){
      Summary_D[j,2] <- 0    #naQ_L
      Summary_D[j,3] <- 0    #R2_aQ_L
      Summary_D[j,4] <- 0    #R2j_aQ_L
      Summary_D[j,5] <- 0    #p_aR2.h2_L
      Summary_D[j,6] <- 0    #p_aR2j.h2_L
    }else{
      if(nRP>1){
        a_QF <- paste(vecPheno[j]," ~ Popu + ",rownames(tempQTN)[1],":Popu",sep="")
        if (nrow(tempQTN) > 1) {
          for (a in 2:nrow(tempQTN)){
            a_QF <- paste(a_QF," + ",rownames(tempQTN)[a],":Popu",sep="")
          }
        }
      }else{
        a_QF <- paste(vecPheno[j]," ~ 1 + ",rownames(tempQTN)[1],sep="")
        if (nrow(tempQTN) > 1) {
          for (a in 2:nrow(tempQTN)){
            a_QF <- paste(a_QF," + ",rownames(tempQTN)[a],sep="")
          }
        }
      }
      lma_QF <- lm(as.formula(a_QF),data=TRVMData)
      vecSS <- anova(lma_QF)   #pOP + qTLs + rESIDUAL
      vecSS1 <- vecSS[,2]   #pOP + qTLs + rESIDUAL
      SSt <- sum(vecSS1)             #SS_total
      if(nRP > 1){
        SSp <- vecSS1[1]               #SS_Pop
        vecSSqr <- vecSS1[-1] 
        SSqr <- sum(vecSSqr)          #SS_(QTLs+Residual)
        SSr <- vecSSqr[length(vecSSqr)] #SS_Res.
        vecSSq <- vecSSqr[-length(vecSSqr)]
        SSq <- sum(vecSSq)            #SS_QTLs
        R2_aQ <- SSq/SSqr
      }else{
        vecSSqr <- vecSS1
        SSqr <- sum(vecSSqr)          #SS_(QTLs+Residual)
        SSr <- vecSSqr[length(vecSSqr)] #SS_Res.
        vecSSq <- vecSSqr[-length(vecSSqr)]
        SSq <- sum(vecSSq)            #SS_QTLs
        R2_aQ <- SSq/SSqr
      }
      nPar <- nrow(tempQTN)
      adjR2_aQ <- 1-(1 - R2_aQ)*(nObs-1)/(nObs-nPar-1)
      Summary_D[j,2] <- nrow(tempQTN)     #naQ_L
      Summary_D[j,3] <- round(R2_aQ, digits = 3)     #R2_aQ_L
      Summary_D[j,4] <- round(adjR2_aQ, digits = 3) #R2j_aQ_L
      Summary_D[j,5] <- round(100*R2_aQ/vecH2[j], digits = 1)    #p_aR2.h2_L
      Summary_D[j,6] <- round(100*adjR2_aQ/vecH2[j], digits = 1)    #p_aR2j.h2_L
      #get every QTL's coef. and its P-value of significance
      mtQTLeff <- summary(lma_QF)$coefficients[-1,c(1,4)]
      if(class(mtQTLeff)=="numeric"){
        mtQTLeff <- t(as.matrix(mtQTLeff))
      }
      colnames(mtQTLeff) <- c("eff","prob")
      dfCoeP <- matrix(0,nrow=nrow(tempQTN),ncol=nRP)
      dfCoeP <- as.data.frame(dfCoeP)
      rownames(dfCoeP) <- rownames(tempQTN)
      colnames(dfCoeP) <- paste("Coe_Pop",c(1:nRP),sep="")
      vecqNM <- rownames(tempQTN)
      if(nRP > 1){
        for(k in 1:nrow(tempQTN)) {
          qNM <- vecqNM[k]
          for(p in 1:nRP){
            CENM <- paste("Popu",vecPop[p],":",qNM,sep="")
            if((CENM %in% rownames(mtQTLeff)) && (mtQTLeff[CENM,2] < 0.15)){
              dfCoeP[qNM,p] <- mtQTLeff[CENM,1]
            }#if
          }#p
        }#k
      }else{
        for(k in 1:nrow(tempQTN)) {
          if(mtQTLeff[k,2] < 0.15){
            dfCoeP[k,1] <- mtQTLeff[k,1]
          }#if
        }#k
      }
      SSqtl <- vecSS[-c(1,nrow(vecSS)),]
      if(nRP>1){
        noNASSqtl <- unlist(strsplit(rownames(SSqtl),":"))
        SSqtlName <- noNASSqtl[seq(2,length(noNASSqtl),2)]
        rownames(SSqtl) <- SSqtlName
      }else{
        SSqtlName <- rownames(vecSS)[-c(1,nrow(vecSS))]
      }
      tempQTNname <- rownames(tempQTN)
      tempPR2 <- matrix(0,length(tempQTNname),1)
      rownames(tempPR2) <- tempQTNname
      
      for(nk in 1:length(tempQTNname)){
        if(tempQTNname[nk] %in% SSqtlName){
          p1 <- which(SSqtlName %in% tempQTNname[nk])
          tempPR2[nk] <- SSqtl[p1,2]
        }
      }
      tempPR2 <- tempPR2/SSqr
      EstQTN <- data.frame(tempQTN,dfCoeP,"PartialR2"=tempPR2,"pG"=100*tempPR2/vecH2[j])
      EstQTN <- EstQTN[tempQTN_Ori,]  #if sorted, please release this line
    }#nrow(tempQTN) = 0
    ltaQTN_D[j] <- list(EstQTN)
  }#j
  names(ltaQTN_D) <- vecPheno
  mainDir=getwd()
  subDir="OUTPUT_JM"
  ifelse(!file.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
  if(method=="AM"){
    write.table(t(PrFt_P),file=paste(file.path(mainDir, subDir),"/",paste(vecPheno,collapse=''),"_Pr.xls",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t")
  }else if(method=="LM"){
    write.table(t(LOD_P),file=paste(file.path(mainDir, subDir),"/",paste(vecPheno,collapse=''),"_LOD.xls",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t")
  }
  for(j in 1:nT){
    write.table(ltaQTN_D[j],file=paste(file.path(mainDir, subDir),"/",vecPheno[j],"_Significant.xls",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t")
  }
  return(ltaQTN_D)
}