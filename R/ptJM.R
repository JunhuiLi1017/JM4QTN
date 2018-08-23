ptJM <-
function(vecPheno,PhenoData,GenoData_EST,method,npt=1000,alpha=0.1,selection="stepwise",tolerance=1e-7,Trace="Pillai",select="SBC",sle=0.05,sls=0.05,Choose="SBC"){
  GenoData_EST <- GenoData_EST[-c(1:2),]
  if(nrow(PhenoData) != nrow(GenoData_EST)){
    stop("Sample size must be equal in PhenoData, GenoData_EST\n")
  }
  if(!all(vecPheno %in% colnames(PhenoData))){
    stop("Trait vector must be included in PhenoData\n")
  }
  if(all(colnames(PhenoData)[1:2] != c("Indi","Popu"))){
    stop("Top Four of column in PhenoData must be Indi\tPopu\tMA\tPA\n")
  }
  if( method!="AM" & method != "LM"){
    stop("Method must be AM or LM\n")
  }
  #* ---------------------------
  #* Permutation test with stepwise regression
  #* ---------------------------
  nT <- length(vecPheno)
  pTdata <- cbind(PhenoData[c(1:4)],PhenoData[c(vecPheno)])
  TRMdata <- cbind(pTdata,GenoData_EST)
  ### catch population levels information
  nRP <- nlevels(as.factor(TRMdata$Popu))  #nRP: No. of Real Populations; nRP > nlp
  if(nRP > 1) nest <- TRUE else nest <- FALSE
  vecPop <- levels(as.factor(TRMdata$Popu))  #vecPop: vector of Populations
  nObs <- nrow(TRMdata)
  rownames(TRMdata) <- TRMdata$Indi
  FullIndi <- rownames(TRMdata)
  TRMdata[,"Popu"] <- as.numeric(TRMdata[,"Popu"])
  ptTRMdata <- TRMdata
  #cofactor selection with R package StepReg
  Popu <- as.matrix(ptTRMdata["Popu"])
  vecPrF <- matrix(,npt,nT)
  vecPrF <- as.data.frame(vecPrF)
  colnames(vecPrF) <- paste(vecPheno,"_PrF",sep="")
  vecLOD <- vecPrF
  colnames(vecLOD) <- paste(vecPheno,"_LOD",sep="")
  nM_R <- ncol(GenoData_EST)
  MarkerID <- colnames(GenoData_EST)
  Ftest_P <- matrix(nrow=nT,ncol=nM_R)
  colnames(Ftest_P) <- MarkerID
  rownames(Ftest_P) <- vecPheno
  LOD_S <- Ftest_P
  PrFt_S <- Ftest_P
  dfptSeq <- c(1:nObs)
  for(v in 1:npt){
    ###1st,Every set of phenotypic values permutation within each population
    PerVec <- sample(FullIndi,nObs,replace=FALSE)
    ptTRMdata[1:nObs,vecPheno] <- TRMdata[PerVec,vecPheno]
    ###2nd,cofactors selection by Calling RCpp
    for(j in 1:nT){
      #j=1
      if(nT>1){
        ptTRMdata_t1 <- ptTRMdata[,-(which(!vecPheno %in% vecPheno[j])+4)]
      }else{
        ptTRMdata_t1 <- ptTRMdata
      }
      y <- vecPheno[j]
      Y <- as.matrix(ptTRMdata[y])
      notX <- c(1,3,4)
      include0 <- "Popu"
      Class0 <- include0
      tempCF_P <- stepwise(ptTRMdata_t1, y, notX, include0, Class0, selection, select, sle, sls, tolerance, Trace, Choose)$variate
      CF_P <- tempCF_P[-c(1,2)]
      ###3rd,catch cofactor and then QTL scanning
      for(q in 1:nM_R){
        #make two models for every SNP
        SNPi <- MarkerID[q]
        tempVecSM <- CF_P[!CF_P %in% SNPi]
        if(nest==TRUE){
          if(length(tempVecSM)>0){
            lmr <- lm(Y~Popu+as.matrix(ptTRMdata[tempVecSM]):Popu)
            lmf <- lm(Y~Popu+as.matrix(ptTRMdata[c(tempVecSM,SNPi)]):Popu)  #Full model
          }else{
            lmr <- lm(Y~Popu)
            lmf <- lm(Y~Popu+as.matrix(ptTRMdata[c(SNPi)]):Popu)  #Full model
          }
        }else{
          if(length(tempVecSM)>0){
            lmr <- lm(Y~1+as.matrix(ptTRMdata[tempVecSM]))
            lmf <- lm(Y~1+as.matrix(ptTRMdata[c(tempVecSM,SNPi)]))  #Full model
          }else{
            lmr <- lm(Y~1)
            lmf <- lm(Y~1+as.matrix(ptTRMdata[c(SNPi)]))  #Full model
          }
        }
        #run two models with single trait separatey
        resF <- resid(lmf)
        resR <- resid(lmr)
        #get RSSp and RSSr and then pleiotropic model for Joint and single effect    
        RSSF <- t(resF) %*% resF
        RSSR <- t(resR) %*% resR
        detRSSF <- det(RSSF)
        detRSSR <- det(RSSR) 
        #residual df and Fvalue+Pvalue+JEs+SEs
        resdff <- df.residual(lmf)
        resdfr <- df.residual(lmr)
        PrFt_S[j,q] <- anova(lmf,lmr)[2,"Pr(>F)"]
        LOD_S[j,q] <- 0.5*nObs*log10(detRSSR/detRSSF)
        #LOD_S[j,q] <- (resdff-0.5)*log10(detRSSR/detRSSF)
      }#q
      vecPrF[v,j] <- min(PrFt_S[j,])
      vecLOD[v,j] <- max(LOD_S[j,])
    }#nT
    if(v %% (npt/10)==0){
      cat(v/npt*100,"% have been completed\t\t","@",as.character(Sys.time()),"\n")
    }
  }#v
  if(method=="LM"){
    vecThrVal <- matrix(NA,3,2*nT)
    colnames(vecThrVal) <- c(paste("Pr_",vecPheno,sep=""),paste("LOD_",vecPheno,sep=""))
    rownames(vecThrVal) <- c("0.05","0.1",alpha)
    for(j in 1:nT){
      vecThrVal[,j] <- sort(vecPrF[,j])[c(npt*0.05,npt*0.1,npt*alpha)]
      vecThrVal[,j+nT] <- sort(vecLOD[,j])[c(npt*0.95,npt*0.9,npt*(1-alpha))]  
    }
  }else if(method=="AM"){
    vecThrVal <- matrix(NA,3,nT)
    colnames(vecThrVal) <- paste("Pr_",vecPheno,sep="")
    rownames(vecThrVal) <- c("0.05","0.1",alpha)
    for(j in 1:nT){
      vecThrVal[,j] <- sort(vecPrF[,j])[c(npt*0.05,npt*0.1,npt*alpha)]
    }
  }
  #* --------------------------
  # create and set working dir
  #*---------------------------
  mainDir=getwd()
  subDir="OUTPUT_ptJM"
  ifelse(!file.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
  # create trait JM result in out_JM directory
  vecPrFLOD <- cbind(vecPrF,vecLOD)
  if(method=="AM"){
    write.table(vecPrF,file=paste(file.path(mainDir, subDir),"/Pval_ptJM_",npt,".xls",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t")
    write.table(vecThrVal,file=paste(file.path(mainDir, subDir),"/Pval_ptJM",".xls",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t") 
  }else if (method=="LM"){
    write.table(vecPrFLOD,file=paste(file.path(mainDir, subDir),"/PvalLOD_ptJM_",npt,".xls",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t")
    write.table(vecThrVal,file=paste(file.path(mainDir, subDir),"/PvalLOD_ptJM",".xls",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE,sep="\t")
  }
  return(vecThrVal)
}
