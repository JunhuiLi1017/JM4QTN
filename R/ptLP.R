ptLP <-
function(vecPheno,PhenoData,GenoData_EST,npt=1000,alpha=0.1,selection="stepwise",tolerance=1e-7,Trace="Pillai",select="SBC",sle=0.05,sls=0.05,Choose="SBC"){
  ##the data structure
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
  nT <- length(vecPheno)
  PhenoData[,"Popu"] <- as.numeric(PhenoData[,"Popu"])
  pTdata <- cbind(PhenoData[c(1:4)],PhenoData[c(vecPheno)])
  TRMdata <- cbind(pTdata,GenoData_EST)
  nRP <- nlevels(as.factor(TRMdata$Popu))  #nRP: No. of Real Populations; nRP > nlp
  if(nRP > 1){
    include0 <- "Popu"
    remove.cf <- c(1,2)
  }else{
    include0 <- NULL
    remove.cf <- 1
  }
  notX <- c(1,3,4)
  Class0 <- include0
  nObs <- nrow(TRMdata)
  rownames(TRMdata) <- TRMdata$Indi
  FullIndi <- rownames(TRMdata)
  ptTRMdata <- TRMdata
  Popu <- as.matrix(ptTRMdata[,"Popu"])
  Resid_L <- matrix(NA,nObs,nT)
  vecPrF <- rep(NA,npt)
  vecLOD <- matrix(NA,npt,2+nT)
  colnames(vecLOD) <- c(vecPheno,"LOD_P","LOD_L")
  nM_R <- ncol(GenoData_EST)
  MarkerID <- colnames(GenoData_EST)
  Ftest_P <- matrix(NA,nrow=nT+1,ncol=nM_R)
  colnames(Ftest_P) <- MarkerID
  rownames(Ftest_P) <- c(vecPheno,"Joint")
  LOD_M <- Ftest_P
  PrFt_M <- rep(NA,nM_R)
  names(PrFt_M) <- MarkerID
  dfptSeq <- c(1:nObs)
  for(v in 1:npt){
    ###1st,Every set of phenotypic values permutation within each population
    for(i in 1:nT){
      PerVec <- sample(FullIndi,nObs,replace=FALSE)
      ptTRMdata[1:nObs,vecPheno[i]] <- TRMdata[PerVec,vecPheno[i]]
    }
    tempCF_P <- stepwise(ptTRMdata, vecPheno, notX, include0, Class0, selection, select, sle, sls, tolerance, Trace, Choose)$variate
    CF_P <- tempCF_P[-remove.cf]
    ###3rd,catch cofactor and then QTL scanning
    Y <- data.matrix(ptTRMdata[,vecPheno])
    
    #Single trait model and pleiotropic model
    for(q in 1:nM_R){
      #make two models for every SNP
      SNPi <- MarkerID[q]
      tempVecSM <- CF_P[!CF_P %in% SNPi]
      
      if(nRP > 1){
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
      #Fvalue <- (detRSSR-detRSSF)/(resdfr-resdff)/(detRSSF/resdff)
      #PrFt_M[q] <- 1-pf(Fvalue,resdfr-resdff,resdff)
      PrFt_M[q] <- anova(lmf,lmr)[2,"Pr(>F)"]
      LOD_M[nT+1,q] <- 0.5*nObs*log10(detRSSR/detRSSF)
      #Single Effect for each trait
      for(i in 1:nT){
        LOD_M[i,q] <-0.5*nObs*log10(RSSR[i,i]/RSSF[i,i])
      }
    }#q
    vecPrF[v] <- min(PrFt_M)
    vecLOD[v,nT+1] <- max(LOD_M[nT+1,])
    ## Linked-QTL model
    for(i in 1:nT){
      vecLOD[v,i] <- max(LOD_M[i,])
      ## for Linkage model
      SNPi <- MarkerID[which.max(LOD_M[i,])]
      tempVecSM <- CF_P[!CF_P %in% SNPi]
      if(nRP > 1){
        if(length(tempVecSM)>0){
          lmr <- lm(Y~Popu+as.matrix(ptTRMdata[tempVecSM]):Popu)
          lmf <- lm(Y[,i]~Popu+as.matrix(ptTRMdata[c(tempVecSM,SNPi)]):Popu)  #Full model
        }else{
          lmr <- lm(Y~Popu)
          lmf <- lm(Y[,i]~Popu+as.matrix(ptTRMdata[c(SNPi)]):Popu)  #Full model
        } 
      }else{
        if(length(tempVecSM)>0){
          lmr <- lm(Y~1+as.matrix(ptTRMdata[tempVecSM]))
          lmf <- lm(Y[,i]~1+as.matrix(ptTRMdata[c(tempVecSM,SNPi)]))  #Full model
        }else{
          lmr <- lm(Y~1)
          lmf <- lm(Y[,i]~1+as.matrix(ptTRMdata[c(SNPi)]))  #Full model
        }
      }
      Resid_L[,i] <- resid(lmf)
    }
    ## reduced model
    Resid_R <- resid(lmr)
    RSSR_L <- t(Resid_R) %*% Resid_R
    detRSSR_L <- det(RSSR_L)
    ## full model for linkage
    RSSF_L <- t(Resid_L) %*% Resid_L
    detRSSF_L <- det(RSSF_L)
    
    dfR <- df.residual(lmr)
    dfF <- df.residual(lmf)
    #vecLOD[v,nT+2] <- (df.residual(lmr)-1 - 0.5*(nT-(dfR-dfF)))*log(detRSSR_L/detRSSF_L)
    vecLOD[v,nT+2] <- 0.5*nObs*abs(log10(detRSSR_L/detRSSF_L))
    if(v %% (npt/10)==0){
      cat(v/npt*100,"% permutation test have been completed\t\t","@",as.character(Sys.time()),"\n")
    }
  }#v 
  vecThrVal <- matrix(NA,3,nT+3)
  colnames(vecThrVal) <- c("Pr_Joint",paste("LOD_",vecPheno,sep=""),"LOD_Pleiotropic","LOD_Linked")
  rownames(vecThrVal) <- c("0.05","0.1",alpha)
  vecThrVal[,1] <- sort(vecPrF)[c(npt*0.05,npt*0.1,npt*alpha)]
  for(j in 1:(nT+2)){
    vecThrVal[,j+1] <- sort(vecLOD[,j])[c(npt*0.95,npt*0.9,npt*(1-alpha))]  
  }
  #* --------------------------
  # create and set working dir
  #*---------------------------
  mainDir=getwd()
  subDir="OUTPUT_ptLP"
  ifelse(!file.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
  write.table(cbind(vecPrF,vecLOD),file=paste(file.path(mainDir, subDir),"/PvalLOD_ptLP_",npt,sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE)
  write.table(vecThrVal,file=paste(file.path(mainDir, subDir),"/PvalLOD_ptLP",sep=""),quote=FALSE,col.names=TRUE,row.names=TRUE)
  return(vecThrVal)  
}
