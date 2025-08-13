#' Permutation Test for Joint Mapping Analysis
#' 
#' Performs permutation tests for joint mapping analysis to determine empirical significance thresholds
#' for QTL detection using stepwise regression. This function implements a comprehensive permutation
#' testing framework for both Association Mapping (AM) and Linkage Mapping (LM) methods.
#' 
#' @param vecPheno Character vector specifying the trait names to analyze.
#' @param PhenoData Data frame containing phenotype data with columns:
#'   \itemize{
#'     \item \code{Indi}: Individual identifiers
#'     \item \code{Popu}: Population identifiers
#'     \item \code{MA}: Maternal allele information
#'     \item \code{PA}: Paternal allele information
#'     \item Additional columns: Trait measurements
#'   }
#' @param GenoData_EST Data frame containing genotype data for cofactor selection.
#' @param method Character string specifying the mapping method:
#'   \itemize{
#'     \item \code{"AM"}: Association Mapping - uses P-values for significance testing
#'     \item \code{"LM"}: Linkage Mapping - uses LOD scores for significance testing
#'   }
#' @param npt Numeric value specifying the number of permutations. Default is 1000.
#'   Higher values provide more accurate threshold estimates but require more computation time.
#' @param alpha Numeric value specifying the significance level for threshold calculation.
#'   Must be between 0 and 1. Default is 0.1.
#' @param selection Character string specifying the selection method for stepwise regression.
#'   Default is "stepwise".
#' @param tolerance Numeric value specifying the tolerance for convergence in stepwise regression.
#'   Default is 1e-7.
#' @param Trace Character string specifying the trace method for multivariate analysis.
#'   Default is "Pillai".
#' @param select Character string specifying the selection criterion for stepwise regression.
#'   Default is "SBC" (Schwarz Bayesian Criterion).
#' @param sle Numeric value specifying the significance level for entry in stepwise regression.
#'   Default is 0.05.
#' @param sls Numeric value specifying the significance level for stay in stepwise regression.
#'   Default is 0.05.
#' @param Choose Character string specifying the model selection criterion.
#'   Default is "SBC" (Schwarz Bayesian Criterion).
#' 
#' @return A matrix containing permutation test thresholds with:
#'   \itemize{
#'     \item \strong{Rows}: Significance levels (0.05, 0.1, alpha)
#'     \item \strong{Columns}: Traits and statistics
#'     \item \strong{For AM method}: P-value thresholds only
#'     \item \strong{For LM method}: P-value and LOD score thresholds
#'   }
#' 
#' @details This function implements a comprehensive permutation testing framework:
#' 
#' \strong{Permutation Testing Process:}
#' \enumerate{
#'   \item \strong{Data Validation}: Checks input data consistency and structure
#'   \item \strong{Phenotypic Permutation}: Randomly shuffles phenotypic values within populations
#'   \item \strong{Cofactor Selection}: Applies stepwise regression for each permutation
#'   \item \strong{QTL Scanning}: Calculates F-statistics and LOD scores for each permutation
#'   \item \strong{Threshold Determination}: Establishes empirical significance thresholds
#' }
#' 
#' \strong{Population Handling:}
#' \itemize{
#'   \item \strong{Single Population}: Simple permutation within the population
#'   \item \strong{Multiple Populations}: Permutation within each population separately
#'   \item \strong{Population Effects}: Accounts for population structure in the analysis
#' }
#' 
#' \strong{Statistical Methods:}
#' \itemize{
#'   \item \strong{Stepwise Regression}: Uses the StepReg package for cofactor selection
#'   \item \strong{Multivariate Analysis}: Handles multiple traits simultaneously
#'   \item \strong{Model Comparison}: Compares reduced and full models for significance testing
#'   \item \strong{Multiple Testing Correction}: Accounts for multiple comparisons
#' }
#' 
#' \strong{Output Interpretation:}
#' \itemize{
#'   \item \strong{P-values}: Empirical significance thresholds for association mapping
#'   \item \strong{LOD scores}: Empirical significance thresholds for linkage mapping
#'   \item \strong{Multiple levels}: Provides thresholds at 0.05, 0.1, and user-specified alpha levels
#' }
#' 
#' @examples
#' # Example phenotype data
#' pheno_data <- data.frame(
#'   Indi = paste0("Ind", 1:100),
#'   Popu = rep(c("Pop1", "Pop2"), each = 50),
#'   MA = rep(1:2, 50),
#'   PA = rep(1:2, 50),
#'   Trait1 = rnorm(100, mean = 100, sd = 15),
#'   Trait2 = rnorm(100, mean = 50, sd = 8)
#' )
#' 
#' # Example genotype data
#' geno_data <- matrix(sample(c(0,1,2), 100*50, replace = TRUE), 
#'                     nrow = 100, ncol = 50)
#' colnames(geno_data) <- paste0("M", 1:50)
#' 
#' # Permutation test for association mapping
#' thresholds_am <- permutation_test_joint_mapping(c("Trait1", "Trait2"), pheno_data, geno_data, 
#'                      method = "AM", npt = 100, alpha = 0.1)
#' 
#' # Permutation test for linkage mapping
#' thresholds_lm <- permutation_test_joint_mapping(c("Trait1", "Trait2"), pheno_data, geno_data, 
#'                      method = "LM", npt = 100, alpha = 0.1)
#' 
#' # View results
#' print(thresholds_am)
#' print(thresholds_lm)
#' 
#' # Example with different parameters
#' thresholds_custom <- permutation_test_joint_mapping(c("Trait1", "Trait2"), pheno_data, geno_data,
#'                          method = "AM", npt = 500, alpha = 0.05,
#'                          sle = 0.01, sls = 0.01)
#' 
#' @seealso \code{\link{joint_mapping_analysis}} for joint mapping analysis,
#'          \code{\link{permutation_test_linkage_pleiotropy}} for linkage-pleiotropy permutation tests
#' 
#' @references
#' Churchill, G.A. and Doerge, R.W. (1994). Empirical threshold values for quantitative 
#' trait mapping. Genetics, 138(3), 963-971.
#' 
#' @export
permutation_test_joint_mapping <-
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
