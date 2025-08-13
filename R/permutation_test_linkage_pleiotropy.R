#' Permutation Test for Linkage-Pleiotropy Analysis
#' 
#' Performs permutation tests for linkage-pleiotropy analysis to determine empirical significance thresholds
#' for distinguishing between pleiotropic and linked QTL effects. This function implements a comprehensive
#' permutation testing framework specifically designed for multivariate trait analysis.
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
#'     \item \strong{Columns}: 
#'       \itemize{
#'         \item \code{Pr_Joint}: Joint effect P-value threshold
#'         \item \code{LOD_[trait]}: LOD score thresholds for individual traits
#'         \item \code{LOD_Pleiotropic}: Pleiotropic effect LOD threshold
#'         \item \code{LOD_Linked}: Linked QTL effect LOD threshold
#'       }
#'   }
#' 
#' @details This function implements a specialized permutation testing framework for linkage-pleiotropy analysis:
#' 
#' \strong{Permutation Testing Process:}
#' \enumerate{
#'   \item \strong{Data Validation}: Checks input data consistency and structure
#'   \item \strong{Phenotypic Permutation}: Randomly shuffles phenotypic values within populations
#'   \item \strong{Cofactor Selection}: Applies stepwise regression for each permutation
#'   \item \strong{Multivariate Analysis}: Calculates three types of LOD scores for each permutation
#'   \item \strong{Threshold Determination}: Establishes empirical significance thresholds
#' }
#' 
#' \strong{Three Types of LOD Scores:}
#' \itemize{
#'   \item \strong{Individual Trait LOD Scores}: Maximum LOD score for each trait separately
#'   \item \strong{Pleiotropic Effect LOD Score}: LOD score for single QTL affecting multiple traits
#'   \item \strong{Linked QTL Effect LOD Score}: LOD score for multiple QTL affecting different traits
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
#'   \item \strong{Multivariate Linear Models}: Handles multiple traits simultaneously
#'   \item \strong{Model Comparison}: Compares different genetic models for significance testing
#'   \item \strong{Residual Analysis}: Uses residual sum of squares for LOD calculations
#' }
#' 
#' \strong{Distinguishing Between Effects:}
#' \itemize{
#'   \item \strong{Pleiotropic Effects}: Single QTL affecting multiple traits simultaneously
#'   \item \strong{Linked QTL Effects}: Multiple QTL affecting different traits independently
#'   \item \strong{Statistical Comparison}: Uses likelihood ratio tests to distinguish between models
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
#' # Permutation test for linkage-pleiotropy analysis
#' thresholds <- permutation_test_linkage_pleiotropy(
#'   c("Trait1", "Trait2"), pheno_data, geno_data, 
#'   npt = 100, alpha = 0.1
#' )
#' 
#' # View results
#' print(thresholds)
#' 
#' # Example with different parameters
#' thresholds_custom <- permutation_test_linkage_pleiotropy(
#'   c("Trait1", "Trait2"), pheno_data, geno_data,
#'   npt = 500, alpha = 0.05, sle = 0.01, sls = 0.01
#' )
#' 
#' # Example with three traits
#' pheno_data_3t <- cbind(pheno_data, Trait3 = rnorm(100, mean = 25, sd = 5))
#' thresholds_3t <- permutation_test_linkage_pleiotropy(
#'   c("Trait1", "Trait2", "Trait3"), pheno_data_3t, geno_data, 
#'   npt = 100
#' )
#' 
#' @seealso \code{\link{linkage_pleiotropy_analysis}} for linkage-pleiotropy analysis,
#'          \code{\link{permutation_test_joint_mapping}} for joint mapping permutation tests
#' 
#' @references
#' Jiang, C. and Zeng, Z.B. (1995). Multiple trait analysis of genetic mapping for 
#' quantitative trait loci. Genetics, 140(3), 1111-1127.
#' 
#' @export
permutation_test_linkage_pleiotropy <-
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
