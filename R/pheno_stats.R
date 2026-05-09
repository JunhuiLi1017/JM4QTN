#' Statistical Analysis of Phenotype Data
#' 
#' Performs comprehensive statistical analysis on phenotype data including normality tests, 
#' analysis of variance (ANOVA), and least squares means calculations for genetic studies.
#' 
#' @param phenoData A data frame containing phenotype data with at least 5 columns.
#'   The first 4 columns must be: Environment (E), Block (B), Repetition (R), and Genotype (G).
#'   Additional columns contain trait measurements for statistical analysis.
#' @param defineForm Optional character vector of custom formula strings for statistical analysis.
#'   If NULL, default formulas are automatically generated based on the data structure.
#' @param effNotation Character string specifying the effect notation for least squares means.
#'   Default is "G" for genotype effect.
#' 
#' @return A list containing comprehensive statistical analysis results for each trait:
#'   \item{normality_test}{Results of Shapiro-Wilk normality test}
#'   \item{formula}{Formula used for the statistical model}
#'   \item{ANOVA}{Complete analysis of variance results}
#'   \item{lsmeans}{Least squares means for genotypes with standard errors}
#' 
#' @details This function performs a comprehensive statistical analysis pipeline:
#' 
#' \enumerate{
#'   \item \strong{Normality Test}: Uses Shapiro-Wilk test to assess normality of each trait
#'   \item \strong{Model Selection}: Automatically generates appropriate statistical models based on data structure
#'   \item \strong{ANOVA}: Performs analysis of variance to test significance of effects
#'   \item \strong{Least Squares Means}: Calculates adjusted means for genotypes
#' }
#' 
#' \strong{Automatic Model Generation:}
#' The function automatically generates appropriate formulas based on the experimental design:
#' \itemize{
#'   \item Multiple environments and blocks: \code{trait ~ E + G + E:G + B\%in\%E}
#'   \item Multiple environments only: \code{trait ~ E*G}
#'   \item Custom formulas: User-defined formulas when \code{defineForm} is provided
#' }
#' 
#' \strong{Data Requirements:}
#' The phenotype data must have the following structure:
#' \itemize{
#'   \item Column 1: Environment (E) - factor variable
#'   \item Column 2: Block (B) - factor variable  
#'   \item Column 3: Repetition (R) - factor variable
#'   \item Column 4: Genotype (G) - factor variable
#'   \item Columns 5+: Trait measurements - numeric variables
#' }
#' 
#' @examples
#' \dontrun{
#' # Example with multiple environments and blocks
#' pheno_data <- data.frame(
#'   E = rep(c("Env1", "Env2", "Env3"), each = 60),
#'   B = rep(c("Block1", "Block2"), each = 30, times = 3),
#'   R = rep(1:5, 36),
#'   G = rep(1:12, 15),
#'   Height = rnorm(180, 175, 8),
#'   Weight = rnorm(180, 75, 12)
#' )
#' 
#' # Perform statistical analysis with default formulas
#' results <- pheno_stats(pheno_data)
#' 
#' # View normality test results
#' results$Height$normality_test
#' 
#' # View ANOVA results
#' results$Height$ANOVA
#' 
#' # View least squares means
#' results$Height$lsmeans
#' 
#' # Example with custom formulas
#' custom_formulas <- c(
#'   "Height ~ E + G + E:G + B%in%E",
#'   "Weight ~ E + G + E:G"
#' )
#' results_custom <- pheno_stats(pheno_data, defineForm = custom_formulas)
#' }
#' @importFrom stats shapiro.test lm anova predict as.formula
#' 
#' @references
#' Shapiro, S.S. and Wilk, M.B. (1965). An analysis of variance test for normality. 
#' Biometrika, 52(3-4), 591-611.
#' 
#' @seealso \code{\link{joint_map}} for joint mapping analysis
#' 
#' @export
pheno_stats <-
function(phenoData,defineForm=NULL,effNotation="G"){
  #Environment(E) effect: multiple year and location;  Block(B) effect: block effect nested in every Environment;  Repetition(R):
  if(ncol(phenoData)<5){
    stop("The number of column should be more than 4!")
  }else{
    colnames(phenoData)[1:4] <- c("E","B","R","G")
    phenoData[,"E"] <- as.factor(phenoData[,"E"])
    phenoData[,"B"] <- as.factor(phenoData[,"B"])
    phenoData[,"G"] <- as.factor(phenoData[,"G"])
  }
  Phelist <- list(list())
  #1 normality abalysis for phenotype
  for(n in 1:(ncol(phenoData)-4)){
    normtest <- shapiro.test(phenoData[,4+n])
    if(normtest$"p.value" < 0.05){
      warning(paste("P value of shapiro.test for ",colnames(phenoData)[4+n]," is ",normtest$"p.value",sep=""))
    }
    Phelist[[colnames(phenoData)[4+n]]]["normality_test"] <- list(normtest)
  }
  #2-3 anova and lsmeans
  if(is.null(defineForm)){
    Form <- NULL
    if(nlevels(phenoData[,"E"])>1){
      if(nlevels(phenoData[,"B"])>1){
        for(i in 1:(ncol(phenoData)-4)){
          Form[i] <- as.character(paste(colnames(phenoData)[4+i],"~ E + G + E:G + B%in%E",sep=""))
        }
      }else{
        Form[i] <- as.character(paste(colnames(phenoData)[4+i],"~ E*G",sep=""))
      }
      #}else if(nlevels(phenoData[,"E"]) == 1 && nlevels(phenoData[,"B"])>1){
      #Form[i] <- as.character(paste(colnames(phenoData)[4+i],"~ G + B",sep=""))
    }else{
      stop("Phenotype data is not multiple year-location data!")
    }
    for(j in 1:(ncol(phenoData)-4)){
      lm0 <- lm(as.formula(Form[j]),data=phenoData)
      ano0 <- anova(lm0)
      lms <- as.data.frame(predict(lsmeans::lsmeans(lm0, "G")))
      colnames(lms) <- colnames(phenoData)[4+j]
      Phelist[[colnames(lms)]]["formula"] <- list(defineForm[j])
      Phelist[[colnames(lms)]]["ANOVA"] <- list(ano0)
      Phelist[[colnames(lms)]]["lsmeans"] <- list(lms)
    }
  }else{
    for(j in 1:length(defineForm)){
      lm0 <- lm(as.formula(defineForm[j]),data=phenoData)
      ano0 <- anova(lm0)
      lms <- as.data.frame(predict(lsmeans::lsmeans(lm0,effNotation)))
      #colnames(lms) <- gsub('\\)|\\(', '', defineForm[j][2])
      colnames(lms) <- unlist(strsplit(defineForm[j], "~"))[1]
      Phelist[[colnames(lms)]]["formula"] <- list(defineForm[j])
      Phelist[[colnames(lms)]]["ANOVA"] <- list(ano0)
      Phelist[[colnames(lms)]]["lsmeans"] <- list(lms)
    }
  }
  #4 varcomponent
  return(Phelist)
}
