#' Permutation Test for Stepwise Regression
#' 
#' Performs permutation tests for stepwise regression to determine empirical significance thresholds
#' for QTL detection using stepwise regression. This function implements a comprehensive permutation
#' testing framework for stepwise regression.
#' 
#' @param formula A model formula defining the response and candidate predictors for
#'   permutation testing.
#' @param data A data frame containing all variables referenced in \code{formula},
#'   including response trait(s) and predictor variables.
#' @param n Integer number of permutations to run. Larger values provide more stable
#'   empirical thresholds but increase computation time. Default is \code{1000}.
#' @param alpha Numeric significance level used to extract empirical cutoff values
#'   from the permutation distributions. Must be between 0 and 1. Default is \code{0.1}.
#' @param type Type of regression model.
#'   Typical values include \code{"linear"}, \code{"logistic"}, and \code{"poisson"}.
#'   Default is \code{"linear"}.
#' @param include Optional character vector of variable names that should always be
#'   included during stepwise model selection. If \code{NULL}, inclusion terms can be
#'   inferred from interaction terms in \code{formula}. Default is \code{NULL}.
#' @param strategy Stepwise selection strategy passed to \code{StepReg::stepwise()}.
#'   Typical values include \code{"forward"}, \code{"backward"}, and \code{"bidirection"}.
#'   Default is \code{"bidirection"}.
#' @param metric Model selection metric used inside stepwise regression.
#'   Typical values include \code{"AIC"}, \code{"BIC"}, or \code{"SBC"}.
#'   Default is \code{"SBC"}.
#' @return A list containing the empirical significance thresholds for p-value and LOD score with:
#'   \itemize{
#'     \item \strong{cut_off}: A vector containing the empirical significance thresholds for p-value and LOD score
#'     \item \strong{pvalue}: A vector containing the p-values for each permutation
#'     \item \strong{lod}: A vector containing the LOD scores for each permutation
#'   }
#' 
#' @details This function implements a comprehensive permutation testing framework.
#' 
#' \strong{Output Interpretation:}
#' \itemize{
#'   \item \strong{P-values}: Empirical significance thresholds for p-value
#'   \item \strong{LOD scores}: Empirical significance thresholds for LOD score
#' }
#' 
#' @examples
#' \dontrun{
#' # Example phenotype data
#' pheno_data <- data.frame(
#'   Trait1 = rnorm(200, mean = 100, sd = 15),
#'   Trait2 = rnorm(200, mean = 50, sd = 8),
#'   Popu = rep(c("Pop1", "Pop2"), each = 100)
#' )
#' 
#' # Example genotype data
#' geno_data <- matrix(sample(c(0,1,2), 100*50, replace = TRUE), 
#'                     nrow = 200, ncol = 50)
#' colnames(geno_data) <- paste0("M", 1:50)
#' 
#' data1 <- cbind(pheno_data,geno_data)
#' 
#' terms <- c("Popu", paste0(colnames(geno_data), ":Popu"))
#' formula1 <- reformulate(terms, response = "Trait1")
#' 
#' cut_off_list <- permutation_test(formula1, data1, n=100,
#'                      alpha = 0.1)
#' 
#' formula2 <- reformulate(terms, response = "cbind(Trait1,Trait2)")
#' 
#' cut_off_list <- permutation_test(formula1, data1, n=100,
#'                      alpha = 0.1)
#' }
#' 
#' @importFrom stats resid anova update terms lm det df.residual log10
#' @importFrom StepReg stepwise
#' 
#' @export
permutation_test <- function(formula, data, n=1000, alpha=0.1, include=NULL, strategy="bidirection", metric="SBC", type="linear"){
  
  ## get dependent variable name
  term_form <- terms(formula, data = data)
  vars <- as.character(attr(term_form, "variables"))[ -1 ]
  y_name <- vars[attr(term_form, "response")]
  
  if (startsWith(y_name, "cbind(")) {
    inner_content <- gsub("^cbind\\(|\\)$", "", y_name)
    y_var <- trimws(strsplit(inner_content, ",")[[1]])
  } else {
    y_var <- y_name
  }
  ## get independent variable name
  inlude_var <- include
  x_name <- attr(term_form, "term.labels")
  continuous_check <- grepl(":", x_name)
  if(any(continuous_check)) {
    cont_class_var <- x_name[continuous_check]
    front_var <- gsub(".*:", "", cont_class_var)
    back_var <- gsub(":.*", "", cont_class_var)
    inlude_var <- unique(back_var)
    data[,inlude_var] <- as.factor(data[,inlude_var])
    x_name_exclude <- x_name[!x_name %in% inlude_var]
  }else{
    x_name_exclude <- x_name
  }

  p_value_pt <- 1:n
  lod_pt <- p_value_pt
  data_pt <- data
  nobs <- nrow(data)
  for(v in 1:n){
    data_pt[1:nobs,y_var] <- data[sample(1:nobs,nobs,replace=FALSE),y_var]
    
    stepwise_var <- stepwise(formula, data_pt, type="linear", strategy = strategy, metric=metric, include=inlude_var)
    
    call_best <- stepwise_var[[strategy]][[metric]]
    formula_best <- call_best$call$formula
    
    term_form_best <- terms(call_best, data = data_pt)
    x_name_best <- attr(term_form_best, "term.labels")
    resd_best <- resid(call_best)
    mat_rss_best <- t(resd_best) %*% resd_best
    det_rss_best <- det(mat_rss_best) 
    df_best <- df.residual(call_best)
    
    p_value <- seq_along(x_name_exclude)
    lod <- p_value

    for(i in seq_along(x_name_exclude)){
      x <- x_name_exclude[i]
      if(x %in% x_name_best){
        formula_simple <- update(formula_best, as.formula(paste(". ~ . -", x)))
        call_simple <- lm(formula_simple,data=data_pt)
        resd_simple <- resid(call_simple)
        mat_rss_simple <- t(resd_simple) %*% resd_simple
        det_rss_simple <- det(mat_rss_simple)
        p_value[i] <- anova(call_best,call_simple)[2,"Pr(>F)"]
        lod[i] <- 0.5*nobs*log10(det_rss_best/det_rss_simple)
      }else{
        formula_full <- update(formula_best, as.formula(paste(". ~ . +", x)))
        call_full <- lm(formula_full,data=data_pt)
        resd_full <- resid(call_full)
        mat_rss_full <- t(resd_full) %*% resd_full
        det_rss_full <- det(mat_rss_full)
        df_full <- df.residual(call_full)
        if(df_full == 0){
          p_value[i] <- 1
          lod[i] <- 0
        }else{
          p_value[i] <- anova(call_full,call_simple)[2,"Pr(>F)"]
          lod[i] <- 0.5*nobs*log10(det_rss_full/det_rss_simple)
        }
      }
    }
    p_value_pt[v] <- min(p_value, na.rm = T)
    lod_pt[v] <- max(lod, na.rm = T)
    if(v %% (n/10)==0){
      cat(v/n*100,"% have been completed\t\t","@",as.character(Sys.time()),"\n")
    }
  }
  p_value_sort <- sort(p_value_pt)
  lod_sort <- sort(lod_pt, decreasing = T)
  cut_off <- c(p_value_sort[round(n*alpha)], lod_sort[round(n*alpha)])
  names(cut_off) <- c("pvalue","lod")
  return(list(cut_off = cut_off, pvalue = p_value_sort, lod=lod_sort))
}
