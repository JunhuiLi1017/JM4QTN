#' Joint Mapping Analysis for Multiple Traits
#' 
#' Performs comprehensive joint mapping analysis for multiple traits using either Association Mapping (AM) 
#' or Linkage Mapping (LM) methods to identify QTL affecting multiple traits simultaneously. This function
#' implements sophisticated statistical methods for multivariate QTL analysis with comprehensive
#' cofactor selection and significance testing.
#' 
#' @param formula A model formula defining the full set of candidate linear-model terms
#'   (for example, factor population effects and marker-by-population interactions).
#'   Term labels from this formula are compared to the selected model in \code{skeleton}.
#' @param data A data frame containing all variables used in \code{formula} and in
#'   \code{skeleton$call$formula}.
#' @param skeleton A fitted model object (typically the final model from stepwise
#'   selection) that contains \code{skeleton$call$formula}. Commonly produced by
#'   \code{skeletion_build()} using the same \code{formula}, \code{data},
#'   \code{include}, and \code{cut_off_list}.
#' @param include Optional character vector of predictor names that should be treated
#'   as grouping variables (coerced to factor when \code{formula} contains \code{:}
#'   interactions). If \code{NULL}, interaction-based inclusion may be inferred from
#'   term labels. Matches the \code{include} argument used when fitting \code{skeleton}.
#' @param cut_off_list A list whose component \code{cut_off} includes a named element
#'   \code{"pvalue"} (for example, output from \code{\link{permutation_test}}).
#'   Used upstream with \code{skeletion_build()} to obtain \code{skeleton}; it is
#'   part of the API for a consistent joint-mapping workflow but is not read inside
#'   \code{joint_map()}.
#' 
#' @return A list containing pvalue and LOD for each term in the formula:
#'   \itemize{
#'     \item \code{pvalue}: P-value for each term
#'     \item \code{lod}: LOD score for each term
#'   }
#' 
#' @examples
#' \donttest{
#' # Example phenotype data
#' set.seed(1)
#' pheno_data <- data.frame(
#'   Trait1 = rnorm(100, mean = 100, sd = 15),
#'   Trait2 = rnorm(100, mean = 50, sd = 8),
#'   Popu = rep(c("Pop1", "Pop2"), each = 50)
#' )
#' 
#' # Example genotype data
#' geno_data <- matrix(sample(c(0,1,2), 100*50, replace = TRUE), 
#'                     nrow = 200, ncol = 50)
#' colnames(geno_data) <- paste0("M", 1:50)
#' 
#' data1 <- cbind(pheno_data,geno_data) 
#' 
#' data1$Popu <- as.factor(data1$Popu)
#' 
#' terms <- c("Popu", paste0(colnames(geno_data), ":Popu"))
#' formula1 <- reformulate(terms, response = "Trait1")
#' 
#' cut_off_list <- permutation_test(formula1, data1, n=100,
#'                      alpha = 0.1)
#' 
#' skeleton <- skeletion_build(
#'   formula1, data1, strategy = "bidirection", metric = "SL",
#'   cut_off_list = cut_off_list, include="Popu"
#' )
#' 
#' results <- joint_map(formula1, data1, skeleton, include = "Popu", cut_off_list = cut_off_list)
#' 
#' print(results)
#' 
#' formula2 <- reformulate(terms, response = "cbind(Trait1,Trait2)")
#' 
#' cut_off_list <- permutation_test(formula2, data1, n = 50, alpha = 0.1)
#' 
#' skeleton <- skeletion_build(
#'   formula2, data1, strategy = "bidirection", metric = "SL",
#'   cut_off_list = cut_off_list, include="Popu"
#' )
#' 
#' results <- joint_map(formula2, data1, skeleton, include = "Popu", cut_off_list = cut_off_list)
#' 
#' print(results)
#'
#' }
#' @seealso \code{\link{permutation_test}} for permutation thresholds,
#'   \code{\link{genotype_prob}} for genotype probability calculations
#' 
#' @references
#' Jiang, C. and Zeng, Z.B. (1995). Multiple trait analysis of genetic mapping for 
#' quantitative trait loci. Genetics, 140(3), 1111-1127.
#' 
#' @importFrom stats anova as.formula df.residual lm resid terms update
#' 
#' @export
#' 

joint_map <- function(formula, data, skeleton, include, cut_off_list){
  nobs <- nrow(data)
  call_best <- skeleton
  formula_best <- skeleton$call$formula
  term_form_best <- terms(formula_best, data = data)
  x_name_best <- attr(term_form_best, "term.labels")

  resd_best <- resid(call_best)
  mat_rss_best <- t(resd_best) %*% resd_best
  det_rss_best <- base::det(mat_rss_best)

  x_name <- attr(terms(formula), "term.labels")
  inlude_var <- include
  continuous_check <- grepl(":", x_name)
  if(any(continuous_check)) {
    cont_class_var <- x_name[continuous_check]
    back_var <- gsub(":.*", "", cont_class_var)
    inlude_var <- unique(back_var)
    data[,inlude_var] <- as.factor(data[,inlude_var])
    x_name_exclude <- x_name[!x_name %in% inlude_var]
  }else{
    x_name_exclude <- x_name
  }
  p_value <- seq_along(x_name_exclude)
  names(p_value) <- x_name_exclude
  lod <- p_value

  for(i in seq_along(x_name_exclude)){
    x <- x_name_exclude[i]
    if(x %in% x_name_best){
      formula_simple <- update(formula_best, as.formula(paste(". ~ . -", x)))
      call_simple <- lm(formula_simple,data=data)
      resd_simple <- resid(call_simple)
      mat_rss_simple <- t(resd_simple) %*% resd_simple
      det_rss_simple <- base::det(mat_rss_simple)
      p_value[i] <- anova(call_simple, call_best)[2, "Pr(>F)"]
      lod[i] <- 0.5 * nobs * log10(det_rss_simple / det_rss_best)
    }else{
      formula_full <- update(formula_best, as.formula(paste(". ~ . +", x)))
      call_full <- lm(formula_full,data=data)
      resd_full <- resid(call_full)
      mat_rss_full <- t(resd_full) %*% resd_full
      det_rss_full <- base::det(mat_rss_full)
      df_full <- df.residual(call_full)
      if(df_full == 0){
        p_value[i] <- 1
        lod[i] <- 0
      }else{
        p_value[i] <- anova(call_best, call_full)[2, "Pr(>F)"]
        lod[i] <- 0.5 * nobs * log10(det_rss_best / det_rss_full)
      }
    }
  }
  return(list(p_value = p_value, lod = lod))
}
