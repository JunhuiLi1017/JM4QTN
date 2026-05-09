#' Fit the stepwise “skeleton” model for joint mapping
#'
#' Runs \code{StepReg::stepwise()} with entry/stay levels set from a permutation
#' p-value threshold, and returns the selected best linear model. The result is
#' typically passed to \code{\link{joint_map}()} as the \code{skeleton} argument.
#'
#' @param formula A model formula (same as used in \code{\link{permutation_test}}).
#' @param data A data frame containing all variables in \code{formula}.
#' @param type Model type passed to \code{StepReg::stepwise()}. Default is
#'   \code{"linear"}.
#' @param strategy Stepwise strategy, for example \code{"forward"},
#'   \code{"backward"}, or \code{"bidirection"}. Default is \code{"bidirection"}.
#' @param metric Information criterion or selection metric in StepReg, for
#'   example \code{"SL"}. Default is \code{"SL"}.
#' @param include Optional character vector of terms to keep in the stepwise
#'   search, passed to \code{StepReg::stepwise()}. Default is \code{NULL}.
#' @param cut_off_list A list with a component \code{cut_off} that includes a
#'   named entry \code{"pvalue"} (for example, the return value of
#'   \code{\link{permutation_test}}). This value is used for both \code{sle} and
#'   \code{sls} in the stepwise call.
#'
#' @return The best model object for the chosen \code{strategy} and
#'   \code{metric} (an element of the \code{StepReg::stepwise()} result,
#'   typically with a \code{call} and \code{call$formula}).
#'
#' @seealso \code{\link{joint_map}}, \code{\link{permutation_test}}
#' 
#' @importFrom StepReg stepwise
#' 
#' @export
skeletion_build <- function(formula, data, type = "linear", strategy = "bidirection",
                            metric = "SL", include = NULL, cut_off_list) {
  cut_off_pvalue <- cut_off_list$cut_off["pvalue"]
  stepwise_model <- StepReg::stepwise(
    formula, data, type = type, strategy = strategy, metric = metric,
    include = include, sle = cut_off_pvalue, sls = cut_off_pvalue
  )
  best_model <- stepwise_model[[strategy]][[metric]]
  best_model
}
