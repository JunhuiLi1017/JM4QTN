#' Haldane Mapping Function
#' 
#' Converts genetic distance in Morgan units to recombination fraction using the Haldane mapping function.
#' This function implements the Haldane mapping function which assumes no interference between crossovers.
#' 
#' @param x A numeric value representing the genetic distance in Morgan units (centimorgans/100).
#'   Must be non-negative.
#' 
#' @return A numeric value representing the recombination fraction (r) between 0 and 0.5.
#' 
#' @details The Haldane mapping function is defined as:
#' \deqn{r = \frac{1}{2}(1 - e^{-2x})}
#' where \eqn{x} is the genetic distance in Morgan units and \eqn{r} is the recombination fraction.
#' 
#' This function assumes no interference between crossovers, meaning that the occurrence
#' of one crossover does not affect the probability of other crossovers occurring nearby.
#' 
#' @examples
#' # Convert 0.1 Morgan (10 cM) to recombination fraction
#' haldane_map(0.1)
#' 
#' # Convert 0.5 Morgan (50 cM) to recombination fraction  
#' haldane_map(0.5)
#' 
#' # Convert 1.0 Morgan (100 cM) to recombination fraction
#' haldane_map(1.0)
#' 
#' @aliases haldane_mapping_function
#' 
#' @references
#' Haldane, J.B.S. (1919). The combination of linkage values and the calculation of 
#' distances between the loci of linked factors. Journal of Genetics, 8(3), 299-309.
#' 
#' @seealso \code{\link{genotype_prob}} for genotype probability calculations
#' 
#' @export
haldane_map <-
function(x){
  r <- 0.5*(1-exp(-2*x))
  return(r)
}

