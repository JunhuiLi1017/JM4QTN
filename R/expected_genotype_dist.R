#' Expected Allele Genotype Distribution Calculator
#' 
#' Calculates expected allele genotype distribution probabilities for different marker types
#' and cross types in genetic mapping studies. This function implements sophisticated algorithms
#' to compute genotype probabilities based on flanking marker information and population structure.
#' 
#' @param marType Character string specifying the marker type. Options include:
#'   \itemize{
#'     \item \code{"22", "21", "20", "12", "11", "10", "02", "01", "00"}: Both flanking markers observed
#'     \item \code{"N2", "N1", "N0"}: Only right flanking marker observed
#'     \item \code{"2N", "1N", "0N"}: Only left flanking marker observed
#'   }
#'   These codes represent different combinations of flanking marker genotypes.
#' @param croType Character string specifying the cross type:
#'   \itemize{
#'     \item \code{"Fn"}: F-n generation populations
#'     \item \code{"BCP1"}: Backcross to parent 1
#'     \item \code{"BCP2"}: Backcross to parent 2
#'     \item \code{"F2"}: F2 generation
#'     \item \code{"DH"}: Doubled haploid
#'     \item \code{"RIL"}: Recombinant inbred line
#'   }
#' @param Gn Numeric value specifying the generation number. Must be greater than 0.
#' @param x Numeric value representing recombination fraction between loci A and Q.
#'   Must be between 0 and 0.5.
#' @param y Numeric value representing recombination fraction between loci Q and B.
#'   Must be between 0 and 0.5. Default is 0.
#' 
#' @return A numeric value representing the expected probability of the specified genotype.
#'   The result is typically between -1 and 1, where positive values indicate higher
#'   probability of the target genotype.
#' 
#' @details This function calculates expected genotype probabilities using different approaches:
#' 
#' \strong{Marker Type Classification:}
#' \itemize{
#'   \item \strong{Double flanking markers} (\code{"22", "21", "20", etc.}): Both left and right
#'         flanking markers are observed, providing maximum information
#'   \item \strong{Single right flanking marker} (\code{"N2", "N1", "N0"}): Only the right
#'         flanking marker is observed
#'   \item \strong{Single left flanking marker} (\code{"2N", "1N", "0N"}): Only the left
#'         flanking marker is observed
#' }
#' 
#' \strong{Calculation Methods by Cross Type:}
#' \itemize{
#'   \item \strong{Fn populations}: Uses \code{\link{genotype_freq}} for complex calculations
#'         involving multiple genotype classes
#'   \item \strong{F2, DH, RIL}: Uses analytical formulas based on classical genetics
#'   \item \strong{BCP1, BCP2}: Uses \code{\link{genotype_freq}} with specific genotype indices
#' }
#' 
#' \strong{Mathematical Framework:}
#' The function uses different probability calculations depending on the available information:
#' \itemize{
#'   \item When both flanking markers are observed, it uses conditional probabilities
#'   \item When only one flanking marker is observed, it uses marginal probabilities
#'   \item The calculations incorporate recombination fractions and population-specific parameters
#' }
#' 
#' \strong{Genotype Coding:}
#' The marker types use a coding system where:
#' \itemize{
#'   \item \code{2}: Homozygous for alternative allele
#'   \item \code{1}: Heterozygous
#'   \item \code{0}: Homozygous for reference allele
#'   \item \code{N}: Missing or unknown genotype
#' }
#' 
#' @examples
#' \dontrun{
#' # Calculate probability for F2 population with both flanking markers
#' prob_f2 <- expected_genotype_dist("22", "F2", Gn = 2, x = 0.1, y = 0.2)
#' 
#' # Calculate probability for BCP1 with only right flanking marker
#' prob_bcp1_n2 <- expected_genotype_dist("N2", "BCP1", Gn = 3, x = 0.15, y = 0.25)
#' 
#' # Calculate probability for DH population with only left flanking marker
#' prob_dh_2n <- expected_genotype_dist("2N", "DH", Gn = 2, x = 0.1, y = 0)
#' 
#' # Calculate probability for Fn population with both flanking markers
#' prob_fn <- expected_genotype_dist("21", "Fn", Gn = 4, x = 0.2, y = 0.3)
#' 
#' # Calculate probability for RIL population
#' prob_ril <- expected_genotype_dist("00", "RIL", Gn = 2, x = 0.1, y = 0.2)
#' 
#' # Example with different recombination fractions
#' prob_low_rec <- expected_genotype_dist("22", "F2", Gn = 2, x = 0.05, y = 0.05)
#' prob_high_rec <- expected_genotype_dist("22", "F2", Gn = 2, x = 0.3, y = 0.4)
#' }
#' 
#' @aliases calculate_expected_genotype_distribution
#' @seealso \code{\link{genotype_freq}} for genotype frequency calculations,
#'          \code{\link{genotype_prob}} for missing genotype imputation
#' 
#' @references
#' Haldane, J.B.S. (1919). The combination of linkage values and the calculation of 
#' distances between the loci of linked factors. Journal of Genetics, 8(3), 299-309.
#' 
#' @export
expected_genotype_dist <-
function(marType,croType,Gn=2,x,y=0){   
  if(Gn < 1){
    stop("Gn should > 0")
  }
  z <- x + y -2*x*y
  if(marType == "22"){
    if (croType=="BCP1"){
      rval <- genotype_freq(croType,Gn+1,8,x,y)/(genotype_freq(croType,Gn+1,8,x,y)+genotype_freq(croType,Gn+1,6,x,y))
    }
    if (croType=="DH"){
      rval <- (1-x-y)/(1-z)
    }
    if (croType=="F2"|croType=="RIL"){
      rval <- ((1-x)*(1-x)*(1-y)*(1-y) - x*x*y*y)/((1-z)*(1-z))
    }
    if (croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,1,x,y) - calculate_genotype_frequencies(croType,Gn,3,x,y))/(calculate_genotype_frequencies(croType,Gn,1,x,y) + calculate_genotype_frequencies(croType,Gn,6,x,y) + calculate_genotype_frequencies(croType,Gn,3,x,y))
    }
  }else if(marType == "21"){
    if (croType=="BCP1"){
      rval <- calculate_genotype_frequencies(croType,Gn+1,7,x,y)/(calculate_genotype_frequencies(croType,Gn+1,7,x,y)+calculate_genotype_frequencies(croType,Gn+1,5,x,y))
    }
    if (croType=="F2"){
      rval <- ((1-x)*(1-x)*(1-y)*y - x*x*y*(1-y))/((1-z)*z)
    }
    if (croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,5,x,y) - calculate_genotype_frequencies(croType,Gn,15,x,y))/(calculate_genotype_frequencies(croType,Gn,5,x,y) + calculate_genotype_frequencies(croType,Gn,7,x,y) + calculate_genotype_frequencies(croType,Gn,11,x,y) + calculate_genotype_frequencies(croType,Gn,15,x,y))
    }
  }else if(marType == "20"){
    if (croType=="DH"){
      rval <- (y-x)/z
    }
    if (croType=="F2"|croType=="RIL"){
      rval <- ((1-x)*(1-x)*y*y - x*x*(1-y)*(1-y))/(z*z)
    }
    if (croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,2,x,y) - calculate_genotype_frequencies(croType,Gn,4,x,y))/(calculate_genotype_frequencies(croType,Gn,2,x,y) + calculate_genotype_frequencies(croType,Gn,12,x,y) + calculate_genotype_frequencies(croType,Gn,4,x,y))
    }
  }else if(marType == "12"){
    if (croType=="BCP1"){
      rval <- calculate_genotype_frequencies(croType,Gn+1,4,x,y)/(calculate_genotype_frequencies(croType,Gn+1,4,x,y)+calculate_genotype_frequencies(croType,Gn+1,2,x,y))
    }
    if (croType=="F2"){
      rval <- ((1-y)*(1-y)*x*(1-x) - x*(1-x)*y*y)/(z*(1-z))
    }
    if (croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,8,x,y) - calculate_genotype_frequencies(croType,Gn,14,x,y))/(calculate_genotype_frequencies(croType,Gn,8,x,y) + calculate_genotype_frequencies(croType,Gn,10,x,y) + calculate_genotype_frequencies(croType,Gn,16,x,y) + calculate_genotype_frequencies(croType,Gn,14,x,y))
    }
  }else if(marType == "11"){
    if (croType=="BCP1"){
      rval <- calculate_genotype_frequencies(croType,Gn+1,3,x,y)/(calculate_genotype_frequencies(croType,Gn+1,3,x,y)+calculate_genotype_frequencies(croType,Gn+1,1,x,y))
    }
    if (croType=="BCP2"){
      rval <- -(calculate_genotype_frequencies(croType,Gn+1,6,x,y)/(calculate_genotype_frequencies(croType,Gn+1,8,x,y)+calculate_genotype_frequencies(croType,Gn+1,6,x,y)))
    }
    if (croType=="F2"|croType=="Fn"){
      rval <- 0
    }
  }else if(marType == "10"){
    if (croType=="BCP2"){
      rval <- -(calculate_genotype_frequencies(croType,Gn+1,5,x,y)/(calculate_genotype_frequencies(croType,Gn+1,7,x,y)+calculate_genotype_frequencies(croType,Gn+1,5,x,y)))
    }
    if (croType=="F2"){
      rval <- (x*(1-x)*y*y -x*(1-x)*(1-y)*(1-y))/(z*(1-z))
    }
    if (croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,14,x,y) - calculate_genotype_frequencies(croType,Gn,8,x,y))/(calculate_genotype_frequencies(croType,Gn,14,x,y) + calculate_genotype_frequencies(croType,Gn,10,x,y) + calculate_genotype_frequencies(croType,Gn,16,x,y) + calculate_genotype_frequencies(croType,Gn,8,x,y))
    }
  }else if(marType == "02"){
    if (croType=="DH"){
      rval <- (x-y)/z
    }
    if (croType=="F2"|croType=="RIL"){
      rval <- ((1-y)*(1-y)*x*x - (1-x)*(1-x)*y*y)/(z*z)
    }
    if (croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,4,x,y) - calculate_genotype_frequencies(croType,Gn,2,x,y))/(calculate_genotype_frequencies(croType,Gn,4,x,y) + calculate_genotype_frequencies(croType,Gn,12,x,y) + calculate_genotype_frequencies(croType,Gn,2,x,y))
    }
  }else if(marType == "01"){
    if (croType=="BCP2"){
      rval <- -(calculate_genotype_frequencies(croType,Gn+1,2,x,y)/(calculate_genotype_frequencies(croType,Gn+1,4,x,y)+calculate_genotype_frequencies(croType,Gn+1,2,x,y)))
    }
    if (croType=="F2"){
      rval <- (x*x*(1-y)*y - (1-x)*(1-x)*(1-y)*y)/((1-z)*z)
    }
    if (croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,15,x,y) - calculate_genotype_frequencies(croType,Gn,5,x,y))/(calculate_genotype_frequencies(croType,Gn,5,x,y) + calculate_genotype_frequencies(croType,Gn,7,x,y) + calculate_genotype_frequencies(croType,Gn,11,x,y) + calculate_genotype_frequencies(croType,Gn,15,x,y))
    }
  }else if(marType == "00"){
    if (croType=="BCP2"){
      rval <- -(calculate_genotype_frequencies(croType,Gn+1,1,x,y)/(calculate_genotype_frequencies(croType,Gn+1,3,x,y)+calculate_genotype_frequencies(croType,Gn+1,1,x,y)))
    }
    if (croType=="DH"){
      rval <- (x+y-1)/(1-z)
    }
    if (croType=="F2"|croType=="RIL"){
      rval <- (x*x*y*y - (1-x)*(1-x)*(1-y)*(1-y))/((1-z)*(1-z))
    }
    if (croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,3,x,y) - calculate_genotype_frequencies(croType,Gn,1,x,y))/(calculate_genotype_frequencies(croType,Gn,1,x,y) + calculate_genotype_frequencies(croType,Gn,6,x,y) + calculate_genotype_frequencies(croType,Gn,3,x,y))
    }
  }else if(marType == "2N"){
    if (croType=="BCP1"){
      rval <- calculate_genotype_frequencies(croType,Gn,4,x)/(calculate_genotype_frequencies(croType,Gn,3,x)+calculate_genotype_frequencies(croType,Gn,4,x))
    }
    if (croType=="DH"){
      rval <- 1-2*x
    }
    if (croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,9,x)-calculate_genotype_frequencies(croType,Gn,7,x))/(calculate_genotype_frequencies(croType,Gn,9,x)+calculate_genotype_frequencies(croType,Gn,8,x)+calculate_genotype_frequencies(croType,Gn,7,x))
    }
    if (croType=="F2"){
      rval <- (1-x)*(1-x) - x*x
    }
    if (croType=="RIL"){
      rval <- ((1-x)*(1-x) - x*x)/((1-x)*(1-x) + x*x)
    }
  }else if(marType == "1N"){
    if (croType=="BCP1"){
      rval <- calculate_genotype_frequencies(croType,Gn,2,x)/(calculate_genotype_frequencies(croType,Gn,1,x)+calculate_genotype_frequencies(croType,Gn,2,x))
    }
    if (croType=="BCP2"){
      rval <- -calculate_genotype_frequencies(croType,Gn,3,x)/(calculate_genotype_frequencies(croType,Gn,4,x)+calculate_genotype_frequencies(croType,Gn,3,x))
    }
    if (croType=="F2"|croType=="Fn"){
      rval <- 0
    }
  }else if(marType == "0N"){
    if (croType=="BCP2"){
      rval <- -calculate_genotype_frequencies(croType,Gn,1,x)/(calculate_genotype_frequencies(croType,Gn,1,x)+calculate_genotype_frequencies(croType,Gn,2,x))
    }
    if (croType=="DH"){
      rval <- 2*x-1
    }
    if(croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,3,x)-calculate_genotype_frequencies(croType,Gn,1,x))/(calculate_genotype_frequencies(croType,Gn,3,x)+calculate_genotype_frequencies(croType,Gn,2,x)+calculate_genotype_frequencies(croType,Gn,1,x))
    }
    if (croType=="F2"){
      rval <- x*x - (1-x)*(1-x)
    }
    if(croType=="RIL"){
      rval <- (x*x - (1-x)*(1-x))/((1-x)*(1-x) + x*x)
    }
  }else if(marType == "N2"){
    if (croType=="BCP1"){
      rval <- calculate_genotype_frequencies(croType,Gn,4,x)/(calculate_genotype_frequencies(croType,Gn,4,x)+calculate_genotype_frequencies(croType,Gn,2,x))
    }
    if (croType=="DH"){
      rval <- 1-2*x
    }
    if (croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,9,x)-calculate_genotype_frequencies(croType,Gn,3,x))/(calculate_genotype_frequencies(croType,Gn,9,x)+calculate_genotype_frequencies(croType,Gn,6,x)+calculate_genotype_frequencies(croType,Gn,3,x))
    }
    if (croType=="F2"){
      rval <- (1-x)*(1-x) - x*x
    }
    if (croType=="RIL"){
      rval <- ((1-x)*(1-x) - x*x)/((1-x)*(1-x) + x*x)
    }
  }else if(marType == "N1"){
    if (croType=="BCP1"){
      rval <- calculate_genotype_frequencies(croType,Gn,3,x)/(calculate_genotype_frequencies(croType,Gn,3,x)+calculate_genotype_frequencies(croType,Gn,1,x))
    }
    if (croType=="BCP2"){
      rval <- -calculate_genotype_frequencies(croType,Gn,2,x)/(calculate_genotype_frequencies(croType,Gn,2,x)+calculate_genotype_frequencies(croType,Gn,4,x))
    }
    if (croType=="F2"|croType=="Fn"){
      rval <- 0
    }
  }else if(marType == "N0"){
    if (croType=="BCP2"){
      rval <- -calculate_genotype_frequencies(croType,Gn,1,x)/(calculate_genotype_frequencies(croType,Gn,1,x)+calculate_genotype_frequencies(croType,Gn,3,x))
    }
    if (croType=="DH"){
      rval <- 2*x-1
    }
    if(croType=="Fn"){
      rval <- (calculate_genotype_frequencies(croType,Gn,7,x)-calculate_genotype_frequencies(croType,Gn,1,x))/(calculate_genotype_frequencies(croType,Gn,1,x)+calculate_genotype_frequencies(croType,Gn,4,x)+calculate_genotype_frequencies(croType,Gn,7,x))
      #error inrval <- (x*x - (1-x)*(1-x))/((1-x)*(1-x) + x*x + (2/(2^(Gn-1)-1))*x*(1-x))
    }
    if (croType=="F2"){
      rval <- x*x - (1-x)*(1-x)
    }
    if(croType=="RIL"){
      rval <- (x*x - (1-x)*(1-x))/((1-x)*(1-x) + x*x)
    }
  }else{
    stop("Wrong mark type for argument marType")
  }
  return(rval)
}

# backward-compatible alias
calculate_expected_genotype_distribution <- expected_genotype_dist
