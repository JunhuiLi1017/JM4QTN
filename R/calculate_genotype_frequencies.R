#' Genotype Frequency Calculator for Genetic Populations
#' 
#' Calculates genotype frequencies for different cross types and generations in genetic mapping studies.
#' This function implements complex recurrence relations to compute genotype frequencies in various
#' genetic populations including Fn, backcross, and advanced generation populations.
#' 
#' @param CrosType Character string specifying the cross type:
#'   \itemize{
#'     \item \code{"Fn"}: F-n generation populations (n > 2)
#'     \item \code{"BCP1"}: Backcross to parent 1 populations
#'     \item \code{"BCP2"}: Backcross to parent 2 populations
#'     \item \code{"F2"}: F2 generation populations
#'     \item \code{"DH"}: Doubled haploid populations
#'     \item \code{"RIL"}: Recombinant inbred line populations
#'   }
#' @param Gn Numeric value specifying the generation number. Must be greater than 1.
#' @param i Numeric value specifying the genotype index to calculate frequency for.
#'   The valid range depends on the cross type:
#'   \itemize{
#'     \item \code{"Fn"}: 1-20 (20 different genotype classes)
#'     \item \code{"BCP1/BCP2"}: 1-8 (8 genotype classes)
#'     \item \code{"F2"}: 1-9 (9 genotype classes)
#'   }
#' @param x Numeric value representing recombination fraction between loci A and Q.
#'   Must be between 0 and 0.5.
#' @param y Numeric value representing recombination fraction between loci Q and B.
#'   Must be between 0 and 0.5. Default is 0.
#' 
#' @return A numeric value representing the frequency of the specified genotype in the given generation.
#'   The result is always between 0 and 1.
#' 
#' @details This function calculates genotype frequencies using sophisticated mathematical models:
#' 
#' \strong{For Fn Populations (y > 0):}
#' \itemize{
#'   \item Handles 20 different genotype classes with complex recurrence relations
#'   \item Accounts for recombination between three loci (A, Q, B)
#'   \item Uses the relationship \eqn{z = x + y - 2xy} for flanking marker recombination
#'   \item Implements generation-by-generation frequency calculations
#' }
#' 
#' \strong{For BCP1 Populations:}
#' \itemize{
#'   \item Handles 8 genotype classes
#'   \item Backcross to the first parent (AAQQBB)
#'   \item Specific recurrence relations for each genotype class
#' }
#' 
#' \strong{For BCP2 Populations:}
#' \itemize{
#'   \item Handles 8 genotype classes
#'   \item Backcross to the second parent (aaqqbb)
#'   \item Different initial conditions and recurrence relations
#' }
#' 
#' \strong{For F2 Populations (y = 0):}
#' \itemize{
#'   \item Handles 9 genotype classes
#'   \item Simplified two-locus model
#'   \item Standard F2 generation frequencies
#' }
#' 
#' \strong{Mathematical Foundation:}
#' The function uses recurrence relations to calculate genotype frequencies across generations.
#' For each generation, the frequency of each genotype class is computed based on the
#' frequencies in the previous generation and the recombination fractions between loci.
#' 
#' @examples
#' # Calculate frequency for F2 generation, genotype 1, with recombination fraction 0.1
#' freq_f2 <- calculate_genotype_frequencies("Fn", Gn = 2, i = 1, x = 0.1, y = 0)
#' 
#' # Calculate frequency for BCP1 generation, genotype 5, with recombination fractions
#' freq_bcp1 <- calculate_genotype_frequencies("BCP1", Gn = 3, i = 5, x = 0.15, y = 0.25)
#' 
#' # Calculate frequency for Fn generation, genotype 10, with recombination fractions
#' freq_fn <- calculate_genotype_frequencies("Fn", Gn = 4, i = 10, x = 0.2, y = 0.3)
#' 
#' # Calculate frequency for BCP2 generation, genotype 3
#' freq_bcp2 <- calculate_genotype_frequencies("BCP2", Gn = 2, i = 3, x = 0.1, y = 0.2)
#' 
#' # Example with different generation numbers
#' freq_gen3 <- calculate_genotype_frequencies("Fn", Gn = 3, i = 1, x = 0.1, y = 0)
#' freq_gen5 <- calculate_genotype_frequencies("Fn", Gn = 5, i = 1, x = 0.1, y = 0)
#' 
#' @seealso \code{\link{calculate_expected_genotype_distribution}} for expected allele genotype distribution calculations
#' 
#' @references
#' Haldane, J.B.S. (1919). The combination of linkage values and the calculation of 
#' distances between the loci of linked factors. Journal of Genetics, 8(3), 299-309.
#' 
#' @export
calculate_genotype_frequencies <-
function(CrosType,Gn=2,i,x,y=0){
  if(Gn<=1){
    stop(" Gn should > 0")
  }
  if(y != 0){
    z <- x + y -2*x*y
    if(CrosType == "Fn"){
      mtGeFre <- matrix(data=c(0), nrow = 50, ncol = 20)	##mtGeFre: Matrix of Frequency of the Genotypes ,where nrow is Gn and ncol is 1~20.
      mtGeFre[1,20] <- 1	 ##mtGeFre[1,20]=1 means the Frequency of F1 which comes frome AAQQBB * aaqqbb -> AaQqBb equals to 1.
      for (g in 2:Gn){
        #AQB/AQB or aqb/aqb
        mtGeFre[g,1] <- 0.25*(4*mtGeFre[g-1,1]+mtGeFre[g-1,5]+mtGeFre[g-1,6]+(1-y)*(1-y)*mtGeFre[g-1,7]+mtGeFre[g-1,8]+(1-z)*(1-z)*mtGeFre[g-1,9]+
                                (1-x)*(1-x)*mtGeFre[g-1,10]+y*y*mtGeFre[g-1,11]+z*z*mtGeFre[g-1,13]+x*x*mtGeFre[g-1,16]+ x*x*(1-y)*(1-y)*mtGeFre[g-1,17]+
                                x*x*y*y*mtGeFre[g-1,18]+y*y*(1-x)*(1-x)*mtGeFre[g-1,19]+(1-x)*(1-x)*(1-y)*(1-y)*mtGeFre[g-1,20])
        #AQb/AQb or aqB/aqB
        mtGeFre[g,2] <- 0.25*(4*mtGeFre[g-1,2]+mtGeFre[g-1,5]+y*y*mtGeFre[g-1,7]+z*z*mtGeFre[g-1,9]+(1-x)*(1-x)*mtGeFre[g-1,10]+(1-y)*(1-y)*mtGeFre[g-1,11]+
                                mtGeFre[g-1,12]+(1-z)*(1-z)*mtGeFre[g-1,13]+mtGeFre[g-1,14]+x*x*mtGeFre[g-1,16]+ x*x*y*y*mtGeFre[g-1,17]+
                                x*x*(1-y)*(1-y)*mtGeFre[g-1,18]+(1-y)*(1-y)*(1-x)*(1-x)*mtGeFre[g-1,19]+(1-x)*(1-x)*y*y*mtGeFre[g-1,20])
        #AqB/AqB or aQb/aQb
        mtGeFre[g,3] <- 0.25*(4*mtGeFre[g-1,3]+mtGeFre[g-1,6]+y*y*mtGeFre[g-1,7]+(1-z)*(1-z)*mtGeFre[g-1,9]+x*x*mtGeFre[g-1,10]+
                                (1-y)*(1-y)*mtGeFre[g-1,11]+z*z*mtGeFre[g-1,13]+mtGeFre[g-1,14]+mtGeFre[g-1,15]+(1-x)*(1-x)*mtGeFre[g-1,16]+ (1-x)*(1-x)*y*y*mtGeFre[g-1,17]+
                                (1-x)*(1-x)*(1-y)*(1-y)*mtGeFre[g-1,18]+(1-y)*(1-y)*x*x*mtGeFre[g-1,19]+x*x*y*y*mtGeFre[g-1,20])
        #aQB/aQB or Aqb/Aqb
        mtGeFre[g,4] <- 0.25*(4*mtGeFre[g-1,4]+(1-y)*(1-y)*mtGeFre[g-1,7]+mtGeFre[g-1,8]+z*z*mtGeFre[g-1,9]+x*x*mtGeFre[g-1,10]+
                                y*y*mtGeFre[g-1,11]+mtGeFre[g-1,12]+(1-z)*(1-z)*mtGeFre[g-1,13]+mtGeFre[g-1,15]+(1-x)*(1-x)*mtGeFre[g-1,16]+ (1-x)*(1-x)*(1-y)*(1-y)*mtGeFre[g-1,17]+
                                (1-x)*(1-x)*y*y*mtGeFre[g-1,18]+y*y*x*x*mtGeFre[g-1,19]+x*x*(1-y)*(1-y)*mtGeFre[g-1,20])
        #aqB/aqb or AQB/AQb
        mtGeFre[g,5] <- 0.5*(mtGeFre[g-1,5]+y*(1-y)*mtGeFre[g-1,7]+z*(1-z)*mtGeFre[g-1,9]+y*(1-y)*mtGeFre[g-1,11]+z*(1-z)*mtGeFre[g-1,13]+
                               x*x*y*(1-y)*mtGeFre[g-1,17]+x*x*y*(1-y)*mtGeFre[g-1,18]+y*(1-y)*(1-x)*(1-x)*mtGeFre[g-1,19]+(1-x)*(1-x)*y*(1-y)*mtGeFre[g-1,20])
        #aqb/aQb or AqB/AQB
        mtGeFre[g,6] <- 0.5*(mtGeFre[g-1,6]+y*(1-y)*mtGeFre[g-1,7]+x*(1-x)*mtGeFre[g-1,10]+y*(1-y)*mtGeFre[g-1,11]+x*(1-x)*mtGeFre[g-1,16]+
                               x*y*(1-x)*(1-y)*mtGeFre[g-1,17]+x*y*(1-x)*(1-y)*mtGeFre[g-1,18]+x*y*(1-x)*(1-y)*mtGeFre[g-1,19]+x*y*(1-x)*(1-y)*mtGeFre[g-1,20])
        #aQB/aqb or Aqb/AQB
        mtGeFre[g,7] <- 0.5*((1-y)*(1-y)*mtGeFre[g-1,7]+y*y*mtGeFre[g-1,11]+ x*(1-x)*(1-y)*(1-y)*mtGeFre[g-1,17]+x*y*y*(1-x)*mtGeFre[g-1,18]+
                               x*y*y*(1-x)*mtGeFre[g-1,19]+x*(1-x)*(1-y)*(1-y)*mtGeFre[g-1,20])
        #Aqb/aqb or aQB/AQB
        mtGeFre[g,8] <- 0.5*(mtGeFre[g-1,8]+z*(1-z)*mtGeFre[g-1,9]+x*(1-x)*mtGeFre[g-1,10]+z*(1-z)*mtGeFre[g-1,13]+x*(1-x)*mtGeFre[g-1,16]+
                               x*(1-x)*(1-y)*(1-y)*mtGeFre[g-1,17]+x*y*y*(1-x)*mtGeFre[g-1,18]+x*y*y*(1-x)*mtGeFre[g-1,19]+x*(1-x)*(1-y)*(1-y)*mtGeFre[g-1,20])
        #AqB/aqb or aQb/AQB
        mtGeFre[g,9] <- 0.5*((1-z)*(1-z)*mtGeFre[g-1,9]+z*z*mtGeFre[g-1,13]+ x*(1-x)*(1-y)*y*mtGeFre[g-1,17]+x*y*(1-y)*(1-x)*mtGeFre[g-1,18]+
                               x*y*(1-y)*(1-x)*mtGeFre[g-1,19]+x*(1-x)*(1-y)*y*mtGeFre[g-1,20])
        #AQb/aqb or aqB/AQB
        mtGeFre[g,10] <- 0.5*((1-x)*(1-x)*mtGeFre[g-1,10]+x*x*mtGeFre[g-1,16]+ x*x*(1-y)*y*mtGeFre[g-1,17]+x*y*(1-y)*x*mtGeFre[g-1,18]+
                                (1-x)*y*(1-y)*(1-x)*mtGeFre[g-1,19]+(1-x)*(1-x)*(1-y)*y*mtGeFre[g-1,20])
        #aqB/aQb or AqB/AQb
        mtGeFre[g,11] <- 0.5*(y*y*mtGeFre[g-1,7]+(1-y)*(1-y)*mtGeFre[g-1,11]+x*(1-x)*y*y*mtGeFre[g-1,17]+x*(1-y)*(1-y)*(1-x)*mtGeFre[g-1,18]+
                                x*(1-y)*(1-y)*(1-x)*mtGeFre[g-1,19]+x*y*(1-x)*y*mtGeFre[g-1,20])
        #aqB/aQB or Aqb/AQb
        mtGeFre[g,12] <- 0.5*(y*(1-y)*mtGeFre[g-1,7]+x*(1-x)*mtGeFre[g-1,10]+y*(1-y)*mtGeFre[g-1,11]+mtGeFre[g-1,12]+x*(1-x)*mtGeFre[g-1,16]+
                                x*(1-x)*y*(1-y)*mtGeFre[g-1,17]+x*(1-x)*y*(1-y)*mtGeFre[g-1,18]+x*(1-x)*y*(1-y)*mtGeFre[g-1,19]+x*(1-x)*y*(1-y)*mtGeFre[g-1,20])
        #aqB/Aqb or aQB/AQb
        mtGeFre[g,13] <- 0.5*(z*z*mtGeFre[g-1,9]+(1-z)*(1-z)*mtGeFre[g-1,13]+ x*(1-x)*(1-y)*y*mtGeFre[g-1,17]+x*y*(1-y)*(1-x)*mtGeFre[g-1,18]+
                                x*y*(1-y)*(1-x)*mtGeFre[g-1,19]+x*(1-x)*(1-y)*y*mtGeFre[g-1,20])
        #aqB/AqB or aQb/AQb
        mtGeFre[g,14] <- 0.5*(z*(1-z)*mtGeFre[g-1,9]+x*(1-x)*mtGeFre[g-1,10]+z*(1-z)*mtGeFre[g-1,13]+mtGeFre[g-1,14]+x*(1-x)*mtGeFre[g-1,16]+
                                x*(1-x)*y*y*mtGeFre[g-1,17]+x*(1-y)*(1-y)*(1-x)*mtGeFre[g-1,18]+x*(1-y)*(1-y)*(1-x)*mtGeFre[g-1,19]+x*y*y*(1-x)*mtGeFre[g-1,20])
        #aQb/aQB or Aqb/AqB
        mtGeFre[g,15] <- 0.5*(y*(1-y)*mtGeFre[g-1,7]+z*(1-z)*mtGeFre[g-1,9]+y*(1-y)*mtGeFre[g-1,11]+z*(1-z)*mtGeFre[g-1,13]+mtGeFre[g-1,15]+
                                y*(1-x)*(1-x)*(1-y)*mtGeFre[g-1,17]+y*(1-x)*(1-x)*(1-y)*mtGeFre[g-1,18]+x*x*y*(1-y)*mtGeFre[g-1,19]+x*x*y*(1-y)*mtGeFre[g-1,20])
        #aQb/Aqb or aQB/AqB
        mtGeFre[g,16] <- 0.5*(x*x*mtGeFre[g-1,10]+(1-x)*(1-x)*mtGeFre[g-1,16]+ (1-x)*(1-x)*(1-y)*y*mtGeFre[g-1,17]+(1-x)*(1-x)*(1-y)*y*mtGeFre[g-1,18]+
                                x*x*y*(1-y)*mtGeFre[g-1,19]+x*x*(1-y)*y*mtGeFre[g-1,20])
        #aQB/Aqb 
        mtGeFre[g,17] <- 0.5*((1-x)*(1-x)*(1-y)*(1-y)*mtGeFre[g-1,17]+(1-x)*(1-x)*y*y*mtGeFre[g-1,18]+y*y*x*x*mtGeFre[g-1,19]+x*x*(1-y)*(1-y)*mtGeFre[g-1,20])
        #aQb/AqB
        mtGeFre[g,18] <- 0.5*((1-x)*(1-x)*y*y*mtGeFre[g-1,17]+(1-x)*(1-x)*(1-y)*(1-y)*mtGeFre[g-1,18]+(1-y)*(1-y)*x*x*mtGeFre[g-1,19]+x*x*y*y*mtGeFre[g-1,20])
        #aqB/AQb
        mtGeFre[g,19] <- 0.5*(x*x*y*y*mtGeFre[g-1,17]+x*x*(1-y)*(1-y)*mtGeFre[g-1,18]+(1-y)*(1-y)*(1-x)*(1-x)*mtGeFre[g-1,19]+(1-x)*(1-x)*y*y*mtGeFre[g-1,20])
        #aqb/AQB
        mtGeFre[g,20] <- 0.5*(x*x*(1-y)*(1-y)*mtGeFre[g-1,17]+x*x*y*y*mtGeFre[g-1,18]+y*y*(1-x)*(1-x)*mtGeFre[g-1,19]+(1-x)*(1-x)*(1-y)*(1-y)*mtGeFre[g-1,20])
      }
      rval <- mtGeFre[Gn,i]
    }else if(CrosType == "BCP1"){
      mtGeFre <- matrix(0,nrow=100,ncol=8)	 ##mtGeFre: Matrix of Frequency of the Genotypes ,where nrow is Gn() and ncol is 1~8().
      mtGeFre[1,1] <- 1											##mtGeFre[1,1]=1 means the Frequency of F1 which comes frome AAQQBB * aaqqbb -> AaQqBb equals to 1.
      for (g in 2:Gn){
        # AQB/aqb
        mtGeFre[g,1] <- 0.5*(1-x)*(1-y)*mtGeFre[g-1,1]
        # AQB/aqB
        mtGeFre[g,2] <- 0.5*((1-x)*y*mtGeFre[g-1,1]+(1-x)*mtGeFre[g-1,2])
        # AQB/aQb
        mtGeFre[g,3] <- 0.5*(x*y*mtGeFre[g-1,1]+(1-z)*mtGeFre[g-1,3])
        # AQB/aQB
        mtGeFre[g,4] <- 0.5*(x*(1-y)*mtGeFre[g-1,1]+x*mtGeFre[g-1,2]+z*mtGeFre[g-1,3]+mtGeFre[g-1,4])
        # AQB/Aqb
        mtGeFre[g,5] <- 0.5*(x*(1-y)*mtGeFre[g-1,1]+(1-y)*mtGeFre[g-1,5])
        # AQB/AqB
        mtGeFre[g,6] <- 0.5*(x*y*mtGeFre[g-1,1]+x*mtGeFre[g-1,2]+y*mtGeFre[g-1,5]+mtGeFre[g-1,6])
        # AQB/AQb
        mtGeFre[g,7] <- 0.5*((1-x)*y*mtGeFre[g-1,1]+z*mtGeFre[g-1,3]+y*mtGeFre[g-1,5]+mtGeFre[g-1,7])
        # AQB/AQB
        mtGeFre[g,8] <- 0.5*((1-x)*(1-y)*mtGeFre[g-1,1]+(1-x)*mtGeFre[g-1,2]+(1-z)*mtGeFre[g-1,3]+mtGeFre[g-1,4]+(1-y)*mtGeFre[g-1,5]+mtGeFre[g-1,6]+mtGeFre[g-1,7])+mtGeFre[g-1,8]
      }
      rval <- mtGeFre[Gn,i]
    }else if(CrosType == "BCP2"){
      mtGeFre <- matrix(0,nrow=100,ncol=8)					##mtGeFre: Matrix of Frequency of the Genotypes ,where nrow is Gn() and ncol is 1~8().
      mtGeFre[1,8] <- 1											##mtGeFre[1,1]=1 means the Frequency of F1 which comes frome AAQQBB * aaqqbb -> AaQqBb equals to 1.
      for (g in 2:Gn){
        # aqb/aqb
        mtGeFre[g,1] <- 0.5*((1-x)*(1-y)*mtGeFre[g-1,8]+(1-x)*mtGeFre[g-1,7]+(1-z)*mtGeFre[g-1,6]+mtGeFre[g-1,5]+(1-y)*mtGeFre[g-1,4]+mtGeFre[g-1,3]+mtGeFre[g-1,2])+mtGeFre[g-1,1]
        # aqb/aqB
        mtGeFre[g,2] <- 0.5*((1-x)*y*mtGeFre[g-1,8]+z*mtGeFre[g-1,6]+y*mtGeFre[g-1,4]+mtGeFre[g-1,2])
        # aqb/aQb
        mtGeFre[g,3] <- 0.5*(x*y*mtGeFre[g-1,8]+x*mtGeFre[g-1,7]+y*mtGeFre[g-1,4]+mtGeFre[g-1,3])
        # aqb/aQB
        mtGeFre[g,4] <- 0.5*(x*(1-y)*mtGeFre[g-1,8]+(1-y)*mtGeFre[g-1,4])
        # aqb/Aqb
        mtGeFre[g,5] <- 0.5*(x*(1-y)*mtGeFre[g-1,8]+x*mtGeFre[g-1,7]+z*mtGeFre[g-1,6]+mtGeFre[g-1,5])
        # aqb/AqB
        mtGeFre[g,6] <- 0.5*(x*y*mtGeFre[g-1,8]+(1-z)*mtGeFre[g-1,6])
        # aqb/AQb
        mtGeFre[g,7] <- 0.5*((1-x)*y*mtGeFre[g-1,8]+(1-x)*mtGeFre[g-1,7])
        # aqb/AQB
        mtGeFre[g,8] <- 0.5*(1-x)*(1-y)*mtGeFre[g-1,8]
      }
      rval <- mtGeFre[Gn,i]
    }else{
      stop("argument CrosType must be 'Fn' or 'BCP1' or 'BCP2'")
    }
  }else{
    if(CrosType == "Fn"){
      mtGeFre <- matrix(0,nrow=100,ncol=9)	 ##mtGeFre: Matrix of Frequency of the Genotypes ,where nrow is Gn() and ncol is 1~4().
      mtGeFre[1,5] <- 1											##mtGeFre[1,1]=1 means the Frequency of F1 which comes frome AABB * aabb -> AaBb equals to 1.
      for (g in 2:Gn){
        # ab/ab
        mtGeFre[g,1] <- 0.25*(4*mtGeFre[g-1,1]+mtGeFre[g-1,2]+mtGeFre[g-1,4]+(1-x)^2*mtGeFre[g-1,5])
        # aB/ab
        mtGeFre[g,2] <- 0.25*(2*mtGeFre[g-1,2]+x*(1-x)*mtGeFre[g-1,5])
        # aB/aB
        mtGeFre[g,3] <- 0.25*(2*mtGeFre[g-1,2]+4*mtGeFre[g-1,3]+x^2*mtGeFre[g-1,5]+mtGeFre[g-1,6])
        # Ab/ab
        mtGeFre[g,4] <- 0.25*(2*mtGeFre[g-1,4]+x*(1-x)*mtGeFre[g-1,5])
        # AB/ab or Ab/aB
        mtGeFre[g,5] <- 0.25*((1-x)^2+x^2)*mtGeFre[g-1,5]
        # AB/aB
        mtGeFre[g,6] <- 0.25*(2*mtGeFre[g-1,6]+x*(1-x)*mtGeFre[g-1,5])
        # Ab/Ab
        mtGeFre[g,7] <- 0.25*(2*mtGeFre[g-1,8]+4*mtGeFre[g-1,7]+x^2*mtGeFre[g-1,5]+mtGeFre[g-1,4])
        # AB/Ab
        mtGeFre[g,8] <- 0.25*(2*mtGeFre[g-1,8]+x*(1-x)*mtGeFre[g-1,5])
        # ab/ab
        mtGeFre[g,9] <- 0.25*(4*mtGeFre[g-1,9]+mtGeFre[g-1,8]+mtGeFre[g-1,6]+(1-x)^2*mtGeFre[g-1,5])
        rval <- mtGeFre[Gn,i]
      }
    }else if(CrosType == "BCP1"){
      mtGeFre <- matrix(0,nrow=100,ncol=4)	 ##mtGeFre: Matrix of Frequency of the Genotypes ,where nrow is Gn() and ncol is 1~4().
      mtGeFre[1,1] <- 1											##mtGeFre[1,1]=1 means the Frequency of F1 which comes frome AABB * aabb -> AaBb equals to 1.
      for (g in 2:Gn){
        # AB/ab
        mtGeFre[g,1] <- 0.5*(1-x)*mtGeFre[g-1,1]
        # AB/aB
        mtGeFre[g,2] <- 0.5*(x*mtGeFre[g-1,1]+mtGeFre[g-1,2])
        # AB/Ab
        mtGeFre[g,3] <- 0.5*(x*mtGeFre[g-1,1]+mtGeFre[g-1,3])
        # AB/AB
        mtGeFre[g,4] <- 0.5*((1-x)*mtGeFre[g-1,1]+mtGeFre[g-1,2]+mtGeFre[g-1,3])+mtGeFre[g-1,4]
      }
      rval <- mtGeFre[Gn,i]
    }else if(CrosType == "BCP2"){
      mtGeFre <- matrix(0,nrow=100,ncol=4)	 ##mtGeFre: Matrix of Frequency of the Genotypes ,where nrow is Gn() and ncol is 1~4().
      mtGeFre[1,4] <- 1											##mtGeFre[1,1]=1 means the Frequency of F1 which comes frome AABB * aabb -> AaBb equals to 1.
      for (g in 2:Gn){
        # ab/ab
        mtGeFre[g,1] <- 0.5*(mtGeFre[g-1,2]+mtGeFre[g-1,3]+(1-x)*mtGeFre[g-1,4])+mtGeFre[g-1,1]
        # ab/aB
        mtGeFre[g,2] <- 0.5*(mtGeFre[g-1,2]+x*mtGeFre[g-1,4])
        # ab/Ab
        mtGeFre[g,3] <- 0.5*(mtGeFre[g-1,3]+x*mtGeFre[g-1,4])
        # AB/AB
        mtGeFre[g,4] <- 0.5*(1-x)*mtGeFre[g-1,4]
      }
      rval <- mtGeFre[Gn,i]
    }else{
      stop("argument CrosType must be 'Fn' or 'BCP1' or 'BCP2'")
    }
  }
  return(rval)
}
