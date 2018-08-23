\name{PhenoData}
\alias{PhenoData}
\docType{data}
\title{
Data about all phenotype
}
\description{
A data frame with 647 observations on the following 10 variables, the phenotype value is calculated by best linear unbias estimates.
}
\usage{data("PhenoData")}
\format{
  A data frame with 647 observations on the following 13 variables.
  \describe{
    \item{\code{Indi}}{The ID of individuals, which is a factor with levels of Pop1 with 131 individuals, Pop2 with 120 individuals, Pop3 with 200 individuals and Pop4 with 200 individuals}
    \item{\code{Popu}}{a factor with levels \code{Pop1} \code{Pop2} \code{Pop3} \code{Pop4} for multiple-line cross population}
    \item{\code{Ma}}{a factor with levels \code{A} \code{D} for parent 1}
    \item{\code{Pa}}{a factor with levels \code{B} \code{C} \code{E} \code{F} for parent 2}
    \item{\code{newEC1}}{a numeric vector for trait newEC1}
    \item{\code{newEC2}}{a numeric vector for trait newEC2}
    \item{\code{newEC3}}{a numeric vector for trait newEC3}
    \item{\code{BM1}}{a numeric vector for trait BM1}
    \item{\code{BM2}}{a numeric vector for trait BM2}
    \item{\code{BM3}}{a numeric vector for trait BM3}
	\item{\code{predPH1}}{a numeric vector for trait predPH1}
    \item{\code{predPH2}}{a numeric vector for trait predPH2}
    \item{\code{predPH3}}{a numeric vector for trait predPH3}
  }
}
\details{
Popu is the levels of multiple cross populations with different parent crosses, if it is filled with several levels, then linkage mapping or association mapping analysis is done for multiple cross populations, otherwise just for only one population
}
\references{
Hu, H., Meng, Y., Wang, H., Liu, H., and Chen, S. (2012). Identifying quantitative trait loci and determining closely related, stalk traits for rind penetrometer resistance in a high-oil maize, population. Theoretical and Applied Genetics, 124(8), 1439-1447.

Hu, H., Liu, W., Fu, Z., Homann, L., Technow, F., and Wang, H., et al. (2013). Qtl mapping of stalk bending strength in a recombinant inbred line maize population. Theoretical and Applied Genetics, 126(9), 2257-66.
}
\examples{
data(PhenoData)
}
\keyword{datasets}