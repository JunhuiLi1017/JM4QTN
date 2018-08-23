\name{LP}
\alias{LP}
\title{
Identification of pleiotropic and linked QTL in multiple line cross populations
}
\description{
Identification of pleiotropic and linked QTL mapping in multiple-line cross populations, which is based on linear regression and likelihood ratio test
}
\usage{
LP(vecPheno, PhenoData, GenoData_EST, GenoData_QTL, vecThrVal, CChr, Interval, nPB, alpha)
}

\arguments{
  \item{vecPheno}{
  The vector for interest traits name
}
  \item{PhenoData}{
  The data frame for interest trait
}
  \item{GenoData_EST}{
  Information about all marker genotype, where missing markers are estimated
}
  \item{GenoData_QTL}{
  Information about all marker genotype with determined distance step, where missing markers are estimated
}
  \item{vecThrVal}{
  A vector of threshold value
}
  \item{CChr}{
  Current chromosome for pleiotropic and linked QTL test
}
  \item{Interval}{
  Support interval for pleiotropic and linked QTL test
}
  \item{nPB}{
  The number times of parametric bootstrap, default is 2000
}
  \item{alpha}{
  Significant levels, default is 0.05
}
}
\details{
Some required files must be prepared, including GeneticMap, GenoData, PhenoData and others
}
\value{
Print information "The interest QTLs are linked QTL" or "The interest QTLs are pleiotropic QTL".
Return a data frame for single trait QTL and pleiotropic or linked QTL information.
}
\references{
Liu, W., Reif, J. C., Ranc, N., Porta, G. D., and Wurschum, T. (2012). Comparison of biometrical approaches for qtl detection in multiple segregating families. Theoretical and Applied Genetics, 125(5), 987-998.

Steinhoff, J., Liu, W., Maurer, H. P., Wurschum, T., Friedrich, H. L. C., and Ranc, N., et al. (2011). Multiple-line cross quantitative trait locus mapping in european elite maize. Crop Science, 51(6), 2505.

Knott, S. A., and Haley, C. S. (2000). Multitrait least squares for quantitative trait loci detection. Genetics, 156(2), 899-911.
}
\author{
Junhui Li
}
\examples{
data(GenoData)
data(PhenoData)
data(GenoData_EST)
data(GenoData_S2)
vecPheno <- c("newEC1","newEC2")
croType <- "RIL"
Gn <- 7
GenoData_QTL <- GenoData_S2
alpha <- 0.1
CChr <- 14
vecThrVal <- c(5.210408e-05,4.99,5.05,6.91)
nPB <- 2000
Interval <- c(70.273,70.999,71.21,73.506)
#QTLresult <- LP(vecPheno,PhenoData,GenoData_EST,GenoData_QTL,vecThrVal,CChr,Interval,nPB,alpha)
}
\keyword{LP}