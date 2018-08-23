\name{JM}
\alias{JM}
\title{
Association mapping or linkage mapping analysis in multiple-line cross populations
}
\description{
Joint Mapping for single trait in multiple-line cross populations based on linear regression
}
\usage{
	JM(vecPheno, vecH2, PhenoData, method, vecThrVal, GenoData_EST, GenoData_QTL)
}

\arguments{
  \item{vecPheno}{
  The vector for interest traits name
}
  \item{vecH2}{
  Heritability of interest traits respectively
}
  \item{PhenoData}{
  The data frame for interest trait, and the first fourth column is fixed
}
  \item{method}{
  Method for association mapping(AM) or linkage mapping(LM)
}
  \item{vecThrVal}{
  A vector of threshold value
}
  \item{GenoData_EST}{
  Information of all marker genotype, where missing markers are estimated
}
  \item{GenoData_QTL}{
  Information of all marker genotype with determined distance step, where missing markers are estimated, this dataset is invalid for association mapping(AM)
}
}
\details{
Some required files must be prepared, including GenoData, PhenoData and others
}
\value{
\item{$traitname}{Result of JM for interest trait}
}
\references{
Liu, W., Reif, J. C., Ranc, N., Porta, G. D., and Wurschum, T. (2012). Comparison of biometrical approaches for qtl detection in multiple segregating families. Theoretical and Applied Genetics, 125(5), 987-998.

Steinhoff, J., Liu, W., Maurer, H. P., Wurschum, T., Friedrich, H. L. C., and Ranc, N., et al. (2011). Multiple-line cross quantitative trait locus mapping in european elite maize. Crop Science, 51(6), 2505.
}
\author{
  Junhui Li
}
\examples{
data(GenoData)
data(PhenoData)
data(GenoData_EST)
data(GenoData_S2)
method <- "LM"
vecPheno <- c("newEC1","newEC2")
croType <- "RIL"
GenoData_QTL <- GenoData_S2
#vecThrVal <- c(0.0001,0.0001,0.05,0.05) #for association mapping
vecThrVal <- c(0.0001,0.0001,3.3,3.2)
vecH2 <- c(0.896,0.865)
#QTLresult <- JM(vecPheno, vecH2, PhenoData, method, vecThrVal, GenoData_EST, GenoData_QTL)
}
\keyword{JM}
