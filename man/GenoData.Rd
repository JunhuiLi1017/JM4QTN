\name{GenoData}
\alias{GenoData}
\docType{data}
\title{
Data about molecular markers genotype
}
\description{
A data frame with 647 individual observations on the following 1722 variables, where 2, 1 and 0 stands for different markers genotype AA AB and BB, NA for missing marker genotype
}
\usage{data("GenoData")}
\format{
  A data frame with 647 individual observations on the following 1722 marker genotype variables
}
\details{
Missing markers are subsitituted only by NA, otherwise an error will occur.
}
\references{
Hu, H., Meng, Y., Wang, H., Liu, H., and Chen, S. (2012). Identifying quantitative trait loci and determining closely related, stalk traits for rind penetrometer resistance in a high-oil maize, population. Theoretical and Applied Genetics, 124(8), 1439-1447.

Hu, H., Liu, W., Fu, Z., Homann, L., Technow, F., and Wang, H., et al. (2013). Qtl mapping of stalk bending strength in a recombinant inbred line maize population. Theoretical and Applied Genetics, 126(9), 2257-66.
}
\examples{
data(GenoData)
}
\keyword{datasets}
