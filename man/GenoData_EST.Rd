\name{GenoData_EST}
\alias{GenoData_EST}
\docType{data}
\title{
Data about information of molecular markers chromosome id, genetic distance position and genotype with estimated missing markers
}
\description{
A data frame with 647 individual observations, chromosome id and genetic distance position on the following 1722 molecular maker variables
}
\usage{data("GenoData_EST")}
\format{
  A data frame with 649 observations on the following 1722 variables
  \describe{
    A data frame with 649 observations on the following 1722 maker variables, the rows 'chr' and 'pos' is the chromosome id and genetic distance position of markers, other rows are the marker genotype on 647 individuals, where 0 and 2 stands for genotype AA and BB respectivly, decimal is missing marker conditional probability genotype.
  }
}
\details{
This data frame is combined from geneticMap and genotype data with estimated missing markers, which can be calculated by function calGenoProb(GeneticMap, GenoData, steps=0, croType, Gn)
}
\references{
Hu, H., Meng, Y., Wang, H., Liu, H., and Chen, S. (2012). Identifying quantitative trait loci and determining closely related, stalk traits for rind penetrometer resistance in a high-oil maize, population. Theoretical and Applied Genetics, 124(8), 1439-1447.

Hu, H., Liu, W., Fu, Z., Homann, L., Technow, F., and Wang, H., et al. (2013). Qtl mapping of stalk bending strength in a recombinant inbred line maize population. Theoretical and Applied Genetics, 126(9), 2257-66.
}
\examples{
data(GenoData_EST)
}
\keyword{datasets}
