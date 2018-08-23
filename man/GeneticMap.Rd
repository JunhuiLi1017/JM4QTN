\name{GeneticMap}
\alias{GeneticMap}
\docType{data}
\title{
Data about genetic map
}
\description{
A data frame with 1722 observations on the following 3 variables, where the first column is marker name, the second column is chromosome numeric ID and the last column is marker genetic distance position
}
\usage{data("GeneticMap")}
\format{
  A data frame with 1722 observations on the following 3 variables.
  \describe{
    \item{\code{marker}}{a vector for marker name}
    \item{\code{chr}}{a numeric vector for chromosome id}
    \item{\code{pos}}{a numeric vector for marker genetic distance position}
  }
}
\details{
The name of marker must be the same with marker`s name of data in GenoData, and the chromosome ID and genetic distance position are numeric, otherwise an error will occur.
}
\references{
Hu, H., Meng, Y., Wang, H., Liu, H., and Chen, S. (2012). Identifying quantitative trait loci and determining closely related, stalk traits for rind penetrometer resistance in a high-oil maize, population. Theoretical and Applied Genetics, 124(8), 1439-1447.

Hu, H., Liu, W., Fu, Z., Homann, L., Technow, F., and Wang, H., et al. (2013). Qtl mapping of stalk bending strength in a recombinant inbred line maize population. Theoretical and Applied Genetics, 126(9), 2257-66.
}
\examples{
data(GeneticMap)
}
\keyword{datasets}
