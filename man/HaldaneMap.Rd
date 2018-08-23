\name{HaldaneMap}
\alias{HaldaneMap}
\title{
Haldane Map function
}
\description{
It can convert genetic distance in centimorgan into recombination frequency
}
\usage{
HaldaneMap(x)
}
\arguments{
  \item{x}{
	The genetic distance in centimorgan
}
}

\value{
Return recombination frequency of two markers
}
\references{
Haldane, J. B. S. (1918). The combination of linkage values, and the calculation of distance between linked factors. j genet 8:299-309. Journal of Genetics, 8(4), 299-309.
}
\author{
JunhuiLi
}
\note{
arguments x should be more than 0
}

\examples{
x <- 2
HaldaneMap(x)
}
\keyword{Missing Marker}
