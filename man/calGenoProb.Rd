\name{calGenoProb}
\alias{calGenoProb}
\title{
Calculate conditional probability of QTL and missing marker genotype
}
\description{
This function is used to calculate conditional probability of QTL and  missing marker genotype based on flanking markers
}
\usage{
calGenoProb(GeneticMap, GenoData, method, croType = NULL, steps = 0, Gn = 2)
}
\arguments{
  \item{GeneticMap}{
	Genetic map data, where the first column is marker names, the second is chromosome numeric ID and the last is genetic distance position
}
  \item{GenoData}{
	A data frame for all maker information, where the NA represents missing markers and 0/1/2 stands for parent genotype AA/AB/BB
}
  \item{method}{
	Method for association mapping(AM) or linkage mapping(LM)
}
  \item{croType}{
	The type of cross populations including Fn, BC1P1, BCP2, F2, DH, and RIL
}
  \item{steps}{
	Maximal genetic distance for computing QTL genotype conditional probabilities, if this argument set to 0, only missing markers genotype are estimated
}
  \item{Gn}{
	The number of generations for Fn, BC1Fn and BC2Fn populations, this parameter is invalid for DH, F2 and RIL populations
}
}
\details{
	This function is to calculate conditional probability for QTL genotype and estimate missing marker genotype based on flanking markers. Here Haldane`s mapping function is used to convert genetic distances in centiMorgans into recombination fractions. Only BCP1, BCP2, DH, Fn and RIL populations are considered in this function, and argument "Gn" is invalid for cross populations F2, DH and RIL. Here, numeric 0 and 2 represent parent polymorphic alleles genotype AA and BB respectively, 1 stands for hybridization genotype AB, and NA represents missing markers.
}
\value{
	Return a data frame about all markers information with some genetic distance, including linkage groups, genetic distances and different individuals genotype.
}
\references{
  Haley, C. S., and Knott, S. A. (1992). A simple regression method for mapping quantitative trait loci in line crosses using flanking markers. Heredity, 69(4), 315.
  APAHaldane, J. B. S. (1918). The combination of linkage values, and the calculation of distance between linked factors. j genet 8:299-309. Journal of Genetics, 8(4), 299-309.
  
  Hospital, F., Dillmann, C., and Melchinger, A. E. (1996). A general algorithm to compute multilocus genotype frequencies under various mating systems. Bioinformatics, 12(6), 455-62.
  
  Martinez, O., and Curnow, R. N. (1994). Missing markers when estimating quantitative trait loci using regression mapping. Heredity, 73(2), 198-206.
  
  Jiang, C., and Zeng, Z. B. (1997). Mapping quantitative trait loci with dominant and missing markers in various crosses from two inbred lines. Genetica, 101(1), 47-58.
}
\author{
	Junhui Li
}
\note{
	This is the first step for linkage mapping with parameter {steps}=0 before permutation test
}
\examples{
data(GeneticMap)
data(GenoData)
GenoData[!GenoData==2 & !GenoData==0] <- NA
steps=0
croType <- "RIL"
Gn=7
method="LM"
#calGenoProb(GeneticMap, GenoData, method, croType, steps, Gn)
}
\keyword{Missing Marker}
