\name{Eagd}
\alias{Eagd}
\title{
	Expectation additive genetic deviation
}
\description{
	Compute expectation for additive genetic deviation based on flanking markers for different cross population
}
\usage{
	Eagd(marType, croType, Gn = 2, x, y = 0)
}
\arguments{
  \item{marType}{
	Flanking markers genotype information, 2/1/0/N represents AA/Aa/aa/missing markers
}
  \item{croType}{
	Cross population type, including Fn, F2, DH, RIL and BCP1, BCP2
}
  \item{Gn}{
	The number of generation
}
  \item{x}{
	Genetic distance from left nearest flanking marker
}
  \item{y}{
	Genetic distance from right nearest flanking marker
}
}
\details{
	If argument y is set up to 0, then frequency of genotype is computed only from left or right nearest flanking marker. For more information, you can see Hospital reference and JunhuiLi reference supplements.

	This function is called by functions "calGenoProb".
}
\value{
	Return expectation for additive genetic deviation based on flanking markers
}
\references{
	Hospital, F., Dillmann, C., & Melchinger, A. E. (1996). A general algorithm to compute multilocus genotype frequencies under various mating systems. Bioinformatics, 12(6), 455-462.

	Junhui Li, Haixiao Hu, Yujie Meng, Kun Cheng, Guoliang Li, Wenxin Liu, and Shaojiang Chen.(2016)Pleiotropic QTL detection for stalk traits in maize and related R package programming. Journal of China Agricultural University. DOI 10.11841/j.issn.1007-4333.2016.06.00(in chinese)
}
\author{
	JunhuiLi
}
\examples{
	Eagd("22", "RIL", Gn = 2, x=1.05, y=0.9) 
}
\keyword{ Eagd }
