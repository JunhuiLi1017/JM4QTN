\name{FreGeno}
\alias{FreGeno}
\title{
	Compute frequency of genotype
}
\description{
	Compute frequency of genotype under selfing and back cross within generation n based on flanking markers
}
\usage{
	FreGeno(CrosType, Gn = 2, i, x, y = 0)
}
\arguments{
  \item{CrosType}{
	Cross type, including Fn, BCP1 and BCP2
}
  \item{Gn}{
	The number of generation
}
  \item{i}{
	Numeric value for genotype
}
  \item{x}{
	Genetic distance from left nearest flanking marker
}
  \item{y}{
	Genetic distance from right nearest flanking marker
}
}
\details{
	If argument y is set up to 0, then frequency of genotype is computed only from left or right nearest marker. For more information, you can see Hospital reference and JunhuiLi reference supplements.
	
	This function is called by functions "calGenoProb".
}
\value{
	Frequency of genotype for selfing or back cross population for generation n
}
\references{
Hospital, F., Dillmann, C., & Melchinger, A. E. (1996). A general algorithm to compute multilocus genotype frequencies under various mating systems. Bioinformatics, 12(6), 455-462.

Junhui Li, Haixiao Hu, Yujie Meng, Kun Cheng, Guoliang Li, Wenxin Liu, and Shaojiang Chen.(2016)Pleiotropic QTL detection for stalk traits in maize and related R package programming. Journal of China Agricultural University. DOI 10.11841/j.issn.1007-4333.2016.06.00(in chinese)
}
\author{
JunhuiLi
}
\examples{
FreGeno("Fn",5, 1, 1, 1)
}
\keyword{Missing Marker}