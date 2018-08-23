\name{StatPheno}
\alias{StatPheno}
\title{
Statistic analysis for phenotype data
}
\description{
Statistic analysis for phenotype data, including normality test, analysis of variance and best linear unbiased estimation
}
\usage{
StatPheno(phenoData, defineForm = NULL, effNotation = "G")
}
\arguments{
  \item{phenoData}{
	Phenotype data
}
  \item{defineForm}{
	Formula in linear model defined by user
}
  \item{effNotation}{
	Effect nanotation for best linear unbiased estimation
}
}
\details{
Environment effect, block effect, replication and genotype effect are in order for the fouth column of phenoData
}
\value{
	\item{normality_test }{The result of shapiro.test}
	\item{formula }{The formula of the linear model}
	\item{ANOVA }{The result of ANOVA}
	\item{lsmeans }{The least squart means for best linear unbiased estimation}
}
\author{
Junhui Li
}
\examples{
phenoData <- matrix(NA,180,6)
phenoData <- as.data.frame(phenoData)
colnames(phenoData) <- c("E","B","R","G","Height","Weight")
phenoData[,1] <- as.factor(c(rep(c("A"),60),rep(c("B"),60),rep(c("C"),60)))
phenoData[,2] <- as.factor(c(rep(1,30),rep(2,30)))
phenoData[,3] <- 1
phenoData[,4] <- as.factor(c(1:30))
phenoData[,5] <- rnorm(180,175,5)
phenoData[,6] <- rnorm(180,75,10)
StatPheno(phenoData)
}
\keyword{ Phenotype }

