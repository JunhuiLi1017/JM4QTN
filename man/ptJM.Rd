\name{ptJM}
\alias{ptJM}
\title{
  Thresholds are determined by Permutation test for QTN association or linkage mapping in multiple-line cross population
}
\description{
 Thresholds are determined by Permutation test for QTN association or linkage mapping in multiple-line cross population. Critical thresholds includes p value for association analysis, while P valude and LOD statistic for linkage mapping analysis in single interest trait
}
\usage{
  ptJM(vecPheno, PhenoData, GenoData_EST, method, npt, alpha,
  selection, tolerance, Trace, select, sle, sls, Choose)
}
\arguments{
  \item{vecPheno}{
  The name vector for interest traits
}
  \item{PhenoData}{
  The data set for trait with best linear unbias estimate(BLUE) or predict(BLUP)
}
  \item{GenoData_EST}{
  Information of chromosome id, genetic distance position and genotype for all marker, where missing markers are estimated
}
  \item{method}{
  Method for association mapping(AM) or linkage mapping(LM)
}
  \item{npt}{
  The Number times of permutation test
}
  \item{alpha}{
  Significant levels, default is 0.1
}
  \item{selection}{
  Model selection method including "forward" and "stepwise",forward selection starts with no effects in the model and adds effects, while stepwise regression is similar to the forward method except that effects already in the model do not necessarily stay there
}
  \item{tolerance}{
  Tolerance value for controlling multicollinearity, default is 1e-7
}
  \item{Trace}{
  Statistic for multiple analysis of variance, including Wilks` lamda, Pillai Trace and Hotelling-Lawley's Trace
}
  \item{select}{
  Specifies the criterion that uses to determine the order in which effects enter and/or leave at each step of the specified selection method including Akaike Information Criterion(AIC), the Corrected form of Akaike Information Criterion(AICc),Bayesian Information Criterion(BIC),Schwarz criterion(SBC),Hannan and Quinn Information Criterion(HQ), Significant Levels(SL)
}
  \item{sle}{
  Specifies the significance level for entry
}
  \item{sls}{
  Specifies the significance level for staying in the model
}
  \item{Choose}{
  Chooses from the list of models at the steps of the selection process the model that yields the best value of the specified criterion. If the optimal value of the specified criterion occurs for models at more than one step, then the model with the smallest number of parameters is chosen. If you do not specify the Choose option, then the model selected is the model at the final step in the selection process
}
}
\details{
Thresholds for Multiple-line cross population QTN association or linkage mapping are determined by Permutation test. Critical thresholds including p value and LOD value are determined to mapping single interest trait. Type of Multiple analysis of variance is invalid for this function
}
\value{
Return a data frame for association mapping(AM) and linkage mapping with 0.05, 0.1 and alpha level
}
\references{
Doerge, R. W., and Churchill, G. A. (1996). Permutation tests for multiple loci affecting a quantitative character. Genetics, 142(1), 285-94.
}
\author{
  Junhui Li
}
\examples{
data(GenoData)
data(PhenoData)
data(GenoData_EST)
vecPheno <- c("newEC1","newEC2")
npt <- 1000
alpha <- 0.15
#vecThrVal <- ptJM(vecPheno, PhenoData, GenoData_EST, method, npt, alpha)
}
\keyword{ Permutation test}