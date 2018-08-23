\name{ptLP}
\alias{ptLP}
\title{
Thresholds are determined by permutation test for joint mapping in multiple-line cross population with multiple trait
}
\description{
  Critical thresholds are determined by permutation test including p value of F test and LOD value with multivariate stepwise regression
}
\usage{
	ptLP(vecPheno, PhenoData, GenoData_EST, npt, alpha,
	selection, tolerance, Trace, select, sle, sls, Choose)
}
\arguments{
  \item{vecPheno}{
  The vector for interest traits name
}
  \item{PhenoData}{
  The data frame for interest trait
}
  \item{GenoData_EST}{
  Information about all markers, where missing markers are estimated
}
  \item{npt}{
  The number times of permutation test
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
  Specifies the criterion that uses to determine the order in which effects enter and/or leave at each step of the specified selection method including BIC BICc SBC AIC AICc and SL(P value)
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
This permutation test is running with multivariate multiple regression, therefore multiple interest traits are need to do this analysis
}
\value{
Return a data frame for critical threshold with 0.05, 0.1 and alpha level
}
\references{
Doerge, R. W., and Churchill, G. A. (1996). Permutation tests for multiple loci affecting a quantitative character. Genetics, 142(1), 285-94.

Knott, S. A., and Haley, C. S. (2000). Multitrait least squares for quantitative trait loci detection. Genetics, 156(2), 899-911.
}
\author{
Junhui Li
}
\examples{
data(GenoData)
data(PhenoData)
data(GenoData_EST)
vecPheno <- c("Mmax","InD")
npt <- 1000
alpha <- 0.15
#vecThrVal <- ptLP(vecPheno, PhenoData, GenoData_EST, npt, alpha)
}
\keyword{ Permutation test}