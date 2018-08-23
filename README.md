# JM4QTN
## 1. Introduction

* What is the package JM4QTN?

Joint Mapping for QTN, it combines composite interval mapping(CIM) and regression analysis to do association mapping and linkage mapping and identify pleiotropic or linked QTL in multiple-traits both in single or multiple line cross populations. It can also compute conditional probability of missing marker genotype using flanking markers. Thresholds are determined by permutation test and parametric bootstrap test. Stepwise regression with multivariate and univariate is carried out with several information criteria to select cofactors across families.

## 2. Statistical and coding details in this package

* What is the main features in SteReg?

** Missing marker genotype estimation: Hidden Markov model based on genotype frequency for different mapping population(RIL, F2, Fn, BCP1, BCP2, DH)

** Linear model selection (stepwise regression) programmed with Cpp: (1)Variable enter direction: forward stepwise selection and hybrid approaches (2)Information criteria: AIC/AICc/BIC/SBC/HQ and Significant Level (3)Multicollinearity: tolerance and VIF (4)Independent variables: continuous effects and continuous-nesting-class effects (5)Dependent variables: univariate and multivariate

** P/LOD value calculation(likelihood ratio test): (1)F test for univariate based on full and reduced linear models (2)Approximate F test for multivariate based on Wilks' lambda, Pillai’s' Trace and Hotelling-Lawley Trace with full and reduced linear models

** Resample:(1)Permutation test (2)bootstrap and (3)cross-validation

## 3. Usage and Examples
	#install.package("JM4QTN")
	
	library(JM4QTN)
	
	#load data
	data(GeneticMap)
	data(GenoData)
	data(PhenoData)
	data(GenoData_EST)
	data(GenoData_S2)
	
	#Estimate missing markers: step=0 and calculate conditional probability

	cross <- "RIL"
	method <- "LM"    #"AM"
	steps <- 0
	vecPheno <- c("newEC1","newEC2")
	Data_calGenoProb <- calGenoProb(GeneticMap,MarkerData,method,cross,0)
	str(Data_calGenoProb);dim(Data_calGenoProb);Data_calGenoProb[1:10,1:10]

	# Calculate QTL conditional probability with max step=2
	
	steps <- 2
	RMdata_2 <- calGenoProb(GeneticMap,MarkerData,method,cross,steps)
	RMdata_2[1:10,1:20];dim(RMdata_2)

	# QTL scan
	
	vecPheno <- c("newEC1","newEC2")
	#vecPheno <- c("newEC1")
	croType <- "RIL"
	GenoData_QTL <- GenoData_S2
	
	#permutation test
	vecThrVal <- ptJM1(vecPheno,PhenoData,GenoData_EST,method,npt=100,alpha=0.1,selection="stepwise",tolerance=1e-7,Trace="Pillai",select="SBC",sle=0.05,sls=0.05,Choose="SBC")
	#vecThrVal <- c(0.0001,0.0001,3.3,3.2)
	vecH2 <- c(0.8,0.8)
	JM(vecPheno,vecH2,PhenoData,method,vecThrVal,GenoData_EST,GenoData_QTL)
	
	# Pleiotropic QTL detection
	
	vecThrVal <- ptLP1(vecPheno,PhenoData,GenoData_EST,npt=100,alpha=0.1,selection="stepwise",tolerance=1e-7,Trace="Pillai",select="SBC",sle=0.05,sls=0.05,Choose="SBC")
	CChr <- 2
	#vecThrVal <- c(0.000206,4.62,4.70,6.64)
	Interval <- c(56.896,57.237,57.237,61.323)
	LP(vecPheno,PhenoData,GenoData_EST,GenoData_QTL,vecThrVal,CChr,Interval,nPB=200,alpha=0.05)
