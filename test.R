library(fda)
library(doParallel)
library(foreach)
library(parallel)
library(matrixcalc)
library(RSpectra)
source(file = "randFPCA.R") # Source the randomized FPCA class

N = 1e4 # full sample size
C = 1e2 # subsample size
R = 5 # cut-off dimension
lambdaSeq = 2^(50:1) # eigenvalue sequence
FullSample <- randFPCA$new(Sample.size = N,Eigen.value = lambdaSeq,Dim.of.Subspace = R,Leverage.Type = "NU") # generate full sample (nearly uniform type)

# FunPrinSS
sX = FullSample$LeverageSample(Sub.Sample.size = C) # subsampling according to FunPrinSS
LeverageList = FullSample$Compute.subSample.Cov(sX = sX)
Lhs <- FullSample$HS.loss(LeverageList) # HS-norm loss of Covariance Operator
Lop <- FullSample$Operator.loss(LeverageList) # Operator-norm loss of Covariance Operator
Lsubhs <- FullSample$HS.loss.subspace(LeverageList) # HS-norm loss of Projection Operator
Lsubop <-FullSample$Operator.loss.subspace(LeverageList)# Operator-norm loss of Projection Operator
Lpind <- FullSample$Individual.eigen.loss(LeverageList) # Individual eigenfunction L2 loss
