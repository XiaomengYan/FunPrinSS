library(R6)

randFLR <- R6Class(
  "randFLR", 
  portable = FALSE, 
  public = list(
    N = NA, ## All sample size
    C = NA, ## subsample size
    Ns = NA,  ## True dimension of eigen space
    R = NA,   ## The dimension of principal subspace
    EigenValue = numeric(0), ## True eigen values of length Ns
    
    t_seq = numeric(0), ## random function are evaluated at t_seq
    nTXX = numeric(0), ## the length of t_seq
    X = numeric(0), ## generated random function matrix (N x nTXX)
    PHI = numeric(0), ## coefficient function (nTXX x 1)
    Y = numeric(0),  ## Response vector (N x 1)
    hatY = numeric(0), ## prediction 
    hatPHI = numeric(0), ## full sample estimation Phi
    
    X.C = numeric(0),  ## Full sample covariance operator
    X.pca = NA,        ## eigen-decomposition of X.C
    X.lambda = numeric(0), ## eigen-value of X.C
    X.fmat = numeric(0),  ## eigen functions of X.C
    Zscore = numeric(0),
    
    initialize = function(Sample.size, Eigen.value,Phi.coeff,Dim.of.Subspace, Leverage.Type){
      require(RSpectra)
      require(mvtnorm)
      require(fda)
      require(matrixcalc)
      nTXX <<- 2e3
      N <<- Sample.size
      R <<- Dim.of.Subspace
      Leverage.Type <<- Leverage.Type
      EigenValue <<- Eigen.value
      Ns <<- length(EigenValue)
      PhiCoeff <<- Phi.coeff
      GenerateRandFun()
      Cov.all()
    },
    
    UniformSample = function(Sub.Sample.size){
      C <<- Sub.Sample.size
      uniIndex = sample(1:N, size = C, replace = TRUE)
      uniResMat = X[uniIndex,]
      uniResRes = Y[uniIndex,]
      return(list("sX" = uniResMat/sqrt(C*1/N),"sY" = uniResRes/sqrt(C*1/N)))
    }, 
    
    
    WeightedSample = function(Sub.Sample.size, alpha = 1){
      C <<- Sub.Sample.size
      normVec = rowSums(X^2)
      weightVec = normVec/sum(normVec)*alpha + (1-alpha)*1/N
      weightIndex = sample(x = 1:N,size = C, replace = TRUE, prob = weightVec)
      weightResMat = X[weightIndex,] / sqrt( C * weightVec[weightIndex] )
      weightResRes = Y[weightIndex]/ sqrt( C * weightVec[weightIndex] )
      return(list("sX" = weightResMat,"sY" = weightResRes))
    }, 
    
    LeverageSample = function(Sub.Sample.size){
      sXp = self$WeightedSample(Sub.Sample.size = Sub.Sample.size, alpha = 0.5)
      PilotList = self$Compute.subSample.Cov(sXp$sX)
      sX.C = PilotList$sX.C
      sX.pca = PilotList$sX.pca
      sX.lambda = PilotList$sX.lambda
      sX.fmat = PilotList$sX.fmat
      
      psi_np = X %*% sX.fmat/nTXX
      x_residual = X - tcrossprod(psi_np,sX.fmat)
      psi_np = sweep(psi_np, 2, sqrt(sX.lambda), "/")
      Pvec = rowSums(psi_np^2) + rowSums(x_residual^2) / (nTXX * sX.lambda[R])
      
      weightVec = Pvec/sum(Pvec)            ## Normalization
      LeverageIndex = sample(1:N,size = C,replace = TRUE,prob = weightVec)
      LeverageResMat = X[LeverageIndex,] / sqrt( C* weightVec[LeverageIndex] )
      LeverageResRes = Y[LeverageIndex] / sqrt( C* weightVec[LeverageIndex] )
      return(list("sX" = LeverageResMat,"sY" = LeverageResRes))
    },
    
    Compute.subSample.Cov = function(sX_){
      sX.TensorMat = crossprod(sX_)   
      sX.C = sX.TensorMat/N
      if(Leverage.Type != "VN"){
        sX.pca =  eigs_sym(sX.C, k = R, which = "LM")
        sX.lambda = sX.pca$values/nTXX
        sX.fmat = sX.pca$vectors*sqrt(nTXX)
      }else{
        sX.pca = eigen(sX.C,symmetric = TRUE)
        sX.lambda = sX.pca$values[1:R]/nTXX
        sX.fmat = sX.pca$vectors[,1:R]*sqrt(nTXX)
      }
      return(list("sX.C" = sX.C, "sX.pca" = sX.pca, "sX.lambda" = sX.lambda, "sX.fmat" = sX.fmat))
    },
    
    LSestimator = function(sX_,sY_){
      hatZ = 1/N*crossprod(sX_,sY_)
      sX.TensorMat = crossprod(sX_)
      sX.C = sX.TensorMat/N
      if(Leverage.Type!="VN"){
        sX.pca =  eigs_sym(sX.C, k = R, which = "LM")
        sX.lambda = sX.pca$values/nTXX
        sX.fmat = sX.pca$vectors*sqrt(nTXX)
      }else{
        sX.pca = eigen(sX.C,symmetric = TRUE)
        sX.lambda = sX.pca$values[1:R]/nTXX
        sX.fmat = sX.pca$vectors[,1:R]*sqrt(nTXX)
      }
      invsX.C = tcrossprod(sweep(sX.fmat, MARGIN = 2, sX.lambda, "/"),sX.fmat)
      tildePHI = invsX.C%*%hatZ / nTXX
      tildeY = X%*%tildePHI/nTXX
      return(list("tildePHI"=tildePHI,"tildeY" = tildeY))
    },
    
    ## Two comparison metrics
    Pred.Error = function(LSlist){
      tildeY = LSlist[[2]]
      return(sum((hatY - tildeY)^2)/N)
    }, 
    
    Est.Error = function(LSlist){
      tildePHI = LSlist[[1]]
      err = norm(hatPHI-tildePHI, type = "2")/sqrt(nTXX)
      return(err)
    }
    
    
  ), ## public list
  private = list(
    GenerateRandFun = function(){
      t_seq <<- seq(from = 0, to  = 1, length.out = nTXX)
      fmat = fourier(t_seq,nbasis = Ns+1)  # eigenfunction matrix
      fmat = fmat[,-1]
      Phi = sweep(fmat,MARGIN = 2,PhiCoeff,FUN = "*")
      Phi = rowSums(Phi)
      
      if(Leverage.Type == "NU"){
        Z = matrix(rnorm(Ns*N),N, Ns)
      }
      if(Leverage.Type == "MN"){
        Z = matrix(rt(n = Ns*N,df = 3),N,Ns)
      }
      if(Leverage.Type == "VN"){
        Z = matrix(rt(n = Ns*N,df = 1),N,Ns)
      }
      Z = sweep(Z,2, sqrt(EigenValue), "*")
      X <<- Z%*%t(fmat)
      PHI <<- Phi
      Y <<- Z%*%PhiCoeff + rnorm(N)
    }, 
    
    Cov.all = function(){
      X.TensorMat = crossprod(X)
      X.C <<- X.TensorMat/N
      
      if(Leverage.Type!="VN"){
        X.pca <<- eigs_sym(A = X.C, k = R, which = "LM")
        X.lambda <<- X.pca$values/nTXX
        X.fmat <<- X.pca$vectors*sqrt(nTXX)
      }else{
        X.pca <<- eigen(X.C,symmetric = TRUE)
        X.lambda <<- X.pca$values[1:R]/nTXX
        X.fmat <<- X.pca$vectors[,1:R]*sqrt(nTXX)
      }
      
      hatZ = crossprod(X, Y)/N
      invX.C = tcrossprod(sweep(X.fmat,MARGIN = 2, X.lambda, "/"),X.fmat)
      hatPHI <<- invX.C%*%hatZ/nTXX
      hatY <<- X%*%hatPHI/nTXX
    }
  )
)