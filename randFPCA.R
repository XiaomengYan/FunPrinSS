library(R6)
randFPCA <- R6Class(
  "randFPCA", 
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
    
    X.C = numeric(0),  ## Full sample covariance operator
    X.pca = NA,        ## eigen-decomposition of X.C
    X.lambda = numeric(0), ## eigen-value of X.C
    X.fmat = numeric(0),  ## eigen functions of X.C
    
    initialize = function(Sample.size, Eigen.value, Dim.of.Subspace, Leverage.Type){
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
      GenerateRandFun()
      Cov.all()
    },
    
    
    UniformSample = function(Sub.Sample.size){
      C <<- Sub.Sample.size
      uniIndex = sample(1:N, size = C, replace = TRUE)
      uniResMat = X[uniIndex,]
      return(uniResMat/sqrt(C*1/N))
    }, 
    
    
    WeightedSample = function(Sub.Sample.size, alpha = 1){
      C <<- Sub.Sample.size
      normVec = rowSums(X^2)
      weightVec = normVec/sum(normVec)*alpha + (1-alpha)*1/N
      weightIndex = sample(x = 1:N,size = C, replace = TRUE, prob = weightVec)
      weightResMat = X[weightIndex,] / sqrt( C * weightVec[weightIndex] )
      return(weightResMat)
    }, 
    
    LeverageSample = function(Sub.Sample.size){
      sXp = self$WeightedSample(Sub.Sample.size = Sub.Sample.size, alpha = 0.5)
      PilotList = self$Compute.subSample.Cov(sXp)
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
      return(LeverageResMat)
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
    
    ## Several Comparision metrics \|\hC_xx - \tC_xx\|_{HS}
    HS.loss = function(Covlist){
      sX.C = Covlist$sX.C
      return(hilbert.schmidt.norm(X.C-sX.C)/nTXX)
    },
    
    # \|\hC_xx - \tC_xx\|
    Operator.loss = function(Covlist){
      sX.C = Covlist$sX.C
      if(Leverage.Type!= "VN"){
      lambdaSeq_max = eigs_sym(A = X.C-sX.C,k = 1,opts = list(retvec = FALSE))$values
      lambdaSeq_min = eigs_sym(A = X.C-sX.C,k = 1,sigma = 0,opts = list(retvec = FALSE))$values
       return(max(abs(c(lambdaSeq_max,lambdaSeq_min)))/nTXX)
      }else{
        lambdaSeq = eigen(X.C-sX.C,symmetric = TRUE,only.values = TRUE)[[1]]
        return(max(abs(lambdaSeq))/nTXX)
      }
    }, 
    
    HS.loss.subspace = function(Covlist){
      sX.fmat = Covlist$sX.fmat
      ResMat =  tcrossprod(X.fmat) - tcrossprod(sX.fmat)
      return(hilbert.schmidt.norm(ResMat)/nTXX)
    },
    
    # \|\hP_R - \tP_R\|
    Operator.loss.subspace = function(Covlist){
        sX.fmat = Covlist$sX.fmat
        ResMat =  tcrossprod(X.fmat) - tcrossprod(sX.fmat)
        if(Leverage.Type != "VN"){
          lambdaSeq_max = eigs_sym(A = ResMat,k = 1,opts = list(retvec = FALSE))$values
          lambdaSeq_min = eigs_sym(A = ResMat,k = 1,sigma = 0,opts = list(retvec = FALSE))$values
          return(max(abs(c(lambdaSeq_max,lambdaSeq_min)))/nTXX)
        }else{
          lambdaSeq = eigen(ResMat,symmetric = TRUE,only.values = TRUE)[[1]]
          return(max(abs(lambdaSeq))/nTXX)
        }
    }, 
    
    # \|\htheta_r - \ttheta_r\|
    Individual.eigen.loss = function(Covlist){
      sX.fmat = Covlist$sX.fmat
      resVec = rep(0,R)
      for (ns in 1:R) {
        resVec[ns] = norm(X.fmat[,ns]-sign(as.numeric(t(X.fmat[,ns])%*%sX.fmat[,ns]))*sX.fmat[,ns],type = "2")
      }
      return(resVec/sqrt(nTXX))
    }
    
  ),
  
  private = list(
    
    ## Function to generate random functions
    GenerateRandFun = function(){
      t_seq <<- seq(from = 0, to = 1, length.out = nTXX)
      fmat =  fourier(t_seq, nbasis =  Ns + 1)
      fmat = fmat[,-1]
      if(Leverage.Type == "NU"){
        Z = matrix(rnorm(Ns*N),N, Ns)
      }
      if(Leverage.Type == "MN"){
        Z = matrix(rt(n = Ns*N, df = 3), N, Ns)
      }
      if(Leverage.Type == "VN"){
        Z = matrix(rt(n = Ns*N, df = 1), N, Ns)
      }
      Z = sweep(Z,2, sqrt(EigenValue), "*")
      X <<- Z%*%t(fmat)
    },
    
    ## Function to calculate full sample covariance operator X.C and eigen-decomposition of X.C
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
    }
    
  )
)