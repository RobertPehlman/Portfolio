### Functions

projGradient <- function(parameters){
  gammaNotOrtho = array(parameters[1:(p*L)], c(p, L))
  gammaSVD = svd(gammaNotOrtho)
  gammaOrtho = gammaSVD$u %*% t(gammaSVD$v)
  return(c(as.vector(gammaOrtho), parameters[(p*L + 1):(p*(L+4))], abs(parameters[(p*(L+4) + 1):((L+5)*p + 3)])))
} 

#' Orders the parameters from smallest to largest variance
#' 
#' @param parameters parameters to be optimized / sorted
#' @export
orderParams <- function(parameters){
  gammaalphalambda = array(parameters[1:(p*(L+5))], c(p,(L+5)))
  lambda = parameters[((L+4)*p + 1):((L+5)*p)]
  
  sortedpars <- gammaalphalambda[order(lambda),]
  return(c(as.vector(sortedpars), parameters[-(1:((L+5)*p))]))
}


#' Calculate Log-Likelihood of gaussian process.  This one takes in some unused params
#' 
#' @param parameters parameters that are optimized over
#' @param W sparse array with observed data over time points
#' @param a treatment variable
#' @param u treatment time variable
#' @param p number of eigenfunctions to select
#' @param L number of basis functions
#' @export
logLikelihood <- function(parameters, W, a = NULL, u = NULL, p, L, trtFlag = FALSE){
  n = dim(W)[2]
  
  likeGamma0 = array(parameters[1:(p*L)], c(p, L))
  likeAlpha = array(parameters[(p*L + 1):((L+4)* p)], c(p, 4))
  likeLambda0 = parameters[((L+4)* p + 1): ((L+5)*p)]
  likeNu = parameters[((L+5)*p + 1)]
  likeTau = parameters[((L+5)*p + 2)]
  likeH = parameters[((L+5)*p + 3)]
  # Empty array holds all m phi0 functions at nMax2 times and n observations 
  phi0 = array(rep(NA, nMax2 * p * n), c(nMax2, p, n))
  phi1 = array(rep(NA, nMax2 * p * n), c(nMax2, p, n))
  
  rank = rep(0, n)
  
  likely = 0
  
  for (i in 1:n){
    
    times = sum(!is.na(W[,i]))       # number of observed time points
    indices = (!is.na(W[,i])) * 1:nMax2  # Exact value time points on grid
    indices = indices[indices != 0]
    
    if(trtFlag){
      rank[i] = sum(indices < u[i]) + 1    #NEW DEFINITION FOR RANK, INCLUDES TIME U, CHECK ERROR PROPAGATION
    }
    
    B_T=t(OBASIS[indices,])
    phi0temp = t(likeGamma0 %*% B_T)
    
    
    phi0[indices, , i] = phi0temp
    
    if(trtFlag){
      phi1test =  phi0[, , i] * array(rep(pexp(1:nMax2 - u[i], rate = likeH), p), c(nMax2,p))
      
      
      phi1[, , i] = phi1test
      
      
      Alpha01 = array(c(likeAlpha[,1],a[i]*likeAlpha[,2]), c(p,2))
      Alpha23 = array(c(likeAlpha[,3],a[i]*likeAlpha[,4]), c(p,2))
      
      
      PhiStar01 = phi1[, , i] * t(array(rep(Alpha01 %*% c(1,1), nMax2), c(p,nMax2)))
      PhiStar23 = phi1[, , i] * t(array(rep(Alpha23 %*% c(1,1), nMax2), c(p,nMax2)))
      
      B = cbind(phi0[indices, , i] + PhiStar23[indices,], phi1[indices, , i], diag(times))
      B1 = phi0[indices, , i] + PhiStar23[indices,]
      B2 = phi1[indices, , i]
      B3 = diag(times)
      
      Omega = B1 %*% diag(likeLambda0^2) %*% t(B1) + likeNu^2 * B2 %*% diag(p) %*% t(B2) + likeTau^2 * diag(times)
      
      
      
      
      likely = likely + dmvnorm(W[indices,i], mean = (PhiStar01[indices,] %*% rep(1,p)), sigma = Omega, log = TRUE)
    }
    
    
    B1 = phi0[indices, , i]
    
    Omega = B1 %*% diag(likeLambda0^2) %*% t(B1) + likeTau^2 * diag(times)
    
    
    
    
    likely = likely + dmvnorm(W[indices,i], mean = rep(0,times), sigma = Omega, log = TRUE)
    
    
  }
  return(likely)
  #return(c(PhiStar01, PhiStar23))
}

#trueParameters = c(as.vector(Gamma0), rep(0,(p*4)), lambda0, tau, 1, 1)
#logLikelihood(trueParameters, train$w, a = NULL, u = NULL, p, L, trtFlag = FALSE)


#' Calculate Log-Likelihood of gaussian process.  This one only takes in used params
#' 
#' @param parameters parameters that are optimized over
#' @export
simpleLogLikelihood2 <- function(parameters){
  W = in.w
  n = dim(W)[2]
  likeGamma0 = t(tx(parameters[1:(p*L)]))
  
  
  likeLambda0 = parameters[((L)* p + 1): ((L+1)*p)]
  likeTau = parameters[((L+1)*p + 1)]
  
  # Empty array holds all m phi0 functions at nMax2 times and n observations 
  phi0 = array(rep(NA, nMax2 * p * n), c(nMax2, p, n))
  
  likely = 0
  
  for (i in 1:n){
    
    times = sum(!is.na(W[,i]))       # number of observed time points
    indices = (!is.na(W[,i])) * 1:nMax2  # Exact value time points on grid
    indices = indices[indices != 0]
    
    B_T=t(OBASIS[indices,])
    phi0temp = t(likeGamma0 %*% B_T)
    
    phi0[indices, , i] = phi0temp
    
    B1 = phi0[indices, , i]
    
    Omega = B1 %*% diag(likeLambda0^2) %*% t(B1) + likeTau^2 * diag(times)
    
    likely = likely + dmvnorm(W[indices,i], mean = rep(0,times), sigma = Omega, log = TRUE)
  }
  return(-likely)
  
}

#' Calculate gradient of gaussian process.  This one only takes in used params
#' 
#' @param parameters parameters that are optimized over
#' @export
simpleGradient2ndTry2 <- function(parameters){
  W = in.w
  n = dim(W)[2]
  # Gamma = array(parameters[1:(p*L)], c(L, p))
  # Gamma = t(Gamma)
  # 
  # gradGamma0 = array(Gamma[1:(p*L)], c(p, L))
  
  gradGamma0 = t(tx(parameters[1:(p*L)]))
  
  gradLambda0 = parameters[(L* p + 1): ((L+1)*p)]
  gradTau = parameters[((L+1)*p + 1)]
  
  # Empty array holds all m phi0 functions at nMax2 times and n observations 
  phi0 = array(rep(NA, nMax2 * p * n), c(nMax2, p, n))
  gammaGradient = array(rep(0, (p*L)), c(p, L))
  lambdaGradient = numeric(p)
  tauGradient = 0
  
  
  for (i in 1:n){
    times = sum(!is.na(W[,i]))       # number of observed time points
    indices = (!is.na(W[,i])) * 1:nMax2  # Exact value time points on grid
    indices = indices[indices != 0]
    
    O_T=t(OBASIS[indices,])
    phi0temp = t(gradGamma0 %*% O_T)
    
    phi0[indices, , i] = phi0temp
    
    B1 = phi0[indices, , i]
    
    
    Omega = B1 %*% diag(gradLambda0^2) %*% t(B1) + gradTau^2 * diag(times)
    
    #print(Omega)
    
    invOmega = solve(Omega)
    
    dOmega = .5 * (invOmega - invOmega %*% (W[indices,i]) %*% t((W[indices,i])) %*% invOmega)
    
    gammaGradient = gammaGradient - 2 * t(O_T %*% dOmega %*% phi0[indices, , i] %*% diag(gradLambda0^2))
    
    lambdaGradient = lambdaGradient - diag(t(B1) %*% dOmega %*% B1) * 2 * gradLambda0
    
    tauGradient = tauGradient - sum(diag(dOmega)) * 2 * gradTau
  }
  
  return(-as.matrix(c(as.vector(t(gammaGradient)), lambdaGradient,tauGradient)))
}

#' @export
tx <- function(x) { matrix(x, L, p) }

#' @export
formattedOutput = function(parameters, GammaTall = TRUE, paramsLong = FALSE){
  Gamma = t(tx(parameters[1:(p*L)]))
  if(!GammaTall){
    Gamma = array(parameters[1:(p*L)], c(p, L))
  }
  lambda = parameters[(p*L+1):(p*(L+1))]
  if(paramsLong){
    lambda= parameters[(p*(L+4)+1):(p*(L+5))]
  }
  Gamma = Gamma[order(abs(lambda)),]
  lambda = sort(abs(lambda))
  tau = parameters[(p*(L+1)+1)]
  if(paramsLong){
    tau= parameters[(p*(L+5)+2)]
  }
  listOut = list("Gamma" = Gamma, "lambda" = lambda, "tau" = tau)
  return(listOut)
}

#' @export
removeMeanFn = function(data, noTimePoints, n){
  longwArray = cbind(1:noTimePoints,as.vector(data)) #coerces 1:noTimePoints to repeat n times
  longwArray = longwArray[complete.cases(longwArray),]
  dataLowess = lowess(longwArray)
  
  predictlowess = function(x){
    xArray = dataLowess$x
    yArray = dataLowess$y
    if(!(x %in% xArray)){
      return(NA)
    }
    index = min(which(xArray == x))
    return(yArray[index])
  }
  
  lowessmeancolumns = apply(as.matrix(1:noTimePoints),1,predictlowess)
  lowessmeanmatrix = array(rep(lowessmeancolumns, n),c(noTimePoints,n))
  
  mean0wArray = data - lowessmeanmatrix
  
  return(mean0wArray)
}


## Code for reshape input longitudinal
#' @export
convertToMatrix = function(tallData, n, noTimePoints){
  # Talldata format: subjNo, timeNo, Val
  wArray = array(rep(NA, (n * noTimePoints)), c(noTimePoints, n))
  counter = 1
  while(counter <= dim(tallData)[1]){
    wArray[tallData[counter,2],tallData[counter,1]] = tallData[counter,3]
    counter = counter + 1
  }
  return(wArray)
}

sillyData = cbind(c(1,1,1,1,2,2,3,3,3),c(1,2,5,10,2,4,4,5,6),rnorm(9))
convertToMatrix(sillyData, 3, 10) # Works!

#' @export
estimateCovariance = function(estimatedObj, s = NULL, t = NULL, asMat = FALSE){
  ## s, t, should lie in interval [0,1]
  ## Current estimates based on 51 discrete time points and nearest discrete time
  
  
  
  nMax2 = 51;
  maxGrid2 = seq (0, 1, length=nMax2);
  op2 = create.bspline.basis (nbasis=L);
  bMat2 = eval.basis (maxGrid2, op2);
  OBASIS = orthonormalization (bMat2, basis=F);
  OBASIS = OBASIS*sqrt(nMax2);
  
  O_T=t(OBASIS)
  phi0temp = t(estimatedObj$estGamma %*% O_T)
  B1 = phi0temp
  Omega = B1 %*% diag(estimatedObj$estLambda^2) %*% t(B1) + estimatedObj$estTau^2 * diag(51)
  
  
  
  if(!asMat){
    sApprox = round(s*50)+1
    tApprox = round(t*50)+1
    return(Omega[sApprox,tApprox] / 51)
  } else if(asMat){
    return(Omega)
  }
} 

### Big estimation function

#' @export
estimateProcess = function(noTimePoints, data, nBasisFns = 10, nEigFns = 5, meanZeroFlag = FALSE, trtFlag = FALSE, method = "PSD", seedValue = NULL, noIter = NULL, tol = NULL){
  
  L = nBasisFns
  p = nEigFns
  
  
  if(meanZeroFlag){
    W0 = data
  } else{
    W0 = removeMeanFn(data, noTimePoints = noTimePoints, n = n)
  }
  
  # # Testing if variance function looks like data
  #   trueCovarianceFunction = OBASIS %*% t(Gamma0) %*% lambda0 %*% t(OBASIS %*% t(Gamma0) %*% lambda0) + tau^2 * diag(noTimePoints)
  #   estimatedCovarianceFunction = OBASIS %*% t(estGamma) %*% estLambda %*% t(OBASIS %*% t(estGamma) %*% estLambda) + estTau^2 * diag(noTimePoints)
  #   lines(sqrt(diag(trueVarianceFunction)), col = "red")
  #   lines(sqrt(diag(estimatedVarianceFunction)), col = "green")
  #   lines(apply(data, 1, sd, na.rm = TRUE), col = "blue")
  # #
  #   
  # sum((trueCovarianceFunction-estimatedCovarianceFunction)^2)/noTimePoints^2  
  # 
  # 
  
  
  
  ####
  
  nMax2 = noTimePoints;
  maxGrid2 = seq (0, 1, length=nMax2);
  op2 = create.bspline.basis (nbasis=L);
  bMat2 = eval.basis (maxGrid2, op2);
  OBASIS = orthonormalization (bMat2, basis=F);
  OBASIS = OBASIS*sqrt(nMax2);
  
  if(!is.null(seedValue)){
    set.seed(seedValue)
  }
  
  if(is.null(noIter)){
    noIter = round((20/log(n)^2)*600)
  }
  
  if(method == "PSD"){
    library(Rcpp)
    library(RcppArmadillo)
    #sourceCpp("testrcpp11NoStageWorking.cpp")
    
    if(is.null(tol)){
      tol = 0.01
    }
    
    stepsizes = c(0,4^(5:10) / 1000000)
    
    startParameters1 = projGradient(rnorm(((L+5)*p + 3)))
    estCands1 = orderParams(newDescentFullOnlyNoStage(W0, as.matrix(OBASIS), startParameters1, 0, p, noIter, as.matrix(stepsizes), tol))
    
    startParameters2 = projGradient(rnorm(((L+5)*p + 3)))
    estCands2 = orderParams(newDescentFullOnlyNoStage(W0, as.matrix(OBASIS), startParameters2, 0, p, noIter, as.matrix(stepsizes), tol))
    
    startParameters3 = projGradient(rnorm(((L+5)*p + 3)))
    estCands3 = orderParams(newDescentFullOnlyNoStage(W0, as.matrix(OBASIS), startParameters3, 0, p, noIter, as.matrix(stepsizes), tol))
    
    allCands = cbind(estCands1, estCands2, estCands3)
    
    # Choose parameters which minimize negative log-likelihood
    estimatedParameters = allCands[,which.min(allCands[(p*(L+5)+4),])]
    
    estGamma = array(estimatedParameters[1:(L*p)], c(p,L))
    estLambda = estimatedParameters[(p*(L+4)+1):(p*(L+5))]
    estTau = estimatedParameters[(p*(L+5)+2)]
    logLikelihood = estimatedParameters[(p*(L+5)+4)]
  } else if(method == "RTRSR1"){
    library(ManifoldOptim)
    in.w = W0
    print(dim(in.w))
    
    if(is.null(tol)){
      tol = 1
    }
    
    ####
    #### Note, need to define simpleLogLikelihood2 and simpleGradient2ndTry2 within this 
    #### part of statement because need it to recognize scoping of variable in.w
    #### and also need function to be just of parameters for manifoldOptim
    ####
    
    simpleLogLikelihood2 <- function(parameters){
      W = in.w
      n = dim(W)[2]
      likeGamma0 = t(tx(parameters[1:(p*L)]))
      
      
      likeLambda0 = parameters[((L)* p + 1): ((L+1)*p)]
      likeTau = parameters[((L+1)*p + 1)]
      
      # Empty array holds all m phi0 functions at nMax2 times and n observations 
      phi0 = array(rep(NA, nMax2 * p * n), c(nMax2, p, n))
      
      likely = 0
      
      for (i in 1:n){
        
        times = sum(!is.na(W[,i]))       # number of observed time points
        indices = (!is.na(W[,i])) * 1:nMax2  # Exact value time points on grid
        indices = indices[indices != 0]
        
        B_T=t(OBASIS[indices,])
        phi0temp = t(likeGamma0 %*% B_T)
        
        phi0[indices, , i] = phi0temp
        
        B1 = phi0[indices, , i]
        
        Omega = B1 %*% diag(likeLambda0^2) %*% t(B1) + likeTau^2 * diag(times)
        
        likely = likely + dmvnorm(W[indices,i], mean = rep(0,times), sigma = Omega, log = TRUE)
      }
      return(-likely)
      
    }
    
    simpleGradient2ndTry2 <- function(parameters){
      W = in.w
      n = dim(W)[2]
      # Gamma = array(parameters[1:(p*L)], c(L, p))
      # Gamma = t(Gamma)
      # 
      # gradGamma0 = array(Gamma[1:(p*L)], c(p, L))
      
      gradGamma0 = t(tx(parameters[1:(p*L)]))
      
      gradLambda0 = parameters[(L* p + 1): ((L+1)*p)]
      gradTau = parameters[((L+1)*p + 1)]
      
      # Empty array holds all m phi0 functions at nMax2 times and n observations 
      phi0 = array(rep(NA, nMax2 * p * n), c(nMax2, p, n))
      gammaGradient = array(rep(0, (p*L)), c(p, L))
      lambdaGradient = numeric(p)
      tauGradient = 0
      
      
      for (i in 1:n){
        times = sum(!is.na(W[,i]))       # number of observed time points
        indices = (!is.na(W[,i])) * 1:nMax2  # Exact value time points on grid
        indices = indices[indices != 0]
        
        O_T=t(OBASIS[indices,])
        phi0temp = t(gradGamma0 %*% O_T)
        
        phi0[indices, , i] = phi0temp
        
        B1 = phi0[indices, , i]
        
        
        Omega = B1 %*% diag(gradLambda0^2) %*% t(B1) + gradTau^2 * diag(times)
        
        #print(Omega)
        
        invOmega = solve(Omega)
        
        dOmega = .5 * (invOmega - invOmega %*% (W[indices,i]) %*% t((W[indices,i])) %*% invOmega)
        
        gammaGradient = gammaGradient - 2 * t(O_T %*% dOmega %*% phi0[indices, , i] %*% diag(gradLambda0^2))
        
        lambdaGradient = lambdaGradient - diag(t(B1) %*% dOmega %*% B1) * 2 * gradLambda0
        
        tauGradient = tauGradient - sum(diag(dOmega)) * 2 * gradTau
      }
      
      return(-as.matrix(c(as.vector(t(gammaGradient)), lambdaGradient,tauGradient)))
    }
    
    mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
    prob <- new(mod$RProblem, simpleLogLikelihood2, simpleGradient2ndTry2)
    
    
    mani.defn2 <- get.product.defn(get.stiefel.defn(L, p), get.euclidean.defn(n = (p+1), m = 1))
    
    mani.params <- get.manifold.params(IsCheckParams = TRUE)
    solver.params <- get.solver.params(DEBUG = 0, Tolerance = tol,
                                       Max_Iteration = 1000, IsCheckParams = TRUE, Stop_Criterion = 1)
    
    randParameters = projGradient(rnorm(((L+5)*p + 3)))
    x0 = c(as.vector(t(array(randParameters[1:(p*L)], c(p,L)))),
           randParameters[(p*(L+4)+1):(p*(L+5))], randParameters[(p*(L+5)+2)])
    print(simpleLogLikelihood2(trueParams2))
    simpleLogLikelihood2(x0)
    
    #RTRSR1
    
    ptm = proc.time()
    res <- manifold.optim(prob, mani.defn2, x0 = x0, method = "RTRSR1",
                          mani.params = mani.params, solver.params = solver.params)
    print(proc.time()-ptm)
    
    print(simpleLogLikelihood2(trueParams2))
    simpleLogLikelihood2(x0)
    #simpleLogLikelihood2(x0Mat)
    print(simpleLogLikelihood2(res$xopt))
    #res
    
    probConverged = (simpleLogLikelihood2(res$xopt) - 1.01*simpleLogLikelihood2(trueParams2) < 0)
    print("Prob Converged?")
    print(probConverged)
    print("Time?")
    print(res$elapsed)
    formatOut = formattedOutput(res$xopt)
    estGamma = formatOut$Gamma
    estLambda = formatOut$lambda
    estTau = formatOut$tau
    logLikelihood = simpleLogLikelihood2(res$xopt)
    time = res$elapsed
  }
  
  paramList = list("estGamma" = estGamma, "estLambda" = estLambda, "estTau" = estTau, "logLikelihood" = logLikelihood, "time" = time)
  return(paramList)
}


#result = estimateProcess(noTimePoints = 51, data = train$w, nBasisFns = 10, nEigFns = 5, meanZeroFlag = TRUE, trtFlag = FALSE, method = "RTRSR1", seedValue = NULL, noIter = 300, tol = 0.01)
#estimateCovariance(result, s=.2, t=.2, asMat = FALSE)

#covMatEst = estimateCovariance(result, asMat = TRUE)