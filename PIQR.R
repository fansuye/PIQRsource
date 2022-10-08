install.packages("quantreg")
library(quantreg)

PIQR = function(X, Y, M, tau, eps, epsm, maxstep){
  
  p = dim(X)[2]-1
  pm = p/M
  distance = 1
  iteration = 0
  time = 0
  R = Y
  X1 = X[,1:(pm+1)]
  beta_PIQR = beta_delta = rep(0, p+1)
  while((distance > eps)&(iteration < maxstep)){
    timeM = rep(0, M)
    ptm = proc.time()
    result_PIQRM = rq.fit.fnb(X1, R/M, tau, eps = epsm)
    beta_delta[1:(pm+1)] = result_PIQRM$coefficients
    R = R-X1%*%beta_delta[1:(pm+1)]
    beta_PIQR[1:(pm+1)] = beta_PIQR[1:(pm+1)]+beta_delta[1:(pm+1)]
    timeM[1] = (proc.time() - ptm)[3]
    for(m in 2:M){
      Xm = X[,((m-1)*pm+2):(m*pm+1)]
      ptm = proc.time()
      result_PIQRM = rq.fit.fnb(Xm, R/M, tau, eps = epsm)
      beta_delta[((m-1)*pm+2):(m*pm+1)] = result_PIQRM$coefficients
      R = R-Xm%*%beta_delta[((m-1)*pm+2):(m*pm+1)]
      beta_PIQR[((m-1)*pm+2):(m*pm+1)] = beta_PIQR[((m-1)*pm+2):(m*pm+1)]+beta_delta[((m-1)*pm+2):(m*pm+1)]
      timeM[m] = (proc.time() - ptm)[3]
    }
    time = time+max(timeM)
    distance = sum(abs(beta_delta))
    iteration = iteration+1
  }
  return(list(beta = beta_PIQR, iteration = iteration))
  
}