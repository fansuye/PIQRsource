#' Residual Projection for Quantile Regression in Vertically Partitioned Big Data
#'
#' @param X The design matrix (with intercept)
#' @param Y The response vector
#' @param M The number of partitions
#' @param tau The quantile of interest
#' @param eps The tolerance parameter for convergence (the PIQR algorithm)
#' @param epsm The tolerance parameter for convergence (the subproblems for residual extraction)
#' @param maxstep Maximum number of iterations allowed for the PIQR algorithm
#' @return The coefficient estimation and the number of iterations
#' @examples
#' gcov = function(p, rho){
#'   cov = matrix(1, p, p);
#'     for(i in 1:p){
#'       for(j in 1:p){
#'         if(i < j) cov[i,j] = rho^{j-i}
#'         else cov[i,j] = cov[j,i]
#'       }
#'     }
#'   cov
#' }
#'
#' n = 10000
#' p = 500
#' rho = 0.5
#' M = 20
#' tau = 0.7
#' epsM = 1e-03
#' eps = 1e-02
#' maxstep = 5000
#' set.seed(66)
#' X = matrix(rt(n*p, 3, ncp = 0), n, p)
#' cov = gcov(p, rho)
#' X = X%*%chol(cov)
#' Xint = cbind(1, X)
#' beta0 = rnorm(1)
#' beta = rnorm(p)
#' e = rnorm(n)
#' Y = beta0+X%*%beta+e
#' 
#' beta_true = c(beta0+qnorm(tau), beta)
#' result_PIQR = PIQR(Xint, Y, M, tau, eps, epsM, maxstep)
#' beta_PIQR = result_PIQR$beta
#' AE_PIQR = sum(abs(beta_PIQR-beta_true))
#' Iteration_PIQR = result_PIQR$iteration
#' @export
#'


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
