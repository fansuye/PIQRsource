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
#' gcov = function(p, rho, type){
#'   if(type == "exchangeable"){
#'     cov = matrix(rho, p, p)
#'     diag(cov) = rep(1, p)
#'   }
#'   else{
#'     cov = diag(p)
#'     for(i in 1:p){
#'       for(j in 1:p){
#'         if(i < j) cov[i,j] = rho^{j-i}
#'         else cov[i,j] = cov[j,i]
#'       }
#'     }
#'   }
#'   cov
#' }
#'
#' N = 10000
#' p = 100
#' n = 10
#' rep = rep(n, N)
#' nsum = sum(rep)
#' d = 0.75*p
#' rho_X = 0.5
#' rho_e = 0.5
#' tau = 0.75
#' set.seed(999)
#' X = matrix(rnorm(nsum*p), nsum, p)
#' cov_X = gcov(p, rho_X, "ar1")
#' X = X%*%chol(cov_X)
#' for(i in 1:d){
#'   X[,i] = pnorm(X[,i])
#' }
#' set.seed(1000)
#' e = matrix(rt(N*n, 3), N, n)
#' cov_e = gcov(n, rho_e, "ar1")
#' e = as.vector(t(e%*%chol(cov_e)))
#' sigma = 0.5
#' e = sigma*e
#' beta0 = rnorm(1)
#' beta = rnorm(p)
#' Y = beta0+X%*%beta+apply(X[,1:d]*e/d, 1, sum)
#' beta_true = c(beta0, quantile(e/d, tau)+beta[1:d], beta[(d+1):p])
#'
#' WQR = WQRADMM(X, Y, rep, tau, TRUE, "WQR")
#' beta_WQR = WQR$Estimation_WQR
#' AE_WQR = sum(abs(beta_WQR-beta_true))
#' Time_WQR = WQR$Time_WQR
#' Time_WQRADMM = WQR$Time_total
#' @export
#'



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
