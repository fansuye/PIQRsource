# Source code for PIQR algorithm

## Instructions

It can be used to reproduce the simulation studies in the following paper:

**Ye Fan, Jr-Shin Li and Nan Lin**. *Residual Projection for Quantile Regression in Vertically Partitioned Big Data.*

## Code for testing PIQR
```
###function for generating the correlation matrix
gcov = function(p, rho){
  cov = matrix(1, p, p);
  for(i in 1:p){
    for(j in 1:p){
      if(i < j) cov[i,j] = rho^{j-i}
      else cov[i,j] = cov[j,i]
    }
  }
  cov
}


n = 10000
p = 500
rho = 0.5
epsM = 1e-03
eps = 1e-02
maxstep = 5000
M = 20
set.seed(66)
beta0 = rnorm(1)
beta = rnorm(p)
X = matrix(rt(n*p, 3, ncp = 0), n, p)
cov = gcov(p, rho)
X = X%*%chol(cov)
Xint = cbind(1, X)
'''
