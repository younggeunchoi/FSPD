## vignette.R
## a working example to use the proposed method.


## set working directory as that stores covsim.R
setwd("/Users/cyg/Dropbox/Codes/FSPD_public/")
source("./src_R/functions.R") 
source("./src_R/cov_estimators.R") ;

## choose epsilon (cut-off point to do PD-modification)
eps = 10^-2 ;

## sample size and data dimension
n = 100 ;
p = 100 ;


# generate true covariance matrix (M1: linearly tapered)
mat.taper = abs(outer(1:p, 1:p, "-"))/10 ;
mat.taper = 1 - mat.taper ;
mat.taper[mat.taper < 0] = 0 ;
cov.true = mat.taper ;
eig.cov.true = eigen(cov.true) ;
eigvec.cov.true = eig.cov.true$vectors ;
eigval.cov.true = eig.cov.true$values ;
sqrt.cov.true = eigvec.cov.true %*% diag(sqrt(eigval.cov.true)) %*% t(eigvec.cov.true) ;
maxeigval.cov.true = max(eigval.cov.true) ;
mineigval.cov.true = min(eigval.cov.true) ;


# generate data (multivariate t-distn)
set.seed(40) ;
y = matrix(rnorm(n*p), nrow=n, ncol=p) ;
u = rchisq(n, df=5) ;
data.gen = sqrt(3/5) * y / sqrt(u/5) ;
data.gen = data.gen %*% sqrt.cov.true ;
## centering data
data = apply(data.gen, 2, scale, center=TRUE, scale=FALSE) ;
## sample covariance matrix
cov.samp = t(data) %*% data / (n-1) ;
    

## Ex 1. soft thresholding estimator (universal thresholding)
## name: "uni"
## function "ftn.cov.uni" calculates the estimator
uni.list = ftn.cov.uni(datamat = data) ;
uni = uni.list$estimate ;

# the resulting estimate is not PD.
print(tail(eigen(uni, only.values=TRUE)$values))
# apply FSPD.
uniL = subftn.LS(cov.raw = uni, mineigval = eps)
# now it is PD.
print(tail(eigen(uniL, only.values=TRUE)$values))

# compare estimation error
# matrix l1 norm
c(norm(as.matrix(uni - cov.true), type="O"), 
  norm(as.matrix(uniL - cov.true), type="O"))
# spectral norm
c(norm(as.matrix(uni - cov.true), type="2"), 
  norm(as.matrix(uniL - cov.true), type="2"))
# Frobenius norm
c(norm(as.matrix(uni - cov.true), type="F"), 
  norm(as.matrix(uniL - cov.true), type="F"))
