library(matrixcalc) ;
library(Matrix) ;
n = 100 ;
p = 400 ;

### cov : AR(1), i want to test a arbitrary sparsity, too.

# make the true covariance structure
cov.ID = diag(rep(1, p)) ;

rho = 0.7 ;
cov.AR = rho^(abs(outer(1:p, 1:p, "-"))) ;

cov.ARlike = abs(outer(1:p, 1:p, "-"))/10 ;
cov.ARlike = 1 - cov.ARlike ;
cov.ARlike[cov.ARlike < 0] = 0 ;

rind.mat = matrix(rep(1:p, p), nrow=p) ;
cind.mat = t(rind.mat) ;
blocksize = floor(p/20) ;
blocknum = floor(p / blocksize) + 1 ;
blockmat = kronecker(Diagonal(blocknum), matrix(1, nrow=blocksize, ncol=blocksize)) ;
blockmat = blockmat[1:p, 1:p] ;

tempmat = NULL ;
for (k in 1:(blocknum-1)) {
  temp1 = matrix(FALSE, nrow=blocknum*blocksize, ncol=blocksize-1) ;
  temp2 = c( rep(FALSE, k*blocksize) , rep(TRUE, blocksize) , rep(FALSE, blocksize*(blocknum-k-1)) ) ;
  tempmat = cbind(tempmat, temp1, temp2) ;
}
tempmat = tempmat[1:p, 1:p] ;
neardiagmat = tempmat + t(tempmat) ;

cov.BDIAG = diag(rep(0.6, p)) + 0.6*blockmat + 0.4*neardiagmat ;
cov.BDIAG = (cov.BDIAG + t(cov.BDIAG))/2 ;


ftn.errcompar = function(n, cov.true, MAXREPL=100) {
  p = ncol(cov.true) ;
  eig.cov.true = eigen(cov.true) ;
  eigvec.true = eig.cov.true$vectors ;
  eigval.true = eig.cov.true$values ;
  sqrt.true = eigvec.true %*% diag(sqrt(eigval.true)) %*% t(eigvec.true) ; 
  

  
  
  # function for performance comparison
  perf.output = function(cov.test) {
    ## matrix 1-norm
    MatrixL1 = max(colSums(abs(cov.test - cov.true)))
    ## operator norm : maxinum eigenvalue
    MatrixL2 = max(eigen(cov.test - cov.true)$values)
    ## frobenius norm
    Frob = frobenius.norm(cov.test - cov.true)
    ## largest eigenvalue difference (abs)
    eigen.test = eigen(cov.test) ;
    eigvec.test = eigen.test$vectors ;
    eigval.test = eigen.test$values ;
    MaxEigvalDiff = abs( max(eigval.test) - max(eigval.true) ) ;
    ## least eigenvalue difference (abs)
    MinEigvalDiff = abs( min(eigval.test) - min(eigval.true) ) ;
    ## 1st pc cos(angle difference) (abs)
    Cos1stPCAngle = abs ( as.numeric(t(eigvec.test[ ,1]) %*% eigvec.true[ ,1]) ) ;
    return(c(MatrixL1=MatrixL1, MatrixL2=MatrixL2,
             Frob=Frob, MaxEigvalDiff=MaxEigvalDiff,
             MinEigvalDiff=MinEigvalDiff, Cos1stPCAngle=Cos1stPCAngle)) ;
  }
  
  errstack.samp = errstack.thres = errstack.shr = errstack.thres.shr = NULL ;
  time1 = Sys.time() ;
  for (repl in 1:MAXREPL) {
    
    # generate multiv normal data from mean 0
    data = matrix(rnorm(n*p), nrow=n, ncol=p) ;    
    data = data %*% sqrt.true ;
    cov.samp = var(data) ;
    
    
    
    # estimator 1 : give thresholding estimator ####################
    ## (cutpoint estimation)
    ## split the samples
    n.test = floor(n / log(n)) ;
    n.fold = floor(n / n.test) ;
    J = floor(1 / (log(p)/n)) ;
    ## make the stack of errors
    thresseq = log(p)/n * (0:J) ;
    errorstack = rep(0, length(thresseq)) ;
    ### for each training/test set
    
    for (i in 1:n.fold) {
      test.ind = ((i-1)*n.test + 1) : (i*n.test) ;
      data.test = data[test.ind, ] ;
      data.train = data[-test.ind, ] ;
      ### test sample cov matrix and training sample cov matrix
      cov.test = var(data.test) ;
      cov.train = var(data.train) ;
      
      cov.train.old = cov.train ;
      #### for each candidate : grid on the integer multiplication of log(p)/n
      for ( j in 1:length(thresseq) ) {
        #### training (thresholded) estimator
        cov.train.old[abs(cov.train.old) < thresseq[j]] = 0 ;
        #### frobenius norm loss
        error = frobenius.norm(cov.train.old - cov.test)^2 ;
        errorstack[j] = errorstack[j] + error ;
      }
    }
    #time2 = Sys.time() ;
    #print(time2 - time1) ;
    
    
    ## (select the optimal cutpoint)
    thres.opt = thresseq[which.min(errorstack)] ;
    
    ## (store the thresholding estimator)
    cov.thres = cov.samp ;
    cov.thres[abs(cov.thres) < thres.opt] = 0 ;
    
    
    
    
    # estimator 2 : give shrinkage estimator ####################
    ## eigenvalue decomposition
    #ftn.shr = function(cov.samp, data, thres=NULL) {
    #  n = nrow(data) ; p = ncol(data) ;
    eigen.samp = eigen(cov.samp) ;
    eigval.samp = eigen.samp$values ;
    
    m.n = mean(eigval.samp) ;
    d.n.sq = mean(eigval.samp^2) - m.n^2 ;
    
    tempstack = rep(0, n) ;
    for (k in 1:n) {
      temp = as.vector(data[k, ]) ;
      temp2 = temp %*% t(temp) ;
      #    if (!is.null(thres)) temp2[abs(temp2) < thres] = 0 ;
      tempstack[k] = frobenius.norm(temp2 - cov.samp)^2/p ;  
    }
    bar.b.n.sq = mean(tempstack)/n ;
    b.n.sq = min(c(bar.b.n.sq, d.n.sq)) ;
    a.n.sq = d.n.sq - b.n.sq ;
    
    cov.shr = diag(rep(b.n.sq / d.n.sq * m.n, p)) + a.n.sq / d.n.sq * cov.samp ;
    #  return(cov.shr) ;
    #}
    #cov.shr = ftn.shr(cov.samp, data) ;
    
    
    # estimator 3 : linear shrinkage after thresholding ####################
    eigen.thres = eigen(cov.thres) ;
    eigval.thres = eigen.thres$values ;
    mean.eigval.thres = mean(eigval.thres) ;
    min.eigval.thres = min(eigval.thres) ;
    max.eigval.thres = max(eigval.thres) ;
    eps = 10^-2 ;
    cov.thres.shr = cov.thres ;
    #eps = 10^-1 ;
    if (min.eigval.thres < eps) {
      alpha = 1 - (eps - min.eigval.thres)/(mean.eigval.thres - min.eigval.thres) ;
      cov.thres.shr = diag(rep((1-alpha)*mean.eigval.thres, p)) + alpha*cov.thres ;  
      
    }
    
    
    ## stack the errors
    
    errstack.samp = cbind(errstack.samp, perf.output(cov.samp)) ;
    errstack.thres = cbind(errstack.thres, perf.output(cov.thres)) ;
    errstack.shr = cbind(errstack.shr, perf.output(cov.shr)) ;
    errstack.thres.shr = cbind(errstack.thres.shr, perf.output(cov.thres.shr)) ;
    
  }
  
  merr.samp = rowMeans(errstack.samp) ;
  merr.thres = rowMeans(errstack.thres) ;
  merr.shr = rowMeans(errstack.shr) ;
  merr.thres.shr = rowMeans(errstack.thres.shr) ;
  
  sderr.samp = apply(errstack.samp, 1, sd) ;
  sderr.thres = apply(errstack.thres, 1, sd) ;
  sderr.shr = apply(errstack.shr, 1, sd) ;
  sderr.thres.shr = apply(errstack.thres.shr, 1, sd) ;
  
  time3 = Sys.time() ;
  print(time3-time1) ;

  cat("True largest & smallest true eigenvalues\n") ;
  print( c( max(eigval.true), min(eigval.true) ) ) ;
  
  return(data.frame(mean.samp=merr.samp, sd.samp=sderr.samp,
                    mean.shr=merr.shr, sd.shr = sderr.shr,
                    mean.thres=merr.thres, sd.thres = sderr.thres,
                    mean.thres.shr=merr.thres.shr, sd.thres.shr=sderr.thres.shr)
  )
}



round(ftn.errcompar(n=100, cov.true=cov.ID), 2)
round(ftn.errcompar(n=100, cov.true=cov.AR, MAXREPL=50), 2)
round(ftn.errcompar(n=100, cov.true=cov.ARlike, MAXREPL=50), 2)
round(ftn.errcompar(n=100, cov.true=data.matrix(cov.BDIAG), MAXREPL=50), 2)


# frobenius.norm(cov.thres.shr2 - cov.true)
# max(eigen(cov.thres.shr2 - cov.true)$values)
# max(eigen(cov.thres.shr2)$values)
# min(eigen(cov.thres.shr2)$values)

## obtain eigenvalues of thresholding estimator
## eigenvalue average
## shrink : note that the smallest eigenvalue must be >= epsilon





