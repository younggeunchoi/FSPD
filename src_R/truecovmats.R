ftn.mat.OvrlpBdiag = function(p) {

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

  mat.BDIAG = diag(rep(0.6, p)) + 0.4*blockmat + 0.4*neardiagmat ;
  mat.BDIAG = (mat.BDIAG + t(mat.BDIAG))/2 ;

  return(mat.BDIAG) ;
}


## true cov structures
mat.ID = diag(rep(1, p)) ;
rho = 0.7 ; mat.AR = rho^(abs(outer(1:p, 1:p, "-"))) ;

mat.taper = abs(outer(1:p, 1:p, "-"))/10 ;
mat.taper = 1 - mat.taper ;
mat.taper[mat.taper < 0] = 0 ;

mat.OvrlpBdiag = ftn.mat.OvrlpBdiag(p) ;

mat.ARlike2 = abs(outer(1:(p/2), 1:(p/2), "-"))/10 ;
mat.ARlike2 = 1 - mat.ARlike2 ;
mat.ARlike2[mat.ARlike2 < 0] = 0 ;
mat.ARlike2 = bdiag(mat.ARlike2, diag(rep(4, p/2))) ;

set.seed(30) ;
mat.unif = matrix(runif(n = p^2, min=0.3, max=1) ,nrow = p) ;
mat.sparseind = floor(runif(n=200) * p^2) + 1 ;
mat.sparse.unif = mat.unif ;
mat.sparse.unif[-mat.sparseind] = 0 ;
mat.sparse.unif[upper.tri(mat.sparse.unif)] = 0 ;
mat.sparse.unif = mat.sparse.unif + t(mat.sparse.unif) ;
diag(mat.sparse.unif) = 1 ;

delta = max( c( -min(eigen(mat.sparse.unif)$values), 0) )  + 0.01 ;
mat.sparse.unif = mat.sparse.unif + diag(rep(delta, p)) ;
mat.sparse = mat.sparse.unif / max(mat.sparse.unif) ;


rm(mat.sparse.unif, mat.unif) ;