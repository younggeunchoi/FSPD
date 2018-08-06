## FUNCTION LIST
# ftn.eig
# ftn.sparsity
# subftn.LS


####################################
########## FUNCTION PART ###########
####################################


ftn.eig = function(symmat) {
  symmat = (symmat + t(symmat)) / 2 ;
  eigval = eigen(symmat, only.values=TRUE)$values ;
  return( c( eig.max = max(eigval), eig.min=min(eigval) ) ) ;
}
ftn.sparsity = function(symmat) {
  return( sum(abs(symmat) > 10^-6) / length(symmat) ) ;
}



subftn.LS = function(cov.raw, mineigval, method="FSopt", message=FALSE) {
  eps = mineigval ;
  p = nrow(cov.raw) ;
  eigen.raw = eigen(cov.raw, only.values = TRUE) ;
  eigval.raw = eigen.raw$values ;
  lmean = mean(eigval.raw) ;
  lmin = min(eigval.raw) ;
  lmax = max(eigval.raw) ;
  
  # calculate mu
  if (method == "Sopt") target = (lmin + lmax) / 2 ;
  if (method == "Max") target = lmax ;
  if (method == "Fopt") {
    target = lmean +
      var(eigval.raw) * (p-1) / p / (lmean - lmin) ;
  }
  if (method == "FSopt") {
    target1 = (lmin + lmax) / 2 ;
    target2 = lmean +
      var(eigval.raw) * (p-1) / p / (lmean - lmin) ;
    target = max(c(target1, target2)) ;
  }
    
  # conduct shrinkage
  if (lmin < eps) {  
    if (method == "Infty") {
      return(  cov.raw + diag( rep(eps - lmin, p) )  )  ;
    }
  
    alpha = 1 - (eps - lmin)/(target - lmin) ;
    alpha = max( alpha, 0 ) ;
    if (message) print(alpha)
    cov.LS.onemd = diag(rep((1-alpha)*target, p)) + alpha*cov.raw ;
    return(cov.LS.onemd) ;
  }
  
  if (message) cat("Subfunction subftn.LS : no shrinkage\n") ;
  return(cov.raw) ;
}






subftn.sqrt = function(cov.raw, mineigval=NULL) {
  ### Projection onto the PSD cone
  #  if(!isSymmetric(mat)) {
  #    if(isSymmetric(mat, tol=10^-3)) { mat = (mat + t(mat)) / 2 ; } else {
  #    stop("Input is not symmetric")  ; }
  #  }
  eig = eigen(cov.raw) ;
  eigval = eig$values ;
  if (!is.null(mineigval)) eigval[eigval < mineigval] = mineigval ;
  eigval.sqrt = sqrt(eigval) ;
  eigvec = eig$vectors ;
  res = eigvec %*% diag(eigval.sqrt) %*% t(eigvec) ;
  res = (res + t(res)) / 2 ;  # to stabilize symmetricity
  return(res) ;   
}





