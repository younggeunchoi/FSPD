## FUNCTION LIST
# thres.hard
# thres.soft
# thres.scad
# subftn.taper
# d.thres.soft : assigns 1 or -1 for positive or negative elements
# offdiagthres.ftn
# ftn.posi.mat
# subftn.unithres
# ftn.cov.uni
# subftn.adapthres
# ftn.cov.adap
# subftn.Xue
# ftn.cov.Xue
# ftn.cov.Roth




thres.hard = function(obj, cut) {
  ## Note that the threshold can varies upon each element
  res = obj ;
  res[abs(res) < cut] = 0 ;
  return(res) ;
}

thres.soft = function(obj, cut) {
  ## Note that the threshold can varies upon each element
  res = abs(obj) - cut ;
  res[res < 0] = 0 ;
  res = res*sign(obj) ;
  return(res) ;
}

thres.scad = function(obj, cut) {
  M1 = (abs(obj) < 2*cut) * sign(obj) *(abs(obj)-cut) *(abs(obj)>cut);
  M2 = (abs(obj) < 3.7*cut) * (abs(obj)>=2*cut) *((3.7-1)*obj-3.7*sign(obj) *cut)/(3.7-2);
  M3 = (abs(obj) >= 3.7*cut) * obj;
  return(M1 + M2 + M3);
}


subftn.taper = function(sampcov, bandwidth, band=FALSE) {
  p = ncol(sampcov) ;
  h = bandwidth ;
  if (h %% 2 != 0) cat("Warning! bandwidth recommended even\n") ;
  if (h > p) cat("Error! bandwidth cannot be broder than p\n") ;
  # weightvec : 1, 1, 1 (by k-th), ..lineary decay by h-th, 0, 0, 0 ....
  if (!band) {
    k = floor(h / 2) ;
    decayvec = seq(from = 1, to = 0, length = h - k + 2) ;
    decayvec = decayvec[-c(1, length(decayvec))] ;
    weightvec = c(rep(1, k), decayvec, rep(0, p - h)) ;
  }
    if (band) {
    weightvec = c(rep(1, h), rep(0, p - h)) ;

  }

  weightmat = toeplitz(weightvec) ;
  return(sampcov * weightmat) ;

}


## function :: universal-soft-threshholding opeartor
ftn.cov.taper = function(datamat, bandwidthseq=NULL, is.band=FALSE, 
                       bandnum=100) {
  #bandwidthseq have priority

  ## (cutpoint estimation)
  ## split the samples
  n = nrow(datamat) ; p = ncol(datamat) ;
  ## centerizing
  data = apply(datamat, 2, scale, center=T, scale=F) ;
  cov.samp = t(data) %*% data / (n-1) ;
  
  if (is.null(bandwidthseq)) bandwidthseq = floor(seq(from = p, to = 1, by = bandnum)) ;

  time1 = Sys.time() ;
  n.test = floor(n / 5) ;
  n.fold = 5 ;
#  J = floor(1 / sqrt(log(p)/n)) ;
  ## make the stack of errors
#  bandwidthseq = sqrt(log(p)/n) * (0:J) ;

  errorseq = rep(0, length(bandwidthseq)) ;
  nzseq = rep(0, length(bandwidthseq)) ;
  ### for each training/test set
  ### randomly choose the validation sets
  set.seed(1) ;
  perm.ind = sample.int(n, size = n, replace = FALSE, prob = NULL) ;

  for (i in 1:n.fold) {
    test.ind = ((i-1)*n.test + 1) : (i*n.test) ;
    test.ind.perm = perm.ind[test.ind] ;
    data.test = data[test.ind.perm, ] ;
    data.train = data[-test.ind.perm, ] ;
    ### test sample cov matrix and training sample cov matrix
    sampcov.test = var(data.test) ;
    sampcov.train = var(data.train) ;

    #### for each candidate : grid on the integer multiplication of log(p)/n
    for ( j in 1:length(bandwidthseq) ) {
      #### training (thresholded) estimator
      #function(datamat, sampcov, thres.ftn, cut, diag.thres=FALSE)
      cov.train = subftn.taper(sampcov=sampcov.train, band = is.band,
                              bandwidth = bandwidthseq[j]) ;
      #### frobenius norm loss
      error = sum((cov.train - sampcov.test)^2) ;
      errorseq[j] = errorseq[j] + error ;
      nzseq[j] = nzseq[j] + ftn.sparsity(cov.train) ;
      }
    }
  #time2 = Sys.time() ;
  #print(time2 - time1) ;
  ## (select the optimal cutpoint)
  nzseq = nzseq / n.fold ;
  band.opt = bandwidthseq[which.min(errorseq)] ;

  ## (store the thresholding estimator)
  resmat = subftn.taper(sampcov=cov.samp, band = is.band,
                              bandwidth=band.opt) ;
  time2 = difftime(Sys.time(), time1, units="secs") ;
  return(list(estimate = resmat, lambda.opt = band.opt, cutmat = NA,
              nz.opt = ftn.sparsity(resmat), lambdaseq=bandwidthseq, errorseq = errorseq,
              nzseq = nzseq, time = time2, iter = NA, tol = NA)) ;
}

d.thres.soft = function(obj, cut) {
  res = abs(obj) - cut ;
  res[res > 0] = 1 ;
  res[res < 0] = 0 ;
  res = res*sign(obj) ;
  return(res) ;  
}


offdiagthres.ftn = function(mat, cut, thres.ftn) {

  resmat = thres.ftn(mat, cut = cut) ;
  diag(resmat) = diag(mat) ; 

  return(resmat) ;
}


ftn.posi.mat = function(mat, mineigval) {
### Projection onto the PSD cone
#  if(!isSymmetric(mat)) {
#    if(isSymmetric(mat, tol=10^-3)) { mat = (mat + t(mat)) / 2 ; } else {
#    stop("Input is not symmetric")  ; }
#  }
  
  eig = eigen(mat) ;
  eigval = eig$values ;
  if (sum(eigval < mineigval) == 0) return(list(res=mat, operated=c(FALSE))) ;
  eigval[eigval < mineigval] = mineigval ;
  eigvec = eig$vectors ;
  
  res = eigvec %*% diag(eigval) %*% t(eigvec) ;
  res = (res + t(res)) / 2 ;  # to stabilize symmetricity
  return(list(res=res, operated=c(TRUE))) ; 
}                                               






subftn.unithres = function(datamat=NULL, sampcov, thres.ftn=thres.soft, cut, diag.thres=FALSE) {

  resmat = thres.ftn(obj = sampcov, cut = cut) ;
  if (!diag.thres) diag(resmat) = diag(sampcov) ;
  return(resmat) ;
}

## function :: universal-soft-threshholding opeartor
ftn.cov.uni = function(datamat, thres.ftn=thres.soft, thresseq=NULL,
                       thresnum=100, diag.thres=FALSE) {
  #thresseq have priority

  ## (cutpoint estimation)
  ## split the samples
  n = nrow(datamat) ; p = ncol(datamat) ;
  ## centerizing
  data = apply(datamat, 2, scale, center=T, scale=F) ;
  cov.samp = t(data) %*% data / (n-1) ;
  cov.samp.off = cov.samp ; diag(cov.samp.off) = 0 ;
  max.thres = max(abs(cov.samp.off)) ;
  thresseq.temp = 1.05^seq(from=0, to=thresnum, by=1) - 1 ;
  if (is.null(thresseq)) thresseq = max.thres * thresseq.temp/max(thresseq.temp) ;

  time1 = Sys.time() ;
  n.test = floor(n / 5) ;
  n.fold = 5 ;
#  J = floor(1 / sqrt(log(p)/n)) ;
  ## make the stack of errors
#  thresseq = sqrt(log(p)/n) * (0:J) ;

  errorseq = rep(0, length(thresseq)) ;
  nzseq = rep(0, length(thresseq)) ;
  ### for each training/test set
  ### randomly choose the validation sets
  set.seed(1) ;
  perm.ind = sample.int(n, size = n, replace = FALSE, prob = NULL) ;

  for (i in 1:n.fold) {
    test.ind = ((i-1)*n.test + 1) : (i*n.test) ;
    test.ind.perm = perm.ind[test.ind] ;
    data.test = data[test.ind.perm, ] ;
    data.train = data[-test.ind.perm, ] ;
    ### test sample cov matrix and training sample cov matrix
    sampcov.test = var(data.test) ;
    sampcov.train = var(data.train) ;

    #### for each candidate : grid on the integer multiplication of log(p)/n
    for ( j in 1:length(thresseq) ) {
      #### training (thresholded) estimator
      #function(datamat, sampcov, thres.ftn, cut, diag.thres=FALSE)
      cov.train = subftn.unithres(datamat=NULL, sampcov=sampcov.train, thres.ftn=thres.ftn,
                              cut=thresseq[j], diag.thres=diag.thres) ;
      #### frobenius norm loss
      error = sum((cov.train - sampcov.test)^2) ;
      errorseq[j] = errorseq[j] + error ;
      nzseq[j] = nzseq[j] + ftn.sparsity(cov.train) ;
      }
    }
  #time2 = Sys.time() ;
  #print(time2 - time1) ;
  ## (select the optimal cutpoint)
  nzseq = nzseq / n.fold ;
  thres.opt = thresseq[which.min(errorseq)] ;
  cutmat = matrix(thres.opt, p, p) ;
  if (!diag.thres) diag(cutmat) = 0 ;

  ## (store the thresholding estimator)
  resmat = subftn.unithres(datamat=NULL, sampcov=cov.samp, thres.ftn=thres.ftn,
                              cut=thres.opt, diag.thres=diag.thres) ;
  time2 = difftime(Sys.time(), time1, units="secs") ;
  return(list(estimate = resmat, lambda.opt = thres.opt, cutmat = cutmat,
              nz.opt = ftn.sparsity(resmat), lambdaseq=thresseq, errorseq = errorseq,
              nzseq = nzseq, time = time2, iter = NA, tol = NA)) ;
}
## demo
#res = ftn.cov.uni(datamat = data, thres.ftn = thres.soft)
#frobenius.norm(as.matrix(res - cov.true))
#spectral.norm(as.matrix(res - cov.true))
#ftn.PR(estimate = res, true = cov.true)




## function :: adaptive-soft-threshholding opeartor
## FOR ADAPTIVE THRESHOLDING

subftn.adapthres = function(datamat, sampcov, thres.ftn=thres.soft, cut, diag.thres=FALSE) {
  n = nrow(datamat) ; p = ncol(datamat) ;
  datamat.centered = datamat - rep(colMeans(datamat), rep(n, p)) ;
  datamat.centered = datamat.centered^2 ;
  ## sqrt(theta_ij)
  stdevmat.sampcov = sqrt( t(datamat.centered) %*% datamat.centered / n  - 
                      (n - 2) / n * sampcov^2 ) ;
#  print(stdevmat.sampcov)
  cutmat = cut * stdevmat.sampcov ;
  if(!diag.thres) diag(cutmat) = 0 ;

#  print(cutmat) ;
  resmat = thres.ftn(obj = sampcov, cut = cutmat) ;
  attr(resmat, 'cutmat') <- cutmat ;
  return(resmat) ;
}
# standard cut : 2 * sqrt(log(p)/n)

## function :: adaptive-soft-threshholding opeartor
ftn.cov.adap = function(datamat, thres.ftn=thres.soft,
  thresseq = NULL, thresnum = 100, diag.thres=FALSE) {
  thresseq.temp = 1.05^seq(from=0, to=thresnum, by=1) - 1 ;
  if (is.null(thresseq)) thresseq = 4 * thresseq.temp/max(thresseq.temp) ;

  time1 = Sys.time() ;
  
  ## (cutpoint estimation)
  ## split the samples
  n = nrow(datamat) ; p = ncol(datamat) ;
  n.test = floor(n / 5) ;
  n.fold = 5 ;
  
  ## HERE IS DIFFERENT FROM OTHERS
  ## make the stack of errors
  
  errorseq = rep(0, length(thresseq)) ;
  nzseq = rep(0, length(thresseq)) ;
  ### for each training/test set
  ### randomly choose the validation sets
  set.seed(1) ;
  perm.ind = sample.int(n, size = n, replace = FALSE, prob = NULL) ;
  
  for (i in 1:n.fold) {
    test.ind = ((i-1)*n.test + 1) : (i*n.test) ;
    test.ind.perm = perm.ind[test.ind] ;
    data.test = datamat[test.ind.perm, ] ;
    data.train = datamat[-test.ind.perm, ] ;
    ### test sample cov matrix and training sample cov matrix
    sampcov.test = var(data.test) ;
    sampcov.train = var(data.train) ;
    
    for ( j in 1:length(thresseq) ) {
      #### training (thresholded) estimator
      cov.train = subftn.adapthres(datamat = data.train, sampcov = sampcov.train,
       thres.ftn = thres.ftn, cut = thresseq[j] * sqrt(log(p) / n), diag.thres=diag.thres )  ;
      #### frobenius norm loss
      error = sum((as.numeric(cov.train - sampcov.test))^2) ;
      errorseq[j] = errorseq[j] + error ;
      nzseq[j] = nzseq[j] + ftn.sparsity(cov.train) ;
    }
  }
  #time2 = Sys.time() ;
  #print(time2 - time1) ;
  ## (select the optimal cutpoint)
  nzseq = nzseq / n.fold ;
  thres.opt = thresseq[which.min(errorseq)] ;


  ## (store the thresholding estimator)
  resmat = subftn.adapthres(datamat = datamat, sampcov = var(datamat),
   thres.ftn = thres.ftn, cut = thres.opt * sqrt(log(p) / n), diag.thres=diag.thres)  ;
  time2 = difftime(Sys.time(), time1, units="secs") ;
  return(list(estimate = resmat, lambda.opt = thres.opt, cutmat = attr(resmat, 'cutmat'),
              nz.opt = ftn.sparsity(resmat), lambdaseq=thresseq, errorseq = errorseq,
              nzseq = nzseq, time = time2, iter = NA, tol = NA)) ;                           
}
#res = ftn.cov.adap(datamat = data, thres.ftn = thres.soft)
#frobenius.norm(as.matrix(res - cov.true))



## CHECK whether it is already positive definite or not!
subftn.Xue = function(sampcov, initial.cov.thres, cut, thres.ftn=thres.soft,
                     mineigval=10^-3, diag.thres=FALSE, tolprimal=10^-5,
                     toldual=10^-4, MAXITER=1000) {

  #### MAIN ALGORITHM
  tuningmat.old = 0 ;
  mu = 1 ;
  covmat.old = initial.cov.thres ;
  Soff = ( sum(abs(sampcov)) - sum(abs(diag(sampcov))) ) / 2 ;

  errorvec1 = rep(0, MAXITER) #####  
  errorvec2 = rep(0, MAXITER) #####    
  for (iter in 1:MAXITER) {
    Thetamat.new.list = ftn.posi.mat(mat = covmat.old + mu * tuningmat.old,
              mineigval = mineigval) ;
    # if ((iter == 1) & !(Thetamat.new.list$operated)) {
    #   return(list(est=covmat.old, iter=0, tol=tol)) ; }
    Thetamat.new = Thetamat.new.list$res ;
    
    temp = mu*(sampcov - tuningmat.old) + Thetamat.new ;
    covmat.new = thres.ftn( temp, cut*mu ) / (1 + mu) ;
    if (!diag.thres) diag(covmat.new) = diag(temp) / (1 + mu) ;
        
    tuningmat.new = tuningmat.old - (Thetamat.new - covmat.new)/mu ;
    
    #### convergence criteria
    err.normalizer = Soff ;               
    error1 = sum(abs(covmat.new - covmat.old)) / err.normalizer ;
   error2 = sum(abs(covmat.new - Thetamat.new)) / err.normalizer ;
    errorvec1[iter] = error1 ; #######
    errorvec2[iter] = error2 ; #######    
    #if ( error1  < tol) {
    if ( (error1 < tolprimal) & (error2 < toldual) ) {
#      cat("Funtion subftn.Xue : convergence acheived at step", iter, "\n") ;
#      windows() ; plot.ts(errorvec1, ylim=c(0,10^-5)) ;
      break ;
    } else {
     covmat.old = covmat.new ;
     tuningmat.old = tuningmat.new ;
    }
    if (iter == MAXITER) { 
#      cat("Funtion subftn.Xue : reached to maximum iteration\n") ;
#    windows() ; plot.ts(log10(errorvec1), ylim=c(-6,1)) ;
#      windows() ; plot.ts(log10(errorvec2), ylim=c(-6,1)) ;
#      windows() ; plot.ts(errorvec2, ylim=c(0,10^-3)) ;
    }
  }

  attr(covmat.new, 'iter') <- iter ;
  attr(covmat.new, 'tol') <- tolprimal ;
  return(covmat.new) ;
}



subftn.Roth = function(sampcov, initial.cov.thres, cut, thres.ftn=thres.soft,
 barrier=10^-2, diag.thres=FALSE, tolin = 10^-8, tolprimal=10^-7, toldual=10^-7,
 MAXITIN = 10000, MAXITOUT=1000) {
  tau = barrier ;
  Soff = ( sum(abs(sampcov)) - sum(abs(diag(sampcov))) ) / 2 ;
  p = ncol(sampcov) ;
  if ( is.null(dim(cut)[1])  ) { lammat = matrix(cut, p, p) ; } else { lammat = cut ; }
  if ( !diag.thres ) diag(lammat) = 0 ;
  if ( sum(dim(lammat) == c(p,p)) < 2 ) cat("Warning: no digonal thresholding, but sampcov[1,1] != initial.cov.thres[1,1]\n") 

  if ( !diag.thres & ( abs(sampcov[1,1] - initial.cov.thres[1,1]) > 1e-6 ) ) {
    cat("Warning: no digonal thresholding, but sampcov[1,1] != initial.cov.thres[1,1]\n") 
  }
  # initial value
  oe = eigen(initial.cov.thres, symmetric = TRUE)
  evs = oe$val/2 + sqrt(oe$val^2 + 4 * tau)/2
  s0 = tcrossprod(oe$vec * rep(evs, each = p), oe$vec)
  i0 = tcrossprod(oe$vec * rep(1/evs, each = p), oe$vec)
  Ioff = ( sum(abs(i0)) - sum(abs(diag(i0))) ) / 2 ;

  #### MAIN ALGORITHM
  covmat.old = s0 ;
  precmat.old = i0 ;
  
  errorvec1 = rep(0, MAXITOUT) #####  
  errorvec2 = rep(0, MAXITOUT)
  dyn.load("d:/LSPD/subftn_lasso.dll") ;
  for (iter in 1:MAXITOUT) {
    print(iter) ;
    covmat.new = covmat.old ;
    precmat.new = precmat.old ;

    for (j in 1:p) {
      ##### step 1
      covmat.new[j,j] = sampcov[j,j] + tau * precmat.old[j,j] - lammat[j,j] ;  

      betainitvec = covmat.new[-j,j] ;
      Vmat = diag(1,p-1) + precmat.old[-j,-j] * (tau / covmat.new[j,j]) ;
      uvec = sampcov[-j,j] ;
      lamvec = lammat[-j,j] ;
      # all size (p-1)
      res = .C("subftn_lasso", betavec = as.numeric(betainitvec), Vmat = as.numeric(Vmat), 
       uvec = as.numeric(uvec), lamvec = as.numeric(lamvec), veclen = as.integer(p-1), 
       reltol = as.numeric(tolin), iter = as.integer(MAXITIN)) ;
      resbeta = res[[1]] ;
      resiter = res[[7]] ;
      covmat.new[-j,j] = resbeta ;
      covmat.new[j,-j] = resbeta ;
      #void subftn_lasso(double* betavec, double* Vmat, double* uvec, 
      # double* lamvec, int* veclen, double* reltol, int* iter)

      ##### step 2
      precmat.new[-j,j] = precmat.old[-j,-j] %*% (resbeta / covmat.new[j,j]) ;
      precmat.new[j,-j] = t(precmat.new[-j,j]) ;

      ##### step 3
      precmat.new[j,j] = (1 - sum(resbeta * as.numeric(precmat.new[-j,j]))) / covmat.new[j,j] ;

    }

    #### convergence criteria
    
    error1 = sum(abs(covmat.new - covmat.old)) / Soff ;
    error2 = sum(abs(precmat.new - precmat.old)) / Ioff ;
    errorvec1[iter] = error1 ; #######
    errorvec2[iter] = error2 ; #######
    print(c(error1, error2)) ;    
    #if ( error1  < tol) {
    if ( (error1 < tolprimal) & (error2 < toldual) ) {
      #      cat("Funtion subftn.Xue : convergence acheived at step", iter, "\n") ;
      #      windows() ; plot.ts(errorvec1, ylim=c(0,10^-5)) ;
      break ;
    } else {
       covmat.old = covmat.new ;
       precmat.old = precmat.new ;
    }
    if (iter == MAXITOUT) { 
      cat("Funtion subftn.Roth : reached to maximum iteration\n") ;
      #    windows() ; plot.ts(log10(errorvec1), ylim=c(-6,1)) ;
      #      windows() ; plot.ts(log10(errorvec2), ylim=c(-6,1)) ;
      #      windows() ; plot.ts(errorvec2, ylim=c(0,10^-3)) ;
    }
  # iter ends  
  }

  attr(covmat.new, 'iter') <- iter ;
  attr(covmat.new, 'tol') <- c(tolin, tolprimal, toldual) ;
  dyn.unload("d:/LSPD/subftn_lasso.dll") ;
  return(covmat.new) ;
}



subftn.Roth2 = function(sampcov, initial.cov.thres, cut, thres.ftn=thres.soft,
 barrier=10^-2, diag.thres=FALSE, tolin = 10^-8, tolprimal=10^-7, toldual=10^-7,
 MAXITIN = 10000, MAXITOUT=1000) {
  tau = barrier ;
  Soff = ( sum(abs(sampcov)) - sum(abs(diag(sampcov))) ) / 2 ;
  p = ncol(sampcov) ;
  if ( is.null(dim(cut)[1])  ) { lammat = matrix(cut, p, p) ; } else { lammat = cut ; }
  if ( !diag.thres ) diag(lammat) = 0 ;
  if ( sum(dim(lammat) == c(p,p)) < 2 ) cat("Warning: no digonal thresholding, but sampcov[1,1] != initial.cov.thres[1,1]\n") 

  if ( !diag.thres & ( abs(sampcov[1,1] - initial.cov.thres[1,1]) > 1e-6 ) ) {
    cat("Warning: no digonal thresholding, but sampcov[1,1] != initial.cov.thres[1,1]\n") 
  }
  # initial value
  oe = eigen(initial.cov.thres, symmetric = TRUE)
  evs = oe$val/2 + sqrt(oe$val^2 + 4 * tau)/2
  s0 = tcrossprod(oe$vec * rep(evs, each = p), oe$vec)
  i0 = tcrossprod(oe$vec * rep(1/evs, each = p), oe$vec)
  Ioff = ( sum(abs(i0)) - sum(abs(diag(i0))) ) / 2 ;

  #### MAIN ALGORITHM
  covmat.old = s0 ;
  precmat.old = i0 ;
  
  errorvec1 = rep(0, MAXITOUT) #####  
  errorvec2 = rep(0, MAXITOUT)
  dyn.load("d:/LSPD/subftn_lasso.dll") ;
  for (iter in 1:MAXITOUT) {
    print(iter) ;
    covmat.temp = sampcov + tau * precmat.old ;
    covmat.new = thres.ftn(covmat.temp, lammat) ;
    precmat.new = solve(covmat.new) ;

    #### convergence criteria
    
    error1 = sum(abs(covmat.new - covmat.old)) / Soff ;
    error2 = sum(abs(precmat.new - precmat.old)) / Ioff ;
    errorvec1[iter] = error1 ; #######
    errorvec2[iter] = error2 ; #######
    print(c(error1, error2)) ;    
    #if ( error1  < tol) {
    if ( (error1 < tolprimal) & (error2 < toldual) ) {
      #      cat("Funtion subftn.Xue : convergence acheived at step", iter, "\n") ;
      #      windows() ; plot.ts(errorvec1, ylim=c(0,10^-5)) ;
      break ;
    } else {
       covmat.old = covmat.new ;
       precmat.old = precmat.new ;
    }
    if (iter == MAXITOUT) { 
      cat("Funtion subftn.Roth : reached to maximum iteration\n") ;
      #    windows() ; plot.ts(log10(errorvec1), ylim=c(-6,1)) ;
      #      windows() ; plot.ts(log10(errorvec2), ylim=c(-6,1)) ;
      #      windows() ; plot.ts(errorvec2, ylim=c(0,10^-3)) ;
    }
    # iter ends  
    print(ftn.eig(covmat.new)) ;
  }

  attr(covmat.new, 'iter') <- iter ;
  attr(covmat.new, 'tol') <- tolprimal ;
  attr(covmat.new, 'prec') <- precmat.new ;
  return(covmat.new) ;
}

# ## function :: Xue, Ma, Zou (2012, JASA), positve l1-penalty
# ftn.cov.Xue = function(datamat, thres.ftn, thresseq,
#      mineigval, tol = 10^-7, MAXITER = 1000, message=FALSE) {
#   time1 = Sys.time() ;
# ## subftn.Xue = function(datamat, sampcov=var(datamat), cut, thres.ftn,
# ##                     mineigval = 10^-4, tol = 10^-4, MAXITER = 1000)

#   ## split the samples
#   n = nrow(datamat) ; p = ncol(datamat) ;
#   n.test = floor(n / 5) ;
#   n.fold = 5 ;
# #  J = floor(1 / sqrt(log(p)/n)) ;
#   ## make the stack of errors
# #  thresseq = sqrt(log(p)/n) * (0:J) ;
# #  thresseq = seq(from=0, to=4, by=0.1) ;
#   errorstack = rep(0, length(thresseq)) ;
#   ### for each training/test set
#   ### randomly choose the validation sets
#   set.seed(1) ;
#   perm.ind = sample.int(n, size = n, replace = FALSE, prob = NULL) ;
#   if (message) { cat("permutated indices :\n") ; print(perm.ind) ; }
#   ### tuningmat(Lambda)=zero, mu=2, covmat(Sigma)

#   ### split 5-fold
#   for (i in 1:n.fold) {
#     if (message) cat("Test set : fold", i, "\n")  ;    
#     test.ind = ((i-1)*n.test + 1) : (i*n.test) ;
#     test.ind.perm = perm.ind[test.ind] ;
#     data.test = datamat[test.ind.perm, ] ;
#     data.train = datamat[-test.ind.perm, ] ;
#     ### test sample cov matrix and training sample cov matrix
#     sampcov.test = var(data.test) ;
#     sampcov.train = var(data.train) ;

#     ### training estimate    
#     #### for each candidate : grid on the integer multiplication of log(p)/n
#     for ( j in 1:length(thresseq) ) {
#       if (message) cat("threshold", thresseq[j], "\n")  ;
#       cov.train = subftn.Xue(datamat = data.train,
#                     sampcov = sampcov.train,
#                     cut = thresseq[j], thres.ftn = thres.ftn,
#                     mineigval = mineigval, tol = tol, MAXITER = MAXITER)$est ;
      
#       #### frobenius norm loss
#       error = frobenius.norm(cov.train - sampcov.test)^2 ;
#       errorstack[j] = errorstack[j] + error ;
#     }
#   }
#   if (message) print(rbind(thresseq, errorstack)) ;
#   #time2 = Sys.time() ;
#   #print(time2 - time1) ;
#   ## (select the optimal cutpoint)
#   thres.opt = thresseq[which.min(errorstack)] ;

#   ## (store the thresholding estimator)
#   cov.Xue = subftn.Xue(datamat = datamat,
#                     sampcov = cov.samp,
#                     cut = thres.opt, thres.ftn = thres.ftn,
#                     mineigval = mineigval, tol = tol, MAXITER = MAXITER) ;

#   time2 = difftime(Sys.time(), time1, units="secs") ;
#   return(list(est = cov.Xue$est, thres=thres.opt, iter = cov.Xue$iter,
#               tol = cov.Xue$tol, time = time2)) ;  
# }
# ##demo
# #res = ftn.cov.Xue(datamat=data, thres.ftn=thres.soft)
# #frobenius.norm(as.matrix(res - cov.true))
# #spectral.norm(as.matrix(res - cov.true))
# #ftn.PR(estimate = res, true = cov.true)
# #


ftn.cov.Roth = function(datamat, thresseq, tol = 10^-7, MAXITER = 1000) {
  time1 = Sys.time() ;
  roth = pdsoft.cv(x = datamat, lam.vec = thresseq, standard = FALSE,
  	init="soft", nsplits = 5, n.tr = 4, tolout = tol, maxitout=MAXITER,
    quiet=FALSE) ;
  time2 = difftime(Sys.time(), time1, units="secs") ;
  return(list(est = roth$sigma, iter = NA, tol = tol,
  		 thres = roth$best.lam, time = time2)) ;
}