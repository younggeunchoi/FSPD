### FUNCTION LIST
# concord.CV
# concord.IV
# symlasso.IV
# glasso.IV
# clime.IV



## FOR CONCORD ESTIMATION

concord.CV = function(data, lambdaseq = 1.2^seq(from=10, to=-30, by=-1),
                      n.fold=5, faster=TRUE) {
  time1 = Sys.time() ;
  #out all results!
  n = nrow(data) ; p = ncol(data) ;
  BICseq = c() ; sparsity = c() ;

  n.test = floor(n / n.fold) ;

  errorstack = matrix(0, nrow=n.fold, ncol=length(lambdaseq)) ;
  sparsitystack = matrix(0, nrow=n.fold, ncol=length(lambdaseq)) ;
  ### for each training/test set
  ### randomly choose the validation sets
  set.seed(100) ;
  perm.ind = sample.int(n, size = n, replace = FALSE, prob = NULL) ;


  for (lamind in 1:length(lambdaseq)) {

    for (foldind in 1:n.fold) {
      test.ind = ((foldind-1)*n.test + 1) : (foldind*n.test) ;
      test.ind.perm = perm.ind[test.ind] ;
      data.test = data[test.ind.perm, ] ;
      data.train = data[-test.ind.perm, ] ;

      ccrd.train.lambda = concord(data=data.train, lambda=lambdaseq[lamind]) ;
      RSS = 0 ;
        for (j in 1:p) {
          weightvec = ccrd.train.lambda[-j,j]/ccrd.train.lambda[j,j] ;
          predictvec = data.test[ ,-j] %*% weightvec ;
          RSS = RSS +
               sum( (as.numeric(data.test[ ,j]) + as.numeric(predictvec))^2 ) ;
        }
      sparsity = ftn.sparsity(ccrd.train.lambda) ;
      sparsitystack[foldind, lamind] = sparsity ;
      errorstack[foldind, lamind] = RSS ;
    }
    errorseq = as.numeric(colMeans(errorstack)) ;
    
    if (faster) {
        if (lamind > 1) { if (errorseq[1] < errorseq[lamind]) break ; }
    }
    
  }
  errorseq = as.numeric(colMeans(errorstack)) ;
  sparsityseq = as.numeric(colMeans(sparsitystack)) ;

  lambda.opt = lambdaseq[which.min(errorseq[errorseq > 0])] ;
  concord.opt = concord(data=data, lambda=lambda.opt) ;
  time2 = difftime(Sys.time(), time1, units = "secs") ;

  return(list(concord.opt=concord.opt, lambda.opt=lambda.opt,
              lambdaseq=lambdaseq, error=errorseq, sparsity=sparsityseq,
              comp.time=time2)) ;
}


concord.IV = function(data, lambdaseq = 1.2^seq(from=10, to=-30, by=-1)*n/300,
                      popn.sqrtcov, diststr="mvn", seed=175, faster=TRUE) {
  time1 = Sys.time() ;
  #out all results!
  n = nrow(data) ; p = ncol(data) ;
  sparsity = c() ;


  #generate validation data
  std.data.val = gen.datamat(SAMPSIZE=n, DIM=p, DISTSTR=diststr) ;
  data.val = std.data.val %*% popn.sqrtcov ;
  data.val = apply(data.val, 2, scale, center=TRUE, scale=FALSE) ;
  #sampcov.val = var(data.val) ;
  
  errorseq = rep(0, length(lambdaseq)) ;
  sparsityseq = rep(0, length(lambdaseq)) ;
  ccrd.train = list() ;

  for (lamind in 1:length(lambdaseq)) {

    ccrd.train[[lamind]] = concord(data=data, lambda=lambdaseq[lamind]) ;
    ccrd.train.lambda = ccrd.train[[lamind]] ;
    RSS = 0 ;
    for (j in 1:p) {
      weightvec = ccrd.train.lambda[-j,j]/ccrd.train.lambda[j,j] ;
      predictvec = data.val[ ,-j] %*% weightvec ;
      RSS = RSS +
           sum( (as.numeric(data.val[ ,j]) + as.numeric(predictvec))^2 ) ;
    }
    sparsity = ftn.sparsity(ccrd.train.lambda) ;
    sparsityseq[lamind] = sparsity ;
    errorseq[lamind] = RSS ;

    if (faster) {
        if (lamind > 10) { if (errorseq[1] < errorseq[lamind]) break ; }
    }
   
  }

  lambda.opt.ind = which.min(errorseq[errorseq > 0]) ;
  lambda.opt = lambdaseq[lambda.opt.ind] ;
  concord.opt = ccrd.train[[lambda.opt.ind]] ;
  nz.opt = sparsityseq[lambda.opt.ind]
  iter.opt = attr(concord.opt, "iter")
  tol.opt = attr(concord.opt, "tol")

  rm(ccrd.train) ;
  time2 = difftime(Sys.time(), time1, units = "secs") ;

  return(list(estimate=concord.opt, lambda.opt=lambda.opt, nz.opt = nz.opt,
              lambdaseq=lambdaseq, errorseq=errorseq, nzseq=sparsityseq,
              time=time2, iter=iter.opt, tol=tol.opt)) ;
}


symlasso.IV = function(data, lambdaseq = 1.2^seq(from=10, to=-30, by=-1)*n/300,
                      popn.sqrtcov, diststr="mvn", seed=175, faster=TRUE) {
  time1 = Sys.time() ;
  #out all results!
  n = nrow(data) ; p = ncol(data) ;
  sparsity = c() ;


  #generate validation data
  std.data.val = gen.datamat(SAMPSIZE=n, DIM=p, DISTSTR=diststr) ;
  data.val = std.data.val %*% popn.sqrtcov ;
  data.val = apply(data.val, 2, scale, center=TRUE, scale=FALSE) ;
  #sampcov.val = var(data.val) ;
  
  errorseq = rep(0, length(lambdaseq)) ;
  sparsityseq = rep(0, length(lambdaseq)) ;
  slasso.train = list() ;

  for (lamind in 1:length(lambdaseq)) {

    slasso.train[[lamind]] = symlasso(data=data, lambda=lambdaseq[lamind]) ;
    slasso.train.lambda = slasso.train[[lamind]] ;
    RSS = 0 ;
    for (j in 1:p) {
      weightvec = slasso.train.lambda[-j,j]/slasso.train.lambda[j,j] ;
      predictvec = data.val[ ,-j] %*% weightvec ;
      RSS = RSS +
           sum( (as.numeric(data.val[ ,j]) + as.numeric(predictvec))^2 ) ;
    }
    sparsity = ftn.sparsity(slasso.train.lambda) ;
    sparsityseq[lamind] = sparsity ;
    errorseq[lamind] = RSS ;

    if (faster) {
        if (lamind > 10) { if (errorseq[1] < errorseq[lamind]) break ; }
    }
   
  }

  lambda.opt.ind = which.min(errorseq[errorseq > 0]) ;
  lambda.opt = lambdaseq[lambda.opt.ind] ;
  symlasso.opt = slasso.train[[lambda.opt.ind]] ;
  nz.opt = sparsityseq[lambda.opt.ind]
  iter.opt = attr(symlasso.opt, "iter")
  tol.opt = attr(symlasso.opt, "tol")

  rm(slasso.train) ;
  time2 = difftime(Sys.time(), time1, units = "secs") ;

  return(list(estimate=symlasso.opt, lambda.opt=lambda.opt, nz.opt = nz.opt,
              lambdaseq=lambdaseq, errorseq=errorseq, nzseq=sparsityseq,
              time=time2, iter=iter.opt, tol=tol.opt)) ;
}



glasso.IV = function(data, lambdaseq = 1.2*((0.01/1.2)^(1/49))^(0:49),
                      popn.sqrtcov, diststr="mvn", seed=175, is.cor=TRUE, is.trace=FALSE) {
  time1 = Sys.time() ;
  #out all results!
  n = nrow(data) ; p = ncol(data) ;
  sparsity = c() ;
  cov.samp = t(data) %*% data / n + diag(rep(1e-7, p)) ;
  if (is.cor) {
    sdvec.samp = sqrt(diag(cov.samp)) ;
    cor.samp = diag(1/sdvec.samp) %*% cov.samp %*% diag(1/sdvec.samp) ;
  }

  #generate validation data
  std.data.val = gen.datamat(SAMPSIZE=n, DIM=p, DISTSTR=diststr) ;
  data.val = std.data.val %*% popn.sqrtcov ;
  data.val = apply(data.val, 2, scale, center=TRUE, scale=FALSE) ;
  cov.samp.val = t(data.val) %*% data.val / n ;
  #sampcov.val = var(data.val) ;
  
  errorseq = rep(0, length(lambdaseq)) ;
  sparsityseq = rep(0, length(lambdaseq)) ;
  if (is.cor) {
    glassop.list = glassopath(s=cor.samp, rholist=lambdaseq, penalize.diagonal=FALSE,
      trace=is.trace, maxit=100, thr=1e-4) ;
  } else {
    glassop.list = glassopath(s=cov.samp, rholist=lambdaseq, penalize.diagonal=FALSE,
                              trace=is.trace, maxit=100, thr=1e-4) ;
  }
  glassop.train = glassop.list$wi ;  
  errflag = glassop.list$errflag ;

  for (lamind in 1:length(lambdaseq)) {

    if (errflag[lamind]) next ;
    glasso.train.lambda = glassop.train[ , ,lamind] ;
    glasso.train.lambda = ( glasso.train.lambda + t(glasso.train.lambda) ) / 2 ;
    if (is.cor) glasso.train.lambda = diag(1/sdvec.samp) %*% glasso.train.lambda %*% diag(1/sdvec.samp) ;

    eigval.full.train = eigen(glasso.train.lambda, only.values=TRUE)$values ;
    eigval.pos.train = eigval.full.train[eigval.full.train > 0] ;
    pred.negloglik = sum(glasso.train.lambda * cov.samp.val) -
                 sum(log(eigval.pos.train)) ;
    if (pred.negloglik == 0) pred.negloglik = Inf ;

    errorseq[lamind] = pred.negloglik ;
    sparsity = ftn.sparsity(glasso.train.lambda) ;
    sparsityseq[lamind] = sparsity ;
  
  }

  lambda.opt.ind = which.min(errorseq[errorseq > 0]) ;
  lambda.opt = lambdaseq[lambda.opt.ind] ;
  grplasso.opt = glassop.train[ , ,lambda.opt.ind] ;
  grplasso.opt = (t(grplasso.opt) + grplasso.opt) /2 ;
  if (is.cor) grplasso.opt = diag(1/sdvec.samp) %*% grplasso.opt %*% diag(1/sdvec.samp) ;

  nz.opt = sparsityseq[lambda.opt.ind]
  iter.opt = NA
  tol.opt = NA

  rm(glassop.train) ;
  time2 = difftime(Sys.time(), time1, units = "secs") ;

  return(list(estimate=grplasso.opt, lambda.opt=lambda.opt, nz.opt = nz.opt,
              lambdaseq=lambdaseq, errorseq=errorseq, nzseq=sparsityseq,
              time=time2, iter=iter.opt, tol=tol.opt)) ;
}



clime.IV = function(data, lambdaseq = 1.2*((0.001/1.2)^(1/99))^(0:99),
                      popn.sqrtcov, diststr="mvn", seed=175, faster=TRUE) {
  time1 = Sys.time() ;
  #out all results!
  n = nrow(data) ; p = ncol(data) ;

  #generate validation data
  std.data.val = gen.datamat(SAMPSIZE=n, DIM=p, DISTSTR=diststr) ;
  data.val = std.data.val %*% popn.sqrtcov ;
  #popn.cov = popn.sqrtcov %*% popn.sqrtcov ;
  sampcov.val = cov(data.val) * ((n-1)/n) ;
  
  errorseq = rep(0, length(lambdaseq)) ;
  sparsityseq = rep(0, length(lambdaseq)) ;
  clime.train = list() ;

  tempcnt = 0 ;
  for (lamind in 1:length(lambdaseq)) {
  
    clime.train[[lamind]] = clime(x = data, lambda = lambdaseq[lamind])$Omegalist[[1]]
  
    #clime.train[[lamind]] = sugm(data = data, lambda = lambdaseq[lamind], nlambda = NULL,
    #      lambda.min.ratio = NULL, rho = NULL, method = "clime", sym = "or",
    #      shrink=NULL, prec = 1e-4, max.ite = 1000, standardize = FALSE,
    #      perturb = FALSE, verbose = FALSE)$icov[[1]] ;
    clime.train.lambda = clime.train[[lamind]] ;

    eigval.full.train = eigen(clime.train.lambda, only.values=TRUE)$values ;
    eigval.pos.train = eigval.full.train[eigval.full.train > 0] ;
    pred.negloglik = sum(clime.train.lambda * sampcov.val) -
                 sum(log(eigval.pos.train)) ;
    if (pred.negloglik == 0) pred.negloglik = Inf ;
    errorseq[lamind] = pred.negloglik ;

    
    sparsity = ftn.sparsity(clime.train.lambda) ;
    sparsityseq[lamind] = sparsity ;
    
    if (faster) { if (lamind > 10) { 
       if (errorseq[1] < errorseq[lamind]) {
          print(lamind) ;
          break ; 
       }
    }}
   
  }

  lambda.opt.ind = which.min(errorseq[errorseq > 0]) ;
  lambda.opt = lambdaseq[lambda.opt.ind] ;
  clime.opt = clime.train[[lambda.opt.ind]] ;
  rm(clime.train) ;
  time2 = difftime(Sys.time(), time1, units = "secs") ;

  return(list(prec.opt=clime.opt, lambda.opt=lambda.opt,
              lambdaseq=lambdaseq, error=errorseq, sparsity=sparsityseq,
              comp.time=time2)) ;
}




