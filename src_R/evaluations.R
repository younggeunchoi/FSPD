ftn.PR = function(estimate, true) {
  num.entry = length(true) ;
  true.P = sum( (abs(true) > 10^-6) ) ;
  true.nonP = num.entry - true.P ;

  FP = sum(  ((abs(estimate) > 10^-6) & (abs(true) <= 10^-6))  ) ;
  TP = sum(  ((abs(estimate) > 10^-6) & (abs(true) > 10^-6))  ) ;

  return(c(FPR = FP / true.nonP, TPR = TP / true.P)) ;
}
#ftn.PR(estimate = matrix(c(1,2,3,4,5,6), nrow=2),
#       true = matrix(c(1,2,3,4,0,0), nrow=2))


## make use of norm(x, type = c("O", "I", "F", "M", "2"))



ftn.store = function(GENINFO=NULL, eps=NULL,
               res=list(estimate=NULL, thres=NULL, iter=NA, tol=NA, time=NULL),
               RESNAME=NULL,
               estimate=res$est, TIME=res$time, thres=res$lambda,
               true, trueeigval = eigen(true)$values,
               trueeigvec = eigen(true)$vectors) {
  ## estimate : list object with covariance matrix estimator / computation time

  norm.frob = norm(as.matrix(estimate - true), type="F")
  norm.spect = max(eigen(as.matrix(estimate - true))$values) ;
  norm.oper = norm(as.matrix(estimate - true), type="O")

  esti.eigen = eigen(estimate) ;
  esti.eigval = as.matrix(t(esti.eigen$values)) ;
  colnames(esti.eigval) = paste("estEig", 1:length(esti.eigval), sep="") ;
  maxeigval = max(esti.eigen$values) - max(trueeigval) ;
  abscos.1stPCangle = as.numeric(abs(esti.eigen$vectors[ ,1] %*% trueeigvec[ ,1])) ;
  
  if (!is.null(GENINFO)) GENINFO = t(GENINFO) ;
  
  return( data.frame(GENINFO, eps, RESNAME, TIME, thres = res$thres, 
            iter = res$iter, tol = res$tol, 
            t(ftn.PR(estimate, true)), 
            frob = norm.frob, spect = norm.spect, oper = norm.oper,
            maxeigval = maxeigval, firstPC = abscos.1stPCangle,
            esti.eigval) ) ;
}

ftn.storeN = function(GENINFO=NULL, eps=NA,
               res=list(estimate=NA, lambda.opt=NA, nz.opt=NA,
                        lambdaseq=NA, errorseq=NA, nzseq=NA,
                        time=NA, iter=NA, tol=NA),
               RESNAME=NULL, DID = NA,
               estimate=res$estimate, TIME=res$time, thres=res$lambda.opt,
               nz = res$nz.opt, iter=res$iter, tol=res$tol,
               true, trueeigval = eigen(true)$values,
               trueeigvec = eigen(true)$vectors) {
  ## estimate : list object with covariance matrix estimator / computation time

  p = ncol(estimate) ;

  norm.frob = norm(as.matrix(estimate - true), type="F")
  norm.spect = max(abs(eigen(as.matrix(estimate - true), only.values=TRUE)$values)) ;
  norm.oper = norm(as.matrix(estimate - true), type="O")

  # diag elements
  t.diag = diag(true) ; e.diag = diag(estimate) ;
  diag.frob = sqrt(sum( (t.diag - e.diag)^2 )) ;
  diag.spect = max(abs( t.diag - e.diag )) ;
  diag.oper = max(abs( t.diag - e.diag )) ;

  # off-diag elements (support)
  ind.off.supp = matrix(as.logical(true), p, p) ;
  ind.off.nsupp = !ind.off.supp ;
  diag(ind.off.supp) = FALSE ;
  t.off.supp = true ; t.off.supp[!ind.off.supp] = 0 ;
  e.off.supp = estimate ; e.off.supp[!ind.off.supp] = 0 ;
  off.supp.frob = norm(as.matrix(e.off.supp - t.off.supp), type="F")
  off.supp.spect = max(abs(eigen(as.matrix(e.off.supp - t.off.supp), only.values=TRUE)$values)) ;
  off.supp.oper = norm(as.matrix(e.off.supp - t.off.supp), type="O")

  # off-diag elements (no-support)
  t.off.nsupp = true ; t.off.nsupp[!ind.off.nsupp] = 0 ;
  e.off.nsupp = estimate ; e.off.nsupp[!ind.off.nsupp] = 0 ;
  off.nsupp.frob = norm(as.matrix(e.off.nsupp - t.off.nsupp), type="F")
  off.nsupp.spect = max(abs(eigen(as.matrix(e.off.nsupp - t.off.nsupp), only.values=TRUE)$values)) ;
  off.nsupp.oper = norm(as.matrix(e.off.nsupp - t.off.nsupp), type="O")  

  esti.eigen = eigen(estimate) ;
  esti.eigval = as.matrix(t(esti.eigen$values)) ;
  esti.eigval.raw = as.numeric(esti.eigval) ;

  mu.S = (max(esti.eigval.raw) + min(esti.eigval.raw)) / 2 ;
  mu.F = mean(esti.eigval.raw) + 
            var(esti.eigval.raw) * (p-1) / p /
            (mean(esti.eigval.raw) - min(esti.eigval.raw)) ;
  colnames(esti.eigval) = paste("estEig", 1:length(esti.eigval), sep="") ;
  maxeigval = max(esti.eigen$values) - max(trueeigval) ;
  abscos.1stPCangle = as.numeric(abs(esti.eigen$vectors[ ,1] %*% trueeigvec[ ,1])) ;
  
  #if (!is.null(GENINFO)) GENINFO = t(GENINFO) ;


  return( data.frame(t(GENINFO), eps, DID, RESNAME, TIME, tuning = thres, 
            iter = iter, tol = tol, nz = nz,
            t(ftn.PR(estimate, true)), 
            frob = norm.frob, spect = norm.spect, oper = norm.oper,
            diag.frob = diag.frob, diag.spect = diag.spect,
            off.supp.frob = off.supp.frob, off.supp.spect = off.supp.spect,
            off.supp.oper = off.supp.oper,
            off.nsupp.frob = off.nsupp.frob, off.nsupp.spect = off.nsupp.spect,
            off.nsupp.oper = off.nsupp.oper,
            mu.S, mu.F, 
            maxeigval = maxeigval, firstPC = abscos.1stPCangle,
            esti.eigval) ) ;
}



## 160428
ftn.storeNN = function(GENINFO=NULL, eps=NA,
               res=list(estimate=NA, lambda.opt=NA, nz.opt=NA,
                        lambdaseq=NA, errorseq=NA, nzseq=NA,
                        time=NA, iter=NA, tol=NA),
               RESNAME=NULL, DID = NA,
               estimate=res$estimate, TIME=res$time, thres=res$lambda.opt,
               nz = res$nz.opt, iter=res$iter, tol=res$tol,
               true, trueeigval = eigen(true)$values,
               trueeigvec = eigen(true)$vectors) {
  ## estimate : list object with covariance matrix estimator / computation time

  p = ncol(estimate) ;

  norm.frob = norm(as.matrix(estimate - true), type="F")
  norm.spect = max(abs(eigen(as.matrix(estimate - true), only.values=TRUE)$values)) ;
  norm.oper = norm(as.matrix(estimate - true), type="O")

  # diag elements
  t.diag = diag(true) ; e.diag = diag(estimate) ;
  diag.frob = sqrt(sum( (t.diag - e.diag)^2 )) ;
  diag.spect = max(abs( t.diag - e.diag )) ;
  diag.oper = max(abs( t.diag - e.diag )) ;

  # off-diag elements (support)
  ind.off.supp = matrix(as.logical(true), p, p) ;
  ind.off.nsupp = !ind.off.supp ;
  diag(ind.off.supp) = FALSE ;
  t.off.supp = true ; t.off.supp[!ind.off.supp] = 0 ;
  e.off.supp = estimate ; e.off.supp[!ind.off.supp] = 0 ;
  off.supp.frob = norm(as.matrix(e.off.supp - t.off.supp), type="F")
  off.supp.spect = max(abs(eigen(as.matrix(e.off.supp - t.off.supp), only.values=TRUE)$values)) ;
  off.supp.oper = norm(as.matrix(e.off.supp - t.off.supp), type="O")

  # off-diag elements (no-support)
  t.off.nsupp = true ; t.off.nsupp[!ind.off.nsupp] = 0 ;
  e.off.nsupp = estimate ; e.off.nsupp[!ind.off.nsupp] = 0 ;
  off.nsupp.frob = norm(as.matrix(e.off.nsupp - t.off.nsupp), type="F")
  off.nsupp.spect = max(abs(eigen(as.matrix(e.off.nsupp - t.off.nsupp), only.values=TRUE)$values)) ;
  off.nsupp.oper = norm(as.matrix(e.off.nsupp - t.off.nsupp), type="O")  

  esti.eigen = eigen(estimate) ;
  esti.eigval = as.matrix(t(esti.eigen$values)) ;
  esti.eigval.raw = as.numeric(esti.eigval) ;

  mu.S = (max(esti.eigval.raw) + min(esti.eigval.raw)) / 2 ;
  mu.F = mean(esti.eigval.raw) + 
            var(esti.eigval.raw) * (p-1) / p /
            (mean(esti.eigval.raw) - min(esti.eigval.raw)) ;
  colnames(esti.eigval) = paste("estEig", 1:length(esti.eigval), sep="") ;
  maxeigval = max(esti.eigen$values) - max(trueeigval) ;
  abscos.1stPCangle = as.numeric(abs(esti.eigen$vectors[ ,1] %*% trueeigvec[ ,1])) ;
  
  #if (!is.null(GENINFO)) GENINFO = t(GENINFO) ;
  esti.eigval.raw.neg = esti.eigval.raw ;
  esti.eigval.raw.neg[ esti.eigval.raw.neg > 0 ] = 0 ;
  bound.frob = sqrt( mean( esti.eigval.raw.neg^2  ) ) ;

  return( data.frame(t(GENINFO), eps, DID, RESNAME, TIME, tuning = thres, 
            iter = iter, tol = tol, nz = nz,
            t(ftn.PR(estimate, true)), 
            frob = norm.frob, spect = norm.spect, oper = norm.oper,
            diag.frob = diag.frob, diag.spect = diag.spect,
            off.supp.frob = off.supp.frob, off.supp.spect = off.supp.spect,
            off.supp.oper = off.supp.oper,
            off.nsupp.frob = off.nsupp.frob, off.nsupp.spect = off.nsupp.spect,
            off.nsupp.oper = off.nsupp.oper,
            mu.S, mu.F, 
            maxeigval = maxeigval, firstPC = abscos.1stPCangle,
            esti.eigval, bound.frob) ) ;
}




ftn.eval.list = function(estimatelist, true, trueeigval = eigen(true)$values,
                 trueeigvec = eigen(true)$vectors) {

  num.esti = length(estimatelist) ;
  resmat = NULL ;
  for (i in 1:num.esti) {
    res = ftn.eval(estimate = estimatelist[[i]], true, trueeigval, trueeigvec) ;
    resmat = rbind(resmat, res) ;
  }
  rownames(resmat) = names(estimatelist) ;

  return(resmat) ;
                 
}