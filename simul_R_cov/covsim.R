## covsim.R
## replicates the results in simulation studies (empirical estimation errors).

## set working directory as that stores covsim.R
setwd("/Users/cyg/Dropbox/Codes/FSPD_public/simul_R_cov/")

## load libraries
library(matrixcalc) ;
library(Matrix) ;
library(PDSCE) ;
library(flare) ; 
library(gtools) ; 

## load auxilary codes (contains functions)
source("../src_R/datagenerate.R") # contains gen.datamat(SAMPSIZE, DIM, DISTSTR)
source("../src_R/functions.R") # dependent on datagenerate.R
source("../src_R/evaluations.R")
source("../src_R/cov_estimators.R") ;
source("../src_R/prec_estimators.R") ;

## how many replications will be done? (default 100)
REPLNUM = 1 ;

## choose epsilon (cut-off point to do PD-modification)
eps = 10^-2 ;

## sample size and data dimension
n = 100 ;
p = 100 ;

## generate true covariance matrices
source("../src_R/truecovmats.r") ;
  
## true covariance matrix : taper (model 1), OvrlpBdiag (model 2)
matstrseq = c("taper", "OvrlpBdiag") ;  
## data generating distribution : MVN, MV-t
diststrseq = c("mvn", "mvt") ;

## set mu-value for linear shrinkage
##  : FSopt ("FSPD(\mu_{SF})"), Infty ("FSPD(\infty)") 
LSseq = c("FSopt", "Infty") ;
  
## empirical error is measured for each combination of (true cov, distribution)
for (MATSTR in matstrseq) { for (DISTSTR in diststrseq) {
  
  ## declare true covarinace matrix  
  eval(parse(text = paste("cov.true = mat.", MATSTR, sep=""))) ;
  
  ## calculate eigenvalue decomposions from the true covariance matrix
  eig.cov.true = eigen(cov.true) ;
  eigvec.cov.true = eig.cov.true$vectors ;
  eigval.cov.true = eig.cov.true$values ;
  sqrt.cov.true = eigvec.cov.true %*% diag(sqrt(eigval.cov.true)) %*% t(eigvec.cov.true) ;
  maxeigval.cov.true = max(eigval.cov.true) ;
  mineigval.cov.true = min(eigval.cov.true) ;

  ## matrix storing results
  resmat = NULL ;

  ## start replication
  for (repl in 1:REPLNUM) {
    
    ## print current replication status
    geninfo = c(repl, MATSTR, DISTSTR, n, p) ;
    names(geninfo) = c("repl", "struct", "dist", "n", "p")
    print(geninfo) ;
    
    ## generate data (X)
    set.seed(repl*10) ;
    data.gen = gen.datamat(SAMPSIZE = n, DIM = p, DISTSTR = DISTSTR) ;
    data.gen = data.gen %*% sqrt.cov.true ;
    ## centering data
    data = apply(data.gen, 2, scale, center=TRUE, scale=FALSE) ;
    ## sample covariance matrix
    cov.samp = t(data) %*% data / (n-1) ;
    
    ## calculate soft thresholding estimator and its information
    ## name: "uni"
    ## function "ftn.cov.uni" calculates the estimator
    uni.list = ftn.cov.uni(datamat = data) ;
    uni = uni.list$estimate ;
    uni.lam = uni.list$lambda.opt ;
    ## function "ftn.eig" calculates the largest/smallest eigenvalues of the estimator
    uni.mineigval = ftn.eig(uni)[2] ;
    ## function "ftn.storeN" extracts information of the estimator
    res.uni = ftn.storeN(GENINFO=geninfo, eps=NA,
                         res = uni.list,
                         RESNAME="uni",
                         true=cov.true, trueeigval = eigval.cov.true, trueeigvec = eigvec.cov.true)
    
    ## calculate FSPD update of soft thresholding estimator
    ## name: "uniL"
    ## function "subftn.LS" is our main procedure, updating first-stage estimator for PDness
    res.uniL = NULL ;
    for (LSoption in LSseq) {
      time1 = Sys.time() ;
      if (eps >= uni.mineigval) {
        uniL = subftn.LS(cov.raw = uni, mineigval = eps, method=LSoption) ; 
        DID = 1 ;
      } else {
        uniL = uni ;
        DID = 0 ;
      }
      time2 = difftime(Sys.time(), time1, units="secs") ;

      res.uniL = rbind(res.uniL, 
        ftn.storeN(GENINFO=geninfo, eps = eps, DID = DID,
          RESNAME= paste("uniL", LSoption, sep=""),
          estimate=uniL, nz = ftn.sparsity(uniL), TIME=time2, thres=NA,
          true=cov.true, trueeigval = eigval.cov.true, trueeigvec = eigvec.cov.true)
        ) ;
    }
    
    ## calculate Xue's eigenvalue constraint estimator
    ## name: "uniX"  
    ## function "subftn.Xue" calculates the estimator
    ## set tolerence level for algorithm convergence
    uniX.tol = 1e-7 ;
    time1 = Sys.time() ;
    if (eps >= uni.mineigval) {
      uniX = subftn.Xue(sampcov = cov.samp, initial.cov.thres = uni, 
                        cut = uni.lam, mineigval=eps, tolprimal=uniX.tol,
                        toldual = 1e-5, MAXITER=1000) ;
      uniX.iter = ifelse(is.null(attr(uniX, 'iter')), NA, attr(uniX, 'iter')) ;
      DID = 1 ;
    } else {
      uniX = uni ;
      uniX.iter = NA ;
      DID = 0 ;
    }
    time2 = difftime(Sys.time(), time1, units="secs") ;
    res.uniX = ftn.storeN(GENINFO=geninfo, eps = eps, DID = DID,
      RESNAME="uniX",
      estimate=uniX, nz = ftn.sparsity(uniX), TIME=time2, thres=NA,
      iter=uniX.iter, tol=uniX.tol,
      true=cov.true, trueeigval = eigval.cov.true, trueeigvec = eigvec.cov.true) ;
    
    ## calculate Rothman's log-det barrier estimator
    ## name: "uniR"  
    ## function "pdsoft" (depends on package PDSCE) calculates the estimator
    ## set tolerence level for algorithm convergence    
    uniR.tol = 1e-7 ;
    ## set tau for the barrier level
    tau = 1e-2 ;
    time1 = Sys.time() ;
    if (eps >= uni.mineigval) {
      uniR.message <- capture.output(uniR.list <- pdsoft(s = cov.samp, lam = uni.lam, tau = tau,
                                                         init="soft", standard=FALSE, tolout=uniR.tol, maxitout=1000, quiet=FALSE))
      time2 = difftime(Sys.time(), time1, units="secs") ;
      uniR.iter = as.numeric(substr(uniR.message, start=25, stop=nchar(uniR.message)))
      uniR = uniR.list$sigma ;
    } else {
      uniR = uni ;
      uniR.iter = NA ;
    }
    res.uniR = ftn.storeN(GENINFO=geninfo, eps = eps, DID = DID,
      RESNAME="uniR",
      estimate=uniR, nz = ftn.sparsity(uniR), TIME=time2, thres=NA,
      iter=uniR.iter, tol=uniR.tol,
      true=cov.true, trueeigval = eigval.cov.true, trueeigvec = eigvec.cov.true) ;
      
    
    ## aggregate results
    resmat.uni = rbind(res.uni, res.uniL, res.uniX, res.uniR) ;
    resmat = rbind(resmat, resmat.uni) ;
    
    
    ## save the results
    storename = paste("cov_", MATSTR, DISTSTR, "size",n,"dim",p, sep="") ;
    
    write.table(resmat, paste(storename, ".txt", sep=""), col.names=TRUE, row.names=FALSE, 
                quote=FALSE, sep="\t") ;
    save.image(paste(storename,".RData", sep=""))
    
    print(Sys.time())
    ## REPL END
  }
  
    
## structure / distribution ends
}}
