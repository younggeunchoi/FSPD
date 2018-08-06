
###############################################################################
###############################################################################
###############################################################################
## COMMON PARTS FOR BOTH PLOT AND RESULT MATRIX
library(Matrix) ; 
library(flare) ;



## set working directory as that stores covsim.R
setwd("/Users/cyg/Dropbox/Codes/LSPD_pub/simul_cov/")

avgstack = NULL;
sestack = NULL ;

n = 100 ;
#p = 100
p = 100 ;
  

#### contains gen.datamat(SAMPSIZE, DIM, DISTSTR)
source("../src_R/datagenerate.R")
source("../src_R/functions.R") # dependent on datagenerate.R
source("../src_R/evaluations.R")
source("../src_R/cov_estimators.R") ;
source("../src_R/prec_estimators.R") ;
source("../src_R/functions.R") ;
source("../src_R/truecovmats.R") ;
source("../src_R/trueprecmats.R") ;


# MATSTR = "taper" ; DISTSTR = "mvn" ;
for (MATSTR in c("taper", "OvrlpBdiag")) { for (DISTSTR in c("mvn", "mvt")) {
  common.name =  paste("cov_", MATSTR, DISTSTR, "size",n,"dim", p, sep="") ;
  
  
  eval(parse(text = paste("cov.true = mat.", MATSTR, sep=""))) ;
  eig.cov.true = eigen(cov.true) ;
  eigvec.cov.true = eig.cov.true$vectors ;
  eigval.cov.true = eig.cov.true$values ;
  
  # collect all uni's
  resmat = read.table(paste(common.name, ".txt", sep=""), header=T, sep="\t") ;
  
  resmat = resmat[ order(resmat$repl, resmat$RESNAME, resmat$eps, resmat$tol) , ]
  resmat[is.na(resmat)] = -1 ;
  
  
  matseq = levels(resmat$struct) ; #covseq =  ;
  distseq = levels(resmat$dist) ; #distseq = distseq[c(2,3,1)] ;
  resnameseq = levels(resmat$RESNAME) ;
  epsseq = 10^c(-2) ;
  
  #covind = "taper" ; distind = "mvn" ; estind="uni.thres" ; epsind = 0 ;
  #precind = "band" ; distind = "mvn" ; estind="concord" ; epsind = 0 ;
  
  ###############################################################################
  ###############################################################################
  ###############################################################################
  
  
  
  
  ######## ERROR INSPECTION :: MEAN, SD of replications
  
  # > names(resmat)[1:30] ;
  # [1] "repl"      "struct"    "dist"      "n"         "p"         "eps"       "DID"       "RESNAME"   "TIME"     
  # [10] "tuning"    "iter"      "tol"       "nz"        "FPR"       "TPR"       "frob"      "spect"     "oper"     
  # [19] "mu.S"      "mu.F"      "maxeigval" "firstPC"   "estEig1"   "estEig2"   "estEig3"   "estEig4"   "estEig5"  
  # [28] "estEig6"   "estEig7"   "estEig8"  
  
  eiginfo = function(eigseq) {
    negprop = sum(eigseq < 0) / length(eigseq) ;
    smallesteig = min(eigseq) ;
    negcase = (smallesteig < 0) ;
    smallcase = (smallesteig < 0.01) ;
    
    eig.diff = as.numeric(abs(eigval.cov.true - eigseq)) ;
    max.eig.diff = max(eig.diff) ;
    last.eig.diff = eig.diff[p] ;
    meansq.eig.diff = sqrt(mean(eig.diff^2)) ;
    
    min.eig = min(eigseq) ;
    bound1 = max(c(0.01 - min.eig, 0)) * var(eigseq) * (length(eigseq) + 1) / sum((eigseq - min.eig)^2) ;
    bound2 = max(c(0.01 - min.eig, 0)) ;
    temp = 0.01 - eigseq ; temp[temp < 0] = 0 ;
    bound3 = sqrt( mean(temp^2) ) ;
    return(c(mineig = smallesteig, negeig.prop = negprop, negeig.case = negcase, smalleig.case=smallcase,
             max.eig.diff = max.eig.diff, last.eig.diff = last.eig.diff, meansq.eig.diff = meansq.eig.diff,
             Frob.lowerbound = bound1, Frob.upperbound= bound2, Frob.optbound = bound3)) ;
  }
  
  eig.summ = apply(resmat[ ,grep("estEig", names(resmat))], 1 ,eiginfo) ;
  
  resmat.error = cbind(resmat[ ,1:29], t(eig.summ)) ;
  # > names(resmat.error)
  # [1] "repl"        "struct"      "dist"        "n"           "p"           "eps"         "DID"        
  # [8] "RESNAME"     "TIME"        "tuning"      "iter"        "tol"         "nz"          "FPR"        
  # [15] "TPR"         "frob"        "spect"       "oper"        "mu.S"        "mu.F"        "maxeigval"  
  # [22] "mineig"      "negeig.prop" "negeig.case" 
  
  avgtarget =  resmat.error[ ,c(1, 9:ncol(resmat.error))] ;
  
  description = resmat.error[resmat.error[ ,1]==1, c(1:8)] ;
  stackmat.sum = matrix(0, nrow=nrow(avgtarget)/100, ncol=ncol(avgtarget)) ;
  stackmat.sqsum = matrix(0, nrow=nrow(avgtarget)/100, ncol=ncol(avgtarget)) ;
  for (replnum in 1:100) {
    temp = avgtarget[avgtarget$repl == replnum, ] ;
    temp$repl = as.numeric(temp$repl) ;
    stackmat.sum = stackmat.sum + temp ;
    stackmat.sqsum = stackmat.sqsum + temp^2 ;  
  }
  stackmat.avg = stackmat.sum / 100 ;
  temp2 = stackmat.sqsum / 100 - stackmat.avg^2 ;
  temp2[abs(temp2) < 10^-8] = 0 ;
  stackmat.se = sqrt(temp2) / sqrt(100) ;
  
  print(head(stackmat.avg), digits=3) ;
  
  
  resmat.avg = cbind(description[ ,-1], stackmat.avg[ ,-1]) ;
  resmat.se = cbind(description[ ,-1], stackmat.se[ ,-1]) ;
  avgstack = rbind(avgstack, resmat.avg) ;
  sestack = rbind(sestack, resmat.se) ;
  
  common.name2 = paste("cov_", MATSTR, DISTSTR, "size",n,"dim", p, sep="") ;    
  write.csv(resmat.avg, paste(common.name2,"summ_avg.csv", sep=""), row.names=FALSE) ;
  write.csv(resmat.se, paste(common.name2,"summ_se.csv", sep=""), row.names=FALSE) ;
  # 
  # save.image("result_organize_150120_p100.RData")
  
  ## diststr, matstr ends
}}

write.csv(avgstack, paste("cov_","summ_avg.csv", sep=""), row.names=FALSE) ;
write.csv(sestack, paste("cov_","summ_se.csv", sep=""), row.names=FALSE) ;











