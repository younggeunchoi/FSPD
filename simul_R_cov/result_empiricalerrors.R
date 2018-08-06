######## RESULT PRINT : empirical errors


library(Matrix) ; 
library(flare) ;

## set working directory as that stores covsim.R
setwd("/Users/cyg/Dropbox/Codes/LSPD_pub/simul_cov/")

errname = c("oper", "spect", "frob") ;
#"maxeigval", "firstPC",
#"FPR", "TPR") ;
prt.errname = c("Matrix $l_1$", "Spect.", "Frob.") ; #"MaxEigval", "1stPCCosAngle",
# "FPR", "TPR") ;

n = 100 ; 
prt.matseq = c("Tapered", "Overlap. block diag.")
matseq = c("taper", "OvrlpBdiag");

  
for (k in 1:length(matseq)) {
  #p = 100 ;
  MATSTR = matseq[k]
  prt.mat = prt.matseq[k] ;
  cat("\\hline \n") ;
  cat("(", prt.mat,") ",sprintf("& %s ", prt.errname), sprintf("& %s ", prt.errname), "\\\\\n", sep="") ;
  cat("\\hline \n") ;
  for (p in c(100)) {
    
    cat(sprintf("\\multicolumn{7}{c}{$p = %d$}", p), "\\\\\n" ) ;  
    #### contains gen.datamat(SAMPSIZE, DIM, DISTSTR)
    source("../src_R/datagenerate.R")
    source("../src_R/functions.R") # dependent on datagenerate.R
    source("../src_R/evaluations.R")
    source("../src_R/cov_estimators.R") ;
    source("../src_R/prec_estimators.R") ;
    source("../src_R/functions.R") ;
    source("../src_R/truecovmats.R") ;
    source("../src_R/trueprecmats.R") ;
    
    eps = 10^c(-2) ;
    eval(parse(text = paste("mat.true = mat.", MATSTR, sep=""))) ;
    tol = 10^-7 ;
    
    
    # names(resmat.avg)
    # [1] "struct"      "dist"        "n"           "p"           "eps"         "DID"         "RESNAME"    
    # [8] "TIME"        "tuning"      "iter"        "tol"         "nz"          "FPR"         "TPR"        
    # [15] "frob"        "spect"       "oper"        "mu.S"        "mu.F"        "maxeigval"   "mineig"     
    # [22] "negeig.prop" "negeig.case"
    
    
    resnameseq = c("uni", "uniLFSopt", "uniLInfty", "uniX", "uniR") ;
    prt.resnameseq = c("Soft thresholding", "+ Linear shrinkage", "+ Diagonal shift",
                       "+ Eigen. constraint", "+ log-det barrier") ;
    
    prt.matname = prt.matseq[1] ;
    
    
    for (j in 1:length(resnameseq)) { 
      
      resname = resnameseq[j] ;
      prt.resname = prt.resnameseq[j] ;
      cat(prt.resname) ;        
      
      for (DISTSTR in c("mvn", "mvt")) {    
        common.name =  paste("cov_", MATSTR, DISTSTR, "size",n,"dim", p, sep="") ;
        resmat.avg = read.csv(paste(common.name,"summ_avg.csv", sep="")) ;
        resmat.se = read.csv(paste(common.name,"summ_se.csv", sep="")) ;
        
        #            if( (eps == eps2) & ((j == 5) | (j == 5) ) ) {
        #              cat(" & & & ") ;
        #              next ;
        #            }
        #           
        
        avg.temp = resmat.avg[ (resmat.avg$RESNAME == resname),  ]
        se.temp = resmat.se[ (resmat.se$RESNAME == resname),  ]
        
        
        if ( nrow(avg.temp) != 1 ) {
          if ( avg.temp$RESNAME[1] != "uniR" ) {
            avg.temp = avg.temp[ abs(avg.temp$eps - eps) < 10^-6 , ] ;
            se.temp = se.temp[ abs(avg.temp$eps - eps) < 10^-6 , ] ;
          } else {
            avg.temp = avg.temp[ (avg.temp$eps == 0.01), ] ;
            se.temp = se.temp[ (se.temp$eps == 0.01), ] ;
          }   
        }
        
        if ( nrow(avg.temp) != 1 ) {
          avg.temp = avg.temp[ (avg.temp$tol == tol), ] ;
          se.temp = se.temp[ (avg.temp$tol == tol), ] ;
        }
        
        if ( avg.temp$RESNAME[1] == "uniR" ) {
          avg.temp = avg.temp[1, ] ;
          se.temp = se.temp[1, ] ;
        }
        
        avg.temp = avg.temp[1, errname] ; se.temp = se.temp[1, errname] ;
        cat(sprintf(" & %.2f (%.2f)", avg.temp, se.temp)) ;
        
        
      }
      cat(" \\\\\n") ;
      
    }
    # p end
  }
  # matstr end
}
# diststr end


