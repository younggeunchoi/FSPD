######## RESULT PRINT : eigenvalue information
## MIN. EIGVAL (SD), NEG. PROP (SD), NEG. CASE
# names(resmat.avg)
# [1] "struct"        "dist"          "n"             "p"             "eps"           "DID"          
# [7] "RESNAME"       "TIME"          "tuning"        "iter"          "tol"           "nz"           
# [13] "FPR"           "TPR"           "frob"          "spect"         "oper"          "mu.S"         
# [19] "mu.F"          "maxeigval"     "mineig"        "negeig.prop"   "negeig.case"   "smalleig.case"


## set working directory as that stores covsim.R
setwd("/Users/cyg/Dropbox/Codes/LSPD_pub/simul_cov/")

prt.matseq = c("Tapered", "Overlap. block diag.")
matseq = c("taper", "OvrlpBdiag");
distseq = c("mvn", "mvt") ;
prt.distseq = c("$\\mathcal{N}$", "$t$") ; 

REPLNUM = 100 ;
resname = "uni" ;



errname = c("mineig", "negeig.prop", "negeig.case")#, "smalleig.case")
prt.errname = c("Min. eig.", "\\#(Neg. eig.)/p",  "\\#(PD)")#, "# of (Min. Eig. < 0.01)") ;

for (temp in 1:1) {
  
  cat("\\hline\n") ;
  cat("& ", sprintf("& \\multicolumn{2}{c}{%s}", prt.matseq), " \\\\\n", sep="") ;
  cat("$p$ & ", sprintf("& %s ", prt.errname), sprintf("& %s ", prt.errname),# sprintf("& %s", prt.errname),
      " \\\\\n", sep="") ;
  cat("\\hline\n") ;
  
  n = 100 ; 
  #p = 100 ;
  for (p in c(100)) {
    
    for (l in 1:length(distseq)) {
      DISTSTR = distseq[l] ;
      if (DISTSTR == "mvn") cat(sprintf(" %d & %s ", p, prt.distseq[l]))     ;
      if (DISTSTR == "mvt") cat(sprintf("     & %s ", prt.distseq[l]))     ;
      
      for (i in 1:2) {
        MATSTR = matseq[i] ;
        
        common.name =  paste("cov_", MATSTR, DISTSTR, "size",n,"dim", p, sep="") ;
        
        resmat.avg = read.csv(paste(common.name,"summ_avg.csv", sep="")) ;
        resmat.se = read.csv(paste(common.name,"summ_se.csv", sep="")) ;
        
        avg.temp = resmat.avg[ (resmat.avg$RESNAME == resname), errname ]
        se.temp = resmat.se[ (resmat.se$RESNAME == resname), errname  ]
        
        for (j in 1:length(errname)) {
          errname.temp = errname[j] ;
          if (errname.temp == "mineig") {
            cat(sprintf(" & %.3f (%.3f)", avg.temp[errname.temp], se.temp[errname.temp]))     ;
         } 
          if (errname.temp == "negeig.prop") {
            #cat(sprintf(" & %.3f", avg.temp[errname.temp]))     ;
            cat(sprintf(" & %.3f (%.3f)", avg.temp[errname.temp], se.temp[errname.temp]))     ;
          }
          if (errname.temp == "negeig.case" | errname.temp == "smalleig.case") {
            cat(sprintf(" & %d/%d", as.integer(100*(1-avg.temp[errname.temp])), as.integer(REPLNUM))) ;
          }
        }
      }
      cat(" \\\\\n") ;
    }
    cat("\\hline\n") ;
  }
  #dim(resmat.error)
  
}
