# require functions.R
# require package 'flare'



set.seed(50) ;
mat.unif = matrix(runif(n = p^2, min= 0.5, max=0.5) ,nrow = p) ;
set.seed(40) ;
#mat.sign = rbinom(n = p^2, size = 1, prob = 0.5) * 2 - 1 ;
#mat.unif = mat.unif * mat.sign ;

## taper
#mat.taper = cov.taper ;
## block diagonal
#mat.OvrlpBdiag = cov.OvrlpBdiag ;
#
#mat.AR = cov.AR ;



### scale-free but random magnititude
#?sugm.generator
list.scfr = sugm.generator(n = n, d = p, graph="scale-free",
                           v = 0.5, u = 0, seed = 50, vis=FALSE, verbose=FALSE)
nz.scfr = list.scfr$theta ;
#ftn.eig(mat.scfr)
mat.sparse.unif = mat.unif ;
mat.sparse.unif[!matrix(as.logical(nz.scfr), p,p)] = 0 ;
mat.sparse.unif[upper.tri(mat.sparse.unif)] = 0 ;
mat.sparse.unif = mat.sparse.unif + t(mat.sparse.unif) ;
diag(mat.sparse.unif) = 1 ;
mat.scfr = subftn.DS(mat.sparse.unif, eig.min=0.00001) ;
rm(list.scfr, nz.scfr) ;
#mat.scfr2 = scalemat %*% mat.scfr %*% scalemat ;


### Erdos-Renyi(random) but random magnititude with 4% sparsity
list.ER = sugm.generator(n = n, d = p, graph="random", v = 0.5, u = 0,
                         prob= 0.15/(p/100), seed = 704, vis=FALSE, verbose=FALSE)
mat.ER = list.ER$omega ;
nz.ER = list.ER$theta ;
mat.sparse.unif = mat.unif ;
mat.sparse.unif[!matrix(as.logical(nz.ER), p,p)] = 0 ;
mat.sparse.unif[upper.tri(mat.sparse.unif)] = 0 ;
mat.sparse.unif = mat.sparse.unif + t(mat.sparse.unif) ;
diag(mat.sparse.unif) = 1 ;
# diagonal shift
#mat.ER = subftn.LS2(mat.sparse.unif, eig.min=0.00001) ;
mat.ER1 = subftn.DS(mat.sparse.unif, eig.min=1e-1) ;
mat.ER3 = subftn.DS(mat.sparse.unif, eig.min=1e-3) ;
mat.ER5 = subftn.DS(mat.sparse.unif, eig.min=1e-5) ;
#ftn.eig(mat.ER)
rm(list.ER, nz.ER) ;


list.ERD = sugm.generator(n = n, d = p, graph="random", v = 0.5, u = 0,
                         prob= 1, seed = 704, vis=FALSE, verbose=FALSE)
mat.ERD = list.ERD$omega ;
nz.ERD = list.ERD$theta ;
mat.sparse.unif = mat.unif ;
mat.sparse.unif[!matrix(as.logical(nz.ERD), p,p)] = 0 ;
mat.sparse.unif[upper.tri(mat.sparse.unif)] = 0 ;
mat.sparse.unif = mat.sparse.unif + t(mat.sparse.unif) ;
diag(mat.sparse.unif) = 1 ;
# diagonal shift
# #mat.ERD = subftn.LS2(mat.sparse.unif, eig.min=0.00001) ;
# mat.ERD1 = subftn.DS(mat.sparse.unif, eig.min=1e-1) ;
# mat.ERD3 = subftn.DS(mat.sparse.unif, eig.min=1e-3) ;
# mat.ERD5 = subftn.DS(mat.sparse.unif, eig.min=1e-5) ;
# #ftn.eig(mat.ERD)
mat.ERD = mat.sparse.unif ;
rm(list.ERD, nz.ERD) ;

### band but random magnititude with 5-nbd
list.band = sugm.generator(n = n, d = p, graph="band", g = 10, v = 0.75, u = 0,
                         seed = 1, vis=FALSE, verbose=FALSE)
mat.band = list.band$omega ;
nz.band = list.band$theta ;
mat.sparse.unif = mat.unif ;
mat.sparse.unif[!matrix(as.logical(nz.band), p,p)] = 0 ;
mat.sparse.unif[upper.tri(mat.sparse.unif)] = 0 ;
mat.sparse.unif = mat.sparse.unif + t(mat.sparse.unif) ;
diag(mat.sparse.unif) = 1 ;
#ftn.eig(mat.sparse.unif)
#mat.band = subftn.LS2(mat.sparse.unif, eig.min=0.00001) ;
mat.band1 = subftn.DS(mat.sparse.unif, eig.min=1e-1) ;
mat.band3 = subftn.DS(mat.sparse.unif, eig.min=1e-3) ;
mat.band5 = subftn.DS(mat.sparse.unif, eig.min=1e-5) ;
#ftn.eig(mat.band)
rm(list.band, nz.band) ;


#mat.TD = triDiag(diagonal=1, upper=0.5, lower=0.5, nrow=p) ;
##mat.TD = subftn.LS2(mat.TD, eig.min=0.00001) ;
#ftn.eig(mat.TD)



### hub but random magnititude
list.hub = sugm.generator(n = n, d = p, graph="hub",
                           v = 0.5, u = 0, seed = 1, vis=FALSE, verbose=FALSE)
#mat.hub = list.hub$omega ;
nz.hub = list.hub$theta ;

mat.sparse.unif = mat.unif ;
mat.sparse.unif[!matrix(as.logical(nz.hub), p,p)] = 0 ;
mat.sparse.unif[upper.tri(mat.sparse.unif)] = 0 ;
mat.sparse.unif = mat.sparse.unif + t(mat.sparse.unif) ;
diag(mat.sparse.unif) = 1 ;
#ftn.eig(mat.sparse.unif)
mat.hub = subftn.DS(mat.sparse.unif, eig.min=0.001) ;
#ftn.eig(mat.hub)
rm(list.hub, nz.hub) ;

rm(mat.unif, mat.sparse.unif) ;

#mat.ER = subftn.LS2(mat.ER, eig.max=5, eig.min=0.005)

# v : off-diagonal, u : added to diagonal
# random sparse
# dense
