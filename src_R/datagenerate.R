

gen.datamat = function(SAMPSIZE, DIM, DISTSTR) {
  ## INPUTS : SAMPSIZE means n, DIM means p
  ## (character) DIST = "gamma", "mvn", "mvt"
  ## (character) COV = "taper", "sparse"

  n = SAMPSIZE ;
  p = DIM ;
  
  if (DISTSTR=="mvn") data = matrix(rnorm(n*p), nrow=n, ncol=p) ;
  if (DISTSTR=="mvt") {
    y = matrix(rnorm(n*p), nrow=n, ncol=p) ;
    u = rchisq(n, df=5) ;
    data = sqrt(3/5) * y / sqrt(u/5) ;
  }
  if (DISTSTR=="gamma") {
    data = matrix(rgamma(n*p, shape=4, scale=0.5), nrow=n, ncol=p) ;
  }
  return(data) ;
}