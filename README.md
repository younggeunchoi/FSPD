# FSPD
Implementation code for the paper of "Fixed support positive-definite modification of covariance matrix estimators via linear shrinkage" by Choi, Roy, Park and Lim (arXiv: 1606.03814)


* Try `vignette.R` to see how FSPD works.

* The simulation study (Chap 4) can be replicated as follows: run the followings consequtively (1 -> 2 -> 3 and 4): 
        1. (repeatedly generate data and calculate estimators and store results) `/simul_R_cov/simul_cov/covsim.R`
        2. (aggregate results over replications) `simul_R_cov/replication_aggregate.R`
        3. (2.1 empirical spectrum) `/simul_R_cov/result_empiricalspectrum.R`
        4. (2.2 error comparison) `/simul_R_cov/result_empiricalerrors.R`
