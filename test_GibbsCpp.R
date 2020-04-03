source("./R/misc.R")
source("./R/Main_Sampler.R")
source("./R/Simu_data_Sampling.R")
require(Matrix)
require(Rcpp)
require(RcppArmadillo)
library(lineprof)
Rcpp::sourceCpp("src/IsingCpp_CFTP_sparse.cpp")

# test Pdet_Ising_single_siteCpp:
thr <- matrix(rnorm(10),5,2)
Z=matrix(1,2,1)
dethis <- 2*(matrix(runif(10),5,2)<0.5)+1
dethis[1,]=NA

sppmat_det <- matrix(c(0,0.3,0.3,0),2,2)

Pdet_Ising_single_siteCpp(thr,Z,dethis,sppmat_det,c(-1L,1L))
Pdet_Ising_single_site(thr,Z,dethis,sppmat_det)

lineprof(lapply(1:1000,function(dummy){ Pdet_Ising_single_siteCpp(thr,Z,dethis,sppmat_det,c(-1L,1L))}))
lineprof(lapply(1:1000,function(dummy){ Pdet_Ising_single_site(thr,Z,dethis,sppmat_det)}))

# test extract_thrCpp
set.seed(42)
thr_list <- list()
thr_list[[1]] <- matrix(runif(6),3,2)
thr_list[[2]] <- matrix(runif(6),3,2)

extract_thr(1,thr_list)
extract_thrCpp(0,thr_list,2,2,3)

lineprof(lapply(1:10000,function(dummy) extract_thr(1,thr_list)))
lineprof(lapply(1:10000,function(dummy) extract_thrCpp(0,thr_list,2,2,3)))

