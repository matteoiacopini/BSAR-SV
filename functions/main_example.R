# Example code to estimate the Bayesian SAR-SV model in
# "Bayesian SAR model with Stochastic Volatility and time-varying weights"
# by Costola, M., Iacopini, M., and Wichers, C. (2024)
# Journal of Financial Econometrics


#### Load packages and functions ####
cat("\014")
rm(list=ls())

# setwd("C:/...")

# packages for MCMC
library(MASS)
library(DirichletReg)
library(stochvol)
library(Rcpp)
library(dplyr)
library(mvtnorm)
library(matrixcalc)
library(gsubfn)
library(invgamma)    # for "rinvgamma"
library(wordspace)   # for "normalize.rows"
library(beepr)

# functions needed for MCMC
source("uni.slice.R")
source("sample_rho_sv.R")
source("log_target_dens_delta_sv.R")
source("log_proposal_delta.R")
source("metropolis_delta_sv.R")
source("norm_rows.R")
source("max_row_norm.R")
source("gibbs_sv.R")
source("sample_eta.R")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### Load the data ####

load("ExampleData.RData")




#### MCMC algorithm ####

# setup
nsave = 2000     # no. iterations to store
nburn = 1500     # no. burn-in iterations (not stored)
nshow = 100      # show time every nshow iters

# choose the hyperparameters
nu = rep(2,2)
rho_params = rep(0.5,2)    # hyperparams of prior Be(a,b)
prior_rho  = "beta_-11"    # Choose:  "beta_01" --> rho in (0,1)  or "beta_-11" --> rho in (-1,1)

# estimation
res = gibbs_sv(nsave,nburn,nshow, ft= data_sarsv$ft, nu, rho_params, y= data_sarsv$yt,
                        min_network= data_sarsv$net1_norm, plus_network= data_sarsv$net2_norm,
                        prior_rho, fixed.effects=TRUE)

beep(1, expr = cat("%---%---%---%---%---%---%\n%   MCMC run finished   %\n%---%---%---%---%---%---%"))



######## Direct, Indirect effects ########
nn   = dim(data_sarsv$yt)[2] 
In   = diag(nn)
iota = matrix(1, nn,1)
direct.store   = array(NA, c(nsave, T,nn))
indirect.store = array(NA, c(nsave, T,nn))
avg.direct.store   = array(NA, c(nsave, T))
avg.indirect.store = array(NA, c(nsave, T))
S0dy.store = array(NA, c(nsave, T))
for (mm in 1:nsave){
   for (t in 1:T){
      Wt.star = (delta.store[1]*data_sarsv$net1_norm[[t]] + delta.store[2]*data_sarsv$net2_norm[[t]])
      rrr = rho.store[mm,]
      rrr[rrr==1] = 0.99
      Vt = solve(In - diag(rrr) %*% Wt.star)
      Vtd = diag(diag(Vt))
      S0dy.store[mm,t] = sum(Vt^2-diag(diag(Vt^2))) / sum(diag(tcrossprod(Vt)))
      direct.store[mm,t,]   = crossprod(iota, Vtd)
      indirect.store[mm,t,] = crossprod(iota, Vt-Vtd)
      avg.direct.store[mm,t]   = (crossprod(iota, Vtd) %*% iota) / nn
      avg.indirect.store[mm,t] = (crossprod(iota, Vt-Vtd) %*% iota) / nn
   }
}
