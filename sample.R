##############################################################################
## Semiparametric least squares regression with doubly-censored data
## Taehwa Choi, Arlene K. H. Kim and Sangbum Choi
## Update: March 5, 2021
##############################################################################
# rm(list=ls())
library(tidyverse)
library(survival)

## Functions
simudata = function(n, beta0=c(1,1)){
  tau = 20
  x1 = rnorm(n)
  x2 = ifelse(rnorm(n) > 0, 1, 0)
  x = cbind(x1, x2)
  
  err = rnorm(n); l = 2; r = 30
  
  T = exp(x %*% beta0 + err)
  L = (1 - 0.5 * x2) * runif(n, -2, l)
  R = L + (1 - 0.5 * x2) * runif(n, 3, r)
  Y = pmin(pmin(R, tau), pmax(T, L))
  delta = case_when(T < L ~ 3,
                    T > R ~ 2,
                    TRUE ~ 1)
  d = data.frame(Y, delta, x1, x2, T, L, R)
  
  return(d[order(Y), ])
}


km_dc = function(Y, delta, weights = rep(1,n)) {
  t = sort(unique(Y))
  m = length(t)
  lambda = lambda_old = rep(1/m,m)
  at_risk = outer(Y, t, ">=")
  ind_eq = outer(Y, t, "==")
  
  maxiter = 20; tol = 1e-5; error = 10; iter = 0
  while (iter < maxiter & error > tol) {
    #Estep
    S = drop(at_risk %*% lambda)
    temp = 1/(1 - exp(-S))
    EW = outer(I(delta==3)*temp, lambda, "*")
    EW = EW * at_risk
    EW = EW + I(delta==1) * ind_eq
    EW = EW + outer(rep(1,n), lambda, "*") * (1-at_risk)
    
    #Mstep
    num = apply(at_risk * EW * weights, 2, sum)
    denom = apply(at_risk * weights, 2, sum)
    lambda = num / denom
    error = sum((lambda - lambda_old)^2)
    iter = iter + 1
    lambda_old = lambda
  }
  t = c(-100, t, 100)
  suv = c(1, exp(-cumsum(lambda)), 0)
  
  return(list(time = t, suv = suv))
}




lss_fn = function(beta, d, weights = rep(1, n)) {
  Y = d$Y; delta = d$delta
  x1 = d$x1; x2 = d$x2; x = cbind(x1, x2)
  xbeta = c(x %*% beta)
  e = log(Y) - xbeta
  es = sort(e)
  ds = delta[order(e)]
  ws = weights[order(e)]
  km = km_dc(es, ds, ws)
  S = approx(km$time, km$suv, es)$y
  F = 1 - S
  lambda = diff(c(0, -log(S)))
  dF = diff(c(0, F))
  denom_R = rev(cumsum(rev(dF)))
  denom_L = F
  num_R = rev(cumsum(rev(es*dF)))
  num_L = cumsum(es*dF)
  ps_L = num_L/pmax(0.001, denom_L)
  ps_R = num_R/pmax(0.001, denom_R)
  ps_L = ps_L[rank(e)] + xbeta
  ps_R = ps_R[rank(e)] + xbeta
  yy = log(Y)*I(delta==1) + ps_R*I(delta==2) + ps_L*I(delta==3)
  beta = lm(yy ~ x1 + x2, weights = weights)$coef[-1]
  
  return(list(beta = beta, time = es, df = F, lambda = lambda))
}


lss_dc = function(d, tp, weights = rep(1,n)) {
  maxiter = 20; error = 10; tol = 1e-5; iter = 0
  beta_old = lm(Y ~ x1 + x2, d)$coef[-1]
  while (iter < maxiter & error > tol) {
    g = lss_fn(beta_old, d, weights = weights)
    beta = g$beta
    error = sum((beta - beta_old) ^ 2)
    iter = iter + 1
    beta_old = beta
  }
  lamb = drop(t(g$lambda) %*% outer(g$time, tp, "<="))
  
  return(list(beta = beta, lamb = lamb, iter = iter))
}


lss_dc_se = function(B = 1, d, tp) {
  out = replicate(B, {
    g = lss_dc(d, tp, weights = rexp(n))
    c(g$beta, g$lamb)
  })
  
  return(apply(out, 1, sd))
}

## Example
beta0 = c(1, 1)
nsim = 1000
B = 200
n = 200
qq = c(0.05, 0.1, 0.2)
par_true = c(beta0, round(-log(1 - qq), 3))
tp = qnorm(qq)

d = simudata(n, beta0 = beta0)
table(d$delta)/n
g = lss_dc(d, tp)
est = c(g$beta, g$lamb)
se = lss_dc_se(B, d, tp)
