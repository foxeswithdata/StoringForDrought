rm(list=ls())

source("code/model.R")
source("code/lagrange_sol_max_ent.R")

library(tidyverse)

# FIND LAGRANGE MULTIPLIER LAMBDA

S0 = 0
M0 = 0.1
T_end = 64
kp = 1
kmax = 0.1
kr = 0.5
tcrit_av = 48;

S0 = 10 #gC
M0 = 100 #gC
T_end = 64
kp_1 = 0.04 #gCgC-1d-1
kmax = 0.1 
kr = 0.02 #gCgC-1d-1

# maxts is the last point at which ts does not produce negative 
maxts = 40

kp = c(rep(kp_1, times=tcrit_av), rep(0, times=T_end-tcrit_av))
#MINIMUM ts = 37 given no restrictions on ST


# GENERATE RESULTS FOR ONLY LIVE PLANTS ts = 30
ts = 30
u = c(rep(kmax,times=ts), rep(0, times=T_end-ts))

out <- simplemodel_sim(M0, S0, T_end, kp, kr,u)
Itot = sum(out$A);
Itot_arr <- c(10, 25, 50, 100, 200)
Itot_arr_long <- seq(from=10, to=250, by=10)
Itot_arr_dead <- c(50,100,200, 400, 600, 800, 1000)
Itot_arr_dead_long <- seq(from=50, to=1000, by=25)

min = -10
max = 10

#Find the roots of the function
root <- lapply(Itot_arr, function (z) {
  Itot <- z
  uniroot(function(x) {
    find_lambda_err_alive(M0,S0,T_end,kp,kr,kmax,maxts,x,Itot)}, 
    interval=c(min,max), tol=.Machine$double.eps)
  })
lambda <- sapply(root, function(x) return(x$root))

root_long <- lapply(Itot_arr_long, function (z) {
  Itot <- z
  uniroot(function(x) {
    find_lambda_err_alive(M0,S0,T_end,kp,kr,kmax,maxts,x,Itot)}, 
    interval=c(min,max), tol=.Machine$double.eps)
})
lambda_long <- sapply(root_long, function(x) return(x$root))


root_dead <- lapply(Itot_arr_dead, function (z) {
  Itot <- z
  uniroot(function(x) {find_lambda_err(M0,S0,T_end,kp,kr,kmax,x,Itot)}, 
    interval=c(min,max), tol=.Machine$double.eps, extendInt = "yes")
})
lambda_dead <- sapply(root_dead, function(x) return(x$root))

root_dead_long <- lapply(Itot_arr_dead_long, function (z) {
  Itot <- z
  uniroot(function(x) {find_lambda_err(M0,S0,T_end,kp,kr,kmax,x,Itot)}, 
          interval=c(min,max), tol=.Machine$double.eps, extendInt = "yes")
})
lambda_dead_long <- sapply(root_dead_long, function(x) return(x$root))


# CALCULATE PK

# pk_dist <- sapply(1:maxts, function(x){pk_res_alive(M0,S0,T_end,kp,kr,kmax, maxts, x, lambda)})
pk_dist_list <- lapply(lambda, function(x){
  lam <- x
  Itot <- Itot_arr[lambda==lam]
  pk_dist <- sapply(1:maxts, function(y){pk_res_alive(M0,S0,T_end,kp,kr,kmax, maxts, y, x)})
  entropy <- -sum(pk_dist*log(pk_dist))
  ts <- 1:maxts
  Ek <- sapply(1:maxts, function(z){
    u = c(rep(kmax,times=z), rep(0, times=T_end-z))
    out=simplemodel_sim(M0, S0, T_end, kp, kr, u);
    return (sum(out$A))
  })
  # df <- data.frame(ts=ts, pk_distribution=pk_dist, Ek=Ek)
  return (list(lambda=lam, Itot=Itot, pk_distribiution=pk_dist, entropy=entropy, ts=ts, Ek=Ek))#, df=df))
})

pk_dist_list_long <- lapply(lambda_long, function(x){
  lam <- x
  Itot <- Itot_arr_long[lambda_long==lam]
  pk_dist <- sapply(1:maxts, function(y){pk_res_alive(M0,S0,T_end,kp,kr,kmax, maxts, y, x)})
  entropy <- -sum(pk_dist*log(pk_dist))
  variance <- var(pk_dist)
  ts <- 1:maxts
  Ek <- sapply(1:maxts, function(z){
    u = c(rep(kmax,times=z), rep(0, times=T_end-z))
    out=simplemodel_sim(M0, S0, T_end, kp, kr, u);
    return (sum(out$A))
  })
  # df <- data.frame(ts=ts, pk_distribution=pk_dist, Ek=Ek)
  return (list(lambda=lam, Itot=Itot, pk_distribiution=pk_dist, entropy=entropy, ts=ts, Ek=Ek, variance=variance))#, df=df))
})

pk_dist_list_dead <- lapply(lambda_dead, function(x){
  lam <- x
  Itot <- Itot_arr_dead[lambda_dead==lam]
  pk_dist <- sapply(1:T_end, function(y){pk_res(M0,S0,T_end,kp,kr,kmax, y, x)})
  entropy <- -sum(pk_dist*log(pk_dist))
  ts <- 1:T_end
  Ek <- sapply(1:T_end, function(z){
    u = c(rep(kmax,times=z), rep(0, times=T_end-z))
    out=simplemodel_sim(M0, S0, T_end, kp, kr, u);
    return (sum(out$A))
  })
  # df <- data.frame(ts=ts, pk_distribution=pk_dist, Ek=Ek)
  return (list(lambda=lam, Itot=Itot, pk_distribiution=pk_dist, entropy=entropy, ts=ts, Ek=Ek))#, df=df))
})

pk_dist_list_dead_long <- lapply(lambda_dead_long, function(x){
  lam <- x
  Itot <- Itot_arr_dead_long[lambda_dead_long==lam]
  pk_dist <- sapply(1:T_end, function(y){pk_res(M0,S0,T_end,kp,kr,kmax, y, x)})
  entropy <- -sum(pk_dist*log(pk_dist))
  ts <- 1:T_end
  Ek <- sapply(1:T_end, function(z){
    u = c(rep(kmax,times=z), rep(0, times=T_end-z))
    out=simplemodel_sim(M0, S0, T_end, kp, kr, u);
    return (sum(out$A))
  })
  variance <- var(pk_dist)
  # df <- data.frame(ts=ts, pk_distribution=pk_dist, Ek=Ek)
  return (list(lambda=lam, Itot=Itot, pk_distribiution=pk_dist, entropy=entropy, ts=ts, Ek=Ek, variance=variance))#, df=df))
})


# FIND 2 ELEMENTS FOR FIGURE

ts = 15
u = c(rep(kmax,times=ts), rep(0, times=T_end-ts))
out <- simplemodel_sim(M0, S0, T_end, kp, kr,u)
Itot = sum(out$A);

Itot1 <- 220
Itot2 <- 260

root1 <- uniroot(function(x) {
  find_lambda_err_alive(M0,S0,T_end,kp,kr,kmax,maxts,x,Itot1)}, 
  interval=c(-0.1,0.1), tol=.Machine$double.eps, extendInt = "yes")

lambda1 <- root1$root

root2 <- uniroot(function(x) {
  find_lambda_err_alive(M0,S0,T_end,kp,kr,kmax,maxts,x,Itot2)}, 
  interval=c(-0.1,0.1), tol=.Machine$double.eps, extendInt = "yes")

lambda2 <- root2$root

pk_dist1 <- sapply(1:maxts, function(y){pk_res_alive(M0,S0,T_end,kp,kr,kmax, maxts, y, lambda1)})
pk_dist2 <- sapply(1:maxts, function(y){pk_res_alive(M0,S0,T_end,kp,kr,kmax, maxts, y, lambda2)})

set1_df <- data.frame(pk_dist = pk_dist1, ts =1:maxts);
set2_df <- data.frame(pk_dist = pk_dist2, ts =1:maxts);

p <- ggplot(set1_df, aes(x=ts, y=pk_dist)) +
  geom_line(color="#6969B3", size=2) + 
  theme_bw() +
  xlab(expression("Allocation switch day " ~ t[s] * " (d)" )) + 
  ylab(expression("Probability of growth switch day " ~ p[t[s]])) + 
  expand_limits(y=c(0,0.05))
p

p <- ggplot(set2_df, aes(x=ts, y=pk_dist)) +
  geom_line(color="#A6D3A0", size=2) + 
  theme_bw() +
  xlab(expression("Allocation switch day " ~ t[s] * " (d)" )) + 
  ylab(expression("Probability of growth switch day " ~ p[t[s]])) +
  expand_limits(y=c(0,0.05))
p


p <- ggplot() +
  geom_line(aes(x=1:T_end, y=kp), color="#656565", size=1) + 
  theme_bw() +
  xlab(expression("day t (d)" )) + 
  ylab(expression("Photosynthesis parameter " ~k[p] ~ "gC" ~ g^-1 ~ d^-1)) +
  expand_limits(y=c(0,0.05))
p
