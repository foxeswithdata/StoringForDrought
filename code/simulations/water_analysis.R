rm(list=ls())

library(tidyverse)
library(plyr)
library(dplyr)
library(ggpubr)


# 
source("code/model.R")
source("code/water_model.R")

# FIND LAGRANGE MULTIPLIER LAMBDA

S0 = 10 #gC
M0 = 100 #gC
Tend = 64
kp = 0.04 #gCgC-1d-1
kmax = 0.1 
kr = 0.02 #gCgC-1d-1
W0 = 400 #kgH20
wue = 1 #kgH20gC-1



###VARY W0

W0_list <- seq(from=10, to=500, by=10)
ts_maxm <- 0
mT_maxm <- 0
sT_maxm <- 0
wT_maxm <- 0
A_maxm <- 0

ts_maxs <- 0
mT_maxs <- 0
sT_maxs <- 0
wT_maxs <- 0
A_maxs <- 0

ts_maxms <- 0
mT_maxms <- 0
sT_maxms <- 0
wT_maxms <- 0
A_maxms <- 0

for(w in W0_list){
  S_out <- vector(length=Tend)
  M_out <- vector(length=Tend)
  A_out <- vector(length=Tend)
  W_out <- vector(length=Tend)
  for (i in 1:Tend){
    u <- c(rep(kmax, i), rep(0, Tend-i))
    out <- simple_water_model_sim(M0, S0, w, Tend, kp, kr, u, wue)
    S_out[i] <- out$S[Tend]
    M_out[i] <- out$M[Tend]
    A_out[i] <- sum(out$A)
    W_out[i] <- out$W[Tend]
  }
  i = which(W0_list==w);
  ts_maxm[i] <- which.max(M_out)
  ts = ts_maxm[i]
  u <- c(rep(kmax, ts), rep(0, Tend-ts))
  out <- simple_water_model_sim(M0, S0, w, Tend, kp, kr, u, wue)
  mT_maxm[i] <- out$M[Tend]
  sT_maxm[i] <- out$S[Tend]
  wT_maxm[i] <- out$W[Tend]
  A_maxm[i] <- sum(out$A)
  
  ts_maxs[i] <- which.max(S_out)
  ts = ts_maxs[i]
  u <- c(rep(kmax, ts), rep(0, Tend-ts))
  out <- simple_water_model_sim(M0, S0, w, Tend, kp, kr, u, wue)
  mT_maxs[i] <- out$M[Tend]
  sT_maxs[i] <- out$S[Tend]
  wT_maxs[i] <- out$W[Tend]
  A_maxs[i] <- sum(out$A)
  
  ts_maxms[i] <- which.max(M_out+S_out)
  ts = ts_maxms[i]
  u <- c(rep(kmax, ts), rep(0, Tend-ts))
  out <- simple_water_model_sim(M0, S0, w, Tend, kp, kr, u, wue)
  mT_maxms[i] <- out$M[Tend]
  sT_maxms[i] <- out$S[Tend]
  wT_maxms[i] <- out$W[Tend]
  A_maxms[i] <- sum(out$A)
}

maxM_strategy_w0 <- data.frame(w0 = W0_list, mT = mT_maxm, ts = ts_maxm, sT = sT_maxm, a = A_maxm, wT = wT_maxm)
maxS_strategy_w0 <- data.frame(w0 = W0_list, mT = mT_maxs, ts = ts_maxs, sT = sT_maxs, a = A_maxs, wT = wT_maxs)
maxMS_strategy_w0 <- data.frame(w0 = W0_list, mT = mT_maxms, ts = ts_maxms, sT = sT_maxms, a = A_maxms, wT = wT_maxms)

water_analysis_figures_w0(maxM_strategy_w0, "figures/oct_maxm_w0")
water_analysis_figures_w0(maxS_strategy_w0, "figures/oct_maxs_w0")
water_analysis_figures_w0(maxMS_strategy_w0, "figures/oct_maxms_w0")

plot(ts_maxm ~ W0_list)
plot(W0_list ~ ts_maxs)
plot(W0_list ~ ts_maxms)

w0_ts_strategies_df <- data.frame(w0 = W0_list, maxm_ts = ts_maxm, maxs_ts = ts_maxs, maxms_ts = ts_maxms)

#plot maxm results

p <- ggplot(w0_ts_strategies_df, aes(x=w0)) +
  geom_line(aes(y=maxm_ts, colour="MaxM")) + 
  geom_line(aes(y=maxs_ts,colour="MaxS")) +
  geom_line(aes(y=maxms_ts, colour="MaxM+S")) +
  xlab("Total rain input [kgH20]") + 
  ylab("Optimal growth switch day ts") +
  scale_colour_manual("", 
                      values = c("MaxM"="#0072B2", "MaxS"="#D55E00", 
                                 "MaxM+S"="#CC79A7")) 
p


### VARY T_END

T_list <- 32:128
ts_maxm <- 0
mT_maxm <- 0
sT_maxm <- 0
wT_maxm <- 0
A_maxm <- 0

ts_maxs <- 0
mT_maxs <- 0
sT_maxs <- 0
wT_maxs <- 0
A_maxs <- 0

ts_maxms <- 0
mT_maxms <- 0
sT_maxms <- 0
wT_maxms <- 0
A_maxms <- 0

w0=240

for(tend in T_list){
  S_out <- vector(length=tend)
  M_out <- vector(length=tend)
  A_out <- vector(length=tend)
  W_out <- vector(length=tend)
  for (i in 1:tend){
    u <- c(rep(kmax, i), rep(0, tend-i))
    out <- simple_water_model_sim(M0, S0, w0, tend, kp, kr, u, wue)
    S_out[i] <- out$S[tend]
    M_out[i] <- out$M[tend]
    A_out[i] <- sum(out$A)
    W_out[i] <- out$W[tend]
  }
  i = which(T_list==tend);
  ts_maxm[i] <- which.max(M_out)
  ts = ts_maxm[i]
  u <- c(rep(kmax, ts), rep(0, tend-ts))
  out <- simple_water_model_sim(M0, S0, w0, tend, kp, kr, u, wue)
  mT_maxm[i] <- out$M[tend]
  sT_maxm[i] <- out$S[tend]
  wT_maxm[i] <- out$W[tend]
  A_maxm[i] <- sum(out$A)
  
  ts_maxs[i] <- which.max(S_out)
  ts = ts_maxs[i]
  u <- c(rep(kmax, ts), rep(0, tend-ts))
  out <- simple_water_model_sim(M0, S0, w0, tend, kp, kr, u, wue)
  mT_maxs[i] <- out$M[tend]
  sT_maxs[i] <- out$S[tend]
  wT_maxs[i] <- out$W[tend]
  A_maxs[i] <- sum(out$A)
  
  ts_maxms[i] <- which.max(M_out+S_out)
  ts = ts_maxms[i]
  u <- c(rep(kmax, ts), rep(0, tend-ts))
  out <- simple_water_model_sim(M0, S0, w0, tend, kp, kr, u, wue)
  mT_maxms[i] <- out$M[tend]
  sT_maxms[i] <- out$S[tend]
  wT_maxms[i] <- out$W[tend]
  A_maxms[i] <- sum(out$A)
}

maxM_strategy_tend <- data.frame(tend = T_list, mT = mT_maxm, ts = ts_maxm, sT = sT_maxm, a = A_maxm, wT = wT_maxm)
maxS_strategy_tend <- data.frame(tend = T_list, mT = mT_maxs, ts = ts_maxs, sT = sT_maxs, a = A_maxs, wT = wT_maxs)
maxMS_strategy_tend <- data.frame(tend = T_list, mT = mT_maxms, ts = ts_maxms, sT = sT_maxms, a = A_maxms, wT = wT_maxms)

water_analysis_figures_tend(maxM_strategy_tend, "figures/oct_maxm_tend")
water_analysis_figures_tend(maxS_strategy_tend, "figures/oct_maxs_tend")
water_analysis_figures_tend(maxMS_strategy_tend, "figures/oct_maxms_tend")

### VARY kp

kp_list <- seq(from=0.03, to=0.1, by=0.005)
ts_maxm <- 0
mT_maxm <- 0
sT_maxm <- 0
wT_maxm <- 0
A_maxm <- 0

ts_maxs <- 0
mT_maxs <- 0
sT_maxs <- 0
wT_maxs <- 0
A_maxs <- 0

ts_maxms <- 0
mT_maxms <- 0
sT_maxms <- 0
wT_maxms <- 0
A_maxms <- 0

w0=240
Tend = 64

for(kp in kp_list){
  S_out <- vector(length=Tend)
  M_out <- vector(length=Tend)
  A_out <- vector(length=Tend)
  W_out <- vector(length=Tend)
  for (i in 1:Tend){
    u <- c(rep(kmax, i), rep(0, Tend-i))
    out <- simple_water_model_sim(M0, S0, w0, Tend, kp, kr, u, wue)
    S_out[i] <- out$S[Tend]
    M_out[i] <- out$M[Tend]
    A_out[i] <- sum(out$A)
    W_out[i] <- out$W[Tend]
  }
  i = which(kp_list == kp);
  ts_maxm[i] <- which.max(M_out)
  ts = ts_maxm[i]
  u <- c(rep(kmax, ts), rep(0, Tend-ts))
  out <- simple_water_model_sim(M0, S0, w0, Tend, kp, kr, u, wue)
  mT_maxm[i] <- out$M[Tend]
  sT_maxm[i] <- out$S[Tend]
  wT_maxm[i] <- out$W[Tend]
  A_maxm[i] <- sum(out$A)
  
  ts_maxs[i] <- which.max(S_out)
  ts = ts_maxs[i]
  u <- c(rep(kmax, ts), rep(0, Tend-ts))
  out <- simple_water_model_sim(M0, S0, w0, Tend, kp, kr, u, wue)
  mT_maxs[i] <- out$M[Tend]
  sT_maxs[i] <- out$S[Tend]
  wT_maxs[i] <- out$W[Tend]
  A_maxs[i] <- sum(out$A)
  
  ts_maxms[i] <- which.max(M_out+S_out)
  ts = ts_maxms[i]
  u <- c(rep(kmax, ts), rep(0, Tend-ts))
  out <- simple_water_model_sim(M0, S0, w0, Tend, kp, kr, u, wue)
  mT_maxms[i] <- out$M[Tend]
  sT_maxms[i] <- out$S[Tend]
  wT_maxms[i] <- out$W[Tend]
  A_maxms[i] <- sum(out$A)
}

maxM_strategy_kp <- data.frame(kp = kp_list, mT = mT_maxm, ts = ts_maxm, sT = sT_maxm, a = A_maxm, wT = wT_maxm)
maxS_strategy_kp <- data.frame(kp = kp_list, mT = mT_maxs, ts = ts_maxs, sT = sT_maxs, a = A_maxs, wT = wT_maxs)
maxMS_strategy_kp <- data.frame(kp = kp_list, mT = mT_maxms, ts = ts_maxms, sT = sT_maxms, a = A_maxms, wT = wT_maxms)

water_analysis_figures_kp(maxM_strategy_kp, "figures/oct_maxm_kp")
water_analysis_figures_kp(maxS_strategy_kp, "figures/oct_maxs_kp")
water_analysis_figures_kp(maxMS_strategy_kp, "figures/oct_maxms_kp")

### Other things



W0_list_maxM <- seq(from=100, to=400, by=50)
W0_list_maxM
maxM_res_list <- lapply(W0_list_maxM, function(x){
  M_out <- vector(length=Tend)
  for (i in 1:Tend){
    u <- c(rep(kmax, i), rep(0, Tend-i))
    out <- simple_water_model_sim(M0, S0, x, Tend, kp, kr, u, wue)
    M_out[i] <- out$M[Tend]
  }
  ts <- which.max(M_out);
  print(ts)
  u <- c(rep(kmax, ts), rep(0, Tend-ts))
  out <- simple_water_model_sim(M0, S0, x, Tend, kp, kr, u, wue)
  return (list(ts=ts, w0 = x, s=out$S, m=out$M, a=out$A, w=out$W, x=1:Tend))
})

maxM_res_list_2 <- lapply(maxM_res_list, function(x) return(l2df(x, Tend)))
maxM_res_df <- ldply(maxM_res_list_2, rbind)
maxM_res_df$w0 <- as.factor(maxM_res_df$w0)

p <- ggplot(maxM_res_df, aes(x=x, y=s, group=w0)) +
  geom_line(aes(linetype=w0, color = w0)) + 
  xlab("t") + 
  ylab("size of storage [gC]")
p

p <- ggplot(maxM_res_df, aes(x=x, y=m, group=w0)) +
  geom_line(aes(linetype=w0, color = w0)) + 
  xlab("t") + 
  ylab("size of biomass [gC]")
p



W0_list_maxS <- seq(from=100, to=400, by=50)
W0_list_maxS
maxS_res_list <- lapply(W0_list_maxS, function(x){
  S_out <- vector(length=Tend)
  for (i in 1:Tend){
    u <- c(rep(kmax, i), rep(0, Tend-i))
    out <- simple_water_model_sim(M0, S0, x, Tend, kp, kr, u, wue)
    S_out[i] <- out$S[Tend]
  }
  ts <- which.max(S_out);
  print(ts)
  u <- c(rep(kmax, ts), rep(0, Tend-ts))
  out <- simple_water_model_sim(M0, S0, x, Tend, kp, kr, u, wue)
  return (list(ts=ts, w0 = x, s=out$S, m=out$M, a=out$A, w=out$W, x=1:Tend))
})

maxS_res_list_2 <- lapply(maxS_res_list, function(x) return(l2df(x, Tend)))
maxS_res_df <- ldply(maxS_res_list_2, rbind)
maxS_res_df$w0 <- as.factor(maxS_res_df$w0)

p <- ggplot(maxS_res_df, aes(x=x, y=s, group=w0)) +
  geom_line(aes(linetype=w0, color = w0)) + 
  xlab("t") + 
  ylab("size of storage [gC]")
p

p <- ggplot(maxS_res_df, aes(x=x, y=m, group=w0)) +
  geom_line(aes(linetype=w0, color = w0)) + 
  xlab("t") + 
  ylab("size of biomass [gC]")
p



#### MAPPING OUT LIVE AND DEAD STRATEGIES 

S0 = 10 #gC
M0 = 100 #gC
Tend = 64
kp = 0.04 #gCgC-1d-1
kmax = 0.1 
kr = 0.02 #gCgC-1d-1
W0 = 400 #kgH20
wue = 1 #kgH20gC-1

W0_list <- seq(from=10, to=500, by=10)
M0_list <- seq(from=10, to=500, by=10)
ts_maxm <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
mT_maxm <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
sT_maxm <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
wT_maxm <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
A_maxm <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))

ts_maxs <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
mT_maxs <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
sT_maxs <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
wT_maxs <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
A_maxs <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))

ts_maxms <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
mT_maxms <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
sT_maxms <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
wT_maxms <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))
A_maxms <- matrix(0, nrow = length(W0_list), ncol = length(M0_list))

for(w in W0_list){
  for(m in M0_list){
    S_out <- vector(length=Tend)
    M_out <- vector(length=Tend)
    A_out <- vector(length=Tend)
    W_out <- vector(length=Tend)
    for (i in 1:Tend){
      u <- c(rep(kmax, i), rep(0, Tend-i))
      out <- simple_water_model_sim(m, S0, w, Tend, kp, kr, u, wue)
      S_out[i] <- out$S[Tend]
      M_out[i] <- out$M[Tend]
      A_out[i] <- sum(out$A)
      W_out[i] <- out$W[Tend]
    }
    i = which(W0_list==w);
    j = which(M0_list==m);
    ts_maxm[i,j] <- which.max(M_out)
    ts = ts_maxm[i,j]
    u <- c(rep(kmax, ts), rep(0, Tend-ts))
    out <- simple_water_model_sim(m, S0, w, Tend, kp, kr, u, wue)
    mT_maxm[i,j] <- out$M[Tend]
    sT_maxm[i,j] <- out$S[Tend]
    wT_maxm[i,j] <- out$W[Tend]
    A_maxm[i,j] <- sum(out$A)
    
    ts_maxs[i,j] <- which.max(S_out)
    ts = ts_maxs[i,j]
    u <- c(rep(kmax, ts), rep(0, Tend-ts))
    out <- simple_water_model_sim(m, S0, w, Tend, kp, kr, u, wue)
    mT_maxs[i,j] <- out$M[Tend]
    sT_maxs[i,j] <- out$S[Tend]
    wT_maxs[i,j] <- out$W[Tend]
    A_maxs[i,j] <- sum(out$A)
    
    ts_maxms[i,j] <- which.max(M_out+S_out)
    ts = ts_maxms[i,j]
    u <- c(rep(kmax, ts), rep(0, Tend-ts))
    out <- simple_water_model_sim(m, S0, w, Tend, kp, kr, u, wue)
    mT_maxms[i,j] <- out$M[Tend]
    sT_maxms[i,j] <- out$S[Tend]
    wT_maxms[i,j] <- out$W[Tend]
    A_maxms[i,j] <- sum(out$A)
  }
}


x <- W0_list
y <- M0_list
data <- expand.grid(X=x, Y=y)
ts_maxm_vec <- ts_maxm;
dim(ts_maxm_vec) <- c(1,nrow(ts_maxm_vec)* ncol(ts_maxm_vec));
data$Z <- t(ts_maxm_vec);

library(ggplot2)
ggplot(data, aes(X, Y, z= Z)) + geom_tile(aes(fill = Z)) + theme_bw() +
  xlab("Initial water input W0 [kgH20]") + 
  ylab("Initial plant biomass")
