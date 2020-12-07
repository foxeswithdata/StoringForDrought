rm(list=ls())

library(tidyverse)
library(plyr)
library(dplyr)
library(ggpubr)

source("code/models/water_model.R")
source("code/parameters.R")

param <- param_list_water_short()

W0_list_maxM <- seq(from=100, to=400, by=50)
W0_list_maxM
maxM_res_list <- lapply(W0_list_maxM, function(x){
  deltat = 0.1
  # Initialise a vector of possible ts values and the length of this
  # vector
  t_seq <- seq(from=0, to=param$Tend, by=deltat)
  simlength <- length(t_seq)
  
  M_out <- vector(length=simlength)
  
  for (i in t_seq){
    out <- simple_water_model_sim_time_breaks(param$M0, param$S0, x, param$Tend, param$kp, param$kr, param$ks, param$kw, i, deltat)
    M_out[i] <- out$M[simlength]
  }
  ts <- which.max(M_out);
  out <- simple_water_model_sim_time_breaks(param$M0, param$S0, x, param$Tend, param$kp, param$kr, param$ks, param$kw, ts, deltat)
  return (list(ts=ts, w0 = x, s=out$S, m=out$M, a=out$A, w=out$W, x=t_seq))
})

maxM_res_list_2 <- lapply(maxM_res_list, function(x) return(l2df(x, length(x$x))))
maxM_res_df <- ldply(maxM_res_list_2, rbind)
maxM_res_df$w0 <- as.factor(maxM_res_df$w0)

p <- ggplot(maxM_res_df, aes(x=x, y=s, group=w0)) +
  geom_line(aes(linetype=w0, color = w0)) +
  xlab("t") +
  ylab("size of storage [gC]")
p





#### TESTING LIVE/DEAD STRATEGIES



W0_list <- seq(from=50, to=500, by=25)
M0_list <- seq(from=50, to=500, by=25)
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

i = 1
for(w in W0_list){
  j = 1
  for(m in M0_list){
    deltat = 0.1
    
    # Initialise a vector of possible ts values and the length of this
    # vector
    t_seq <- seq(from=0, to=param$Tend, by=deltat)
    simlength <- length(t_seq)
    
    S_out <- vector(length=simlength)
    M_out <- vector(length=simlength)
    A_out <- vector(length=simlength)
    W_out <- vector(length=simlength)
    for (l in t_seq){
      out <- simple_water_model_sim_time_breaks(m, param$S0, w, param$Tend, param$kp, param$kr, param$ks, param$kw, l, deltat)
      S_out[l] <- out$S[simlength]
      M_out[l] <- out$M[simlength]
      A_out[l] <- out$A[simlength]
      W_out[l] <- out$W[simlength]
    }
    
    ts_ind <- which.max(M_out)
    ts_maxm[i,j] <- t_seq[ts_ind]
    mT_maxm[i,j] <- M_out[ts_ind]
    sT_maxm[i,j] <- S_out[ts_ind]
    wT_maxm[i,j] <- W_out[ts_ind]
    A_maxm[i,j] <- A_out[ts_ind]
    
    ts_ind <- which.max(S_out)
    ts_maxs[i,j] <- t_seq[ts_ind]
    mT_maxs[i,j] <- M_out[ts_ind]
    sT_maxs[i,j] <- S_out[ts_ind]
    wT_maxs[i,j] <- W_out[ts_ind]
    A_maxs[i,j] <- A_out[ts_ind]
    
    ts_ind <- which.max(M_out+S_out)
    ts_maxms[i,j] <- t_seq[ts_ind]
    mT_maxms[i,j] <- M_out[ts_ind]
    sT_maxms[i,j] <- S_out[ts_ind]
    wT_maxms[i,j] <- W_out[ts_ind]
    A_maxms[i,j] <- A_out[ts_ind]
    j=j+1
  }
  i=i+1
}

x <- W0_list
y <- M0_list
data <- expand.grid(X=x, Y=y)
sT_maxm_vec <- sT_maxm;
dim(sT_maxm_vec) <- c(1,nrow(sT_maxm_vec)* ncol(sT_maxm_vec));
sT_maxm_vec2 <- sT_maxm_vec

sT_maxm_vec2[sT_maxm_vec>=0] = "Live"
sT_maxm_vec2[sT_maxm_vec<0] = "Dead"

data$Z <- t(sT_maxm_vec2);
data$Z <- as.factor(data$Z)

ggplot(data, aes(X, Y, z= Z)) + geom_tile(aes(fill = Z)) + theme_bw() +
  xlab("Initial water input W0 [kgH20]") +
  ylab("Initial plant biomass")



