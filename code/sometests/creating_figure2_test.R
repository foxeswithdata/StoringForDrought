library(ggplot2)
library(wesanderson)
library(ggpubr)


source("code/parameters.R")
source("code/models/water_model.R")

param <- param_list_water_short2()

W0_list <- seq(from=100, to=500, by=10)

deltat = 0.1
# Initialise a vector of possible ts values and the length of this
# vector
t_seq <- seq(from=0, to=param$Tend, by=deltat)
simlength <- length(t_seq)

ts_maxM <- vector(length = length(W0_list))
ts_maxS <- vector(length = length(W0_list))
ts_maxMS <- vector(length = length(W0_list))

j = 1
for(w in W0_list){
  S_out <- vector(length=simlength)
  M_out <- vector(length=simlength)
  k=1
  for (i in t_seq){
    # Find the output values of simulation for the specific ts value
    
    out <- simple_water_model_sim_time_breaks(param, i, deltat)
    
    # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
    S_out[k] <- out$S[simlength]
    M_out[k] <- out$M[simlength]
    
    if(sum(out$S<0) > 0){
      S_out[k] <- NA
      M_out[k] <- NA
    }
    k = k+1
  }
  
  kf = 1
  ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
  ts_ind <- which.max(ms_sum)
  if(length(ts_ind) == 0){
    ts_maxM[j] = NA
  }
  else{
    ts_maxM[j] = t_seq[ts_ind]
  }
  
  
  kf = 0
  ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
  ts_ind <- which.max(ms_sum)
  if(length(ts_ind) == 0){
    ts_maxS[j] = NA
  }
  else{
    ts_maxS[j] = t_seq[ts_ind]
  }
  
  kf = 0.5
  ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
  ts_ind <- which.max(ms_sum)
  if(length(ts_ind) == 0){
    ts_maxMS[j] = NA
  }
  else{
    ts_maxMS[j] = t_seq[ts_ind]
  }
  
  j = j+1
}


