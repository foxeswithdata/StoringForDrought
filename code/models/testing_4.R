x0 <- c(0.25,0.5,0.75,1,2,5,10)
L = 0.04
k = 10

w <- seq(from = 0, to = 5, by = 0.05)
kp <- L/(1 + exp(-k * (w - x0[5])))

plot(kp ~ w)

### Test whether comparable to original model

library(tidyverse)
library(plyr)
library(dplyr)
library(ggpubr)

source("./code/models/water_model.R")
source("./code/helper_functions.R")
source("./code/parameters.R")

param <- param_list_water_short2()

  # Set simulation time tolerance. The smaller this value the 
  # more precise the simulation
  deltat = 0.05
  
  # Initialise a vector of possible ts values and the length of this
  # vector
  t_seq <- seq(from=0, to=param$Tend, by=deltat)
  simlength <- length(t_seq)
  
  # Initialise output vectors, these will contain the 
  # Final outputs of simulations for each ts value
  S_out <- vector(length=simlength)
  M_out <- vector(length=simlength)
  A_out <- vector(length=simlength)
  W_out <- vector(length=simlength)
  
  S_out_sig <- vector(length=simlength)
  M_out_sig <- vector(length=simlength)
  A_out_sig <- vector(length=simlength)
  W_out_sig <- vector(length=simlength)
  
  
  # Iterate through the possible ts values
  k = 1
  for (i in t_seq){
    # Find the output values of simulation for the specific ts value
    out <- simple_water_model_sim_time_breaks(param$M0, param$S0, param$W0, param$Tend, 
                                              param$kp, param$kr, param$ks, param$kw, 
                                              i, deltat)
    
    out_sig <- water_model_log(param$M0, param$S0, param$W0, param$Tend, 
                                              param$kp, param$kr, param$ks, param$kw, 
                                              i, deltat, w_mid = 2, growth_rate = 10)
    
    
    
    # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
    S_out[k] <- out$S[simlength]
    M_out[k] <- out$M[simlength]
    A_out[k] <- out$A[simlength]
    W_out[k] <- out$W[simlength]
    
    S_out_sig[k] <- out_sig$S[simlength]
    M_out_sig[k] <- out_sig$M[simlength]
    A_out_sig[k] <- out_sig$A[simlength]
    W_out_sig[k] <- out_sig$W[simlength]
    

    k = k+1
  }
  
  
  plot(t_seq, S_out_sig, type="l", col="red")
  lines(t_seq, S_out_sig, col="blue")
  
  # Try ts = 10
  
  ts_test = 10
  
  out <- simple_water_model_sim_time_breaks(param$M0, param$S0, param$W0, param$Tend, 
                                            param$kp, param$kr, param$ks, param$kw, 
                                            ts_test, deltat)
  
  out_sig <- water_model_log(param$M0, param$S0, param$W0, param$Tend, 
                             param$kp, param$kr, param$ks, param$kw, 
                             ts_test, deltat, w_mid = 2, growth_rate = 10)

  
  
  plot(out$t, out$W, type = "l", col = "red")  
  lines(out_sig$t, out_sig$W, col = "blue")

  plot(out$S - out_sig$S)  
  
  deltat = 0.05
  ts_test = 20
  
  out <- simple_water_model_sim_time_breaks(param$M0, param$S0, param$W0, param$Tend, 
                                            param$kp, param$kr, param$ks, param$kw, 
                                            ts_test, deltat)
  
  out_sig <- water_model_log(param$M0, param$S0, param$W0, param$Tend, 
                             param$kp, param$kr, param$ks, param$kw, 
                             ts_test, deltat, w_mid = 2, growth_rate = 10)
  
  
  
  plot(out_sig$t, out_sig$W, type = "l", col = "red")  
  lines(out$t, out$W, col = "blue")
  
  plot(abs(out$S - out_sig$S)  )
  