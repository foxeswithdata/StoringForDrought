#### TESTING THE FORWARD-BACKWARDS MODEL

rm(list=ls())

library(ggplot2)

source('code/oct/forward_backward.R')
source("code/models/water_model.R")

S0 = 10;
M0 = 100;
W0 = 250;
X0 = c(M0 ,S0 ,W0);
T_end = 64;
kp = 0.04;
ks = 0.1;
kr = 0.02;
delta = 1;
ubound = c(0, ks);
kpbound = c(0, kp);
kw = 1

w50 = 2
k = 10

output_FB <- forwardback_kp_log_RK(kp, kr, kw, X0, T_end, ubound, k = 10, w50 = 2)

# plot(output$time, output$x[1,], type='l', col="red", ylim=c(-40,W0))
# lines(output$time, output$x[2,], col="green")
# lines(output$time, output$x[3,], col="blue")


# Set simulation time tolerance. The smaller this value the 
# more precise the simulation
deltat = 0.05

# Initialise a vector of possible ts values and the length of this
# vector
t_seq <- seq(from=0, to=T_end, by=deltat)
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
  out <- simple_water_model_sim_time_breaks(M0, S0, W0, T_end, 
                                            kp, kr, ks, kw, 
                                            i, deltat)
  
  out_sig <- water_model_log(M0, S0, W0, T_end, 
                             kp, kr, ks, kw, 
                             i, deltat, w_mid = w50, growth_rate = k)
  
  
  
  # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
  S_out[k] <- out$S[simlength]
  M_out[k] <- out$M[simlength]
  A_out[k] <- out$A[simlength]
  W_out[k] <- out$W[simlength]
  
  S_out_sig[k] <- out_sig$S[simlength]
  M_out_sig[k] <- out_sig$M[simlength]
  A_out_sig[k] <- out_sig$A[simlength]
  W_out_sig[k] <- out_sig$W[simlength]
  
  if(sum(out$S<0) > 0){
    S_out[k] <- NA
    M_out[k] <- NA
    A_out[k] <- NA
    W_out[k] <- NA
  }
  
  if(sum(out_sig$S<0) > 0){
    S_out_sig[k] <- NA
    M_out_sig[k] <- NA
    A_out_sig[k] <- NA
    W_out_sig[k] <- NA
  }
  
  k = k+1
}

ts_sig <- which.max(S_out_sig)

out_sig <- water_model_log(M0, S0, W0, T_end, 
                           kp, kr, ks, kw, 
                           ts_sig, deltat, w_mid = w50, growth_rate = k)

ts <- which.max(S_out)
out <- simple_water_model_sim_time_breaks(M0, S0, W0, T_end, 
                                          kp, kr, ks, kw, 
                                          ts, deltat)


plot(output_FB$time, output_FB$x[1,], type='l', col="red", ylim=c(-40,W0))
lines(output_FB$time, output_FB$x[2,], col="green")
lines(output_FB$time, output_FB$x[3,], col="blue")
lines(out_sig$t, out_sig$M, col="red", lty=4)
lines(out_sig$t, out_sig$S, col="green", lty=4)
lines(out_sig$t, out_sig$W, col="blue", lty=4)
lines(out$t, out$M, col="red", lty=3)
lines(out$t, out$S, col="green", lty=3)
lines(out$t, out$W, col="blue", lty=3)

