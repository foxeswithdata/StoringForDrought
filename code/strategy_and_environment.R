source("code/parameters.R")
source("code/models/water_model.R")

kf_list <- seq(from=0, to=1, by=0.01) # I will change this once I have the method figured out.
W0_list <- seq(from=50, to=500, by=10) # as above

param <- param_list_water_short2()

ts_out <- matrix(nrow=length(kf_list), ncol=length(W0_list))

deltat = 0.001
# Initialise a vector of possible ts values and the length of this
# vector
t_seq <- seq(from=0, to=param$Tend, by=deltat)
simlength <- length(t_seq)

j=1
for(w in W0_list){
  S_out <- vector(length=simlength)
  M_out <- vector(length=simlength)
  k=1
  
  for (i in t_seq){
    # Find the output values of simulation for the specific ts value
    param_2 <- param
    param_2$W0 <- w
    param_2$x0 <- c(param_2$M0, param_2$S0, param_2$W0)
    
    out <- simple_water_model_sim_time_breaks(param_2, i, deltat)
    
    # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
    S_out[k] <- out$S[simlength]
    M_out[k] <- out$M[simlength]
    
    if(sum(out$S<0) > 0){
      S_out[k] <- NA
      M_out[k] <- NA
    }
    k = k+1
  }
  
  i = 1
  for (kf in kf_list){
    ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
    print(ms_sum[floor(simlength/2)])
    ts_ind <- which.max(ms_sum)
    if(length(ts_ind) == 0){
      ts_out[i,j] = NA
    }
    else{
      print("Ts ind")
      print(ts_ind)
      ts_out[i,j] = t_seq[ts_ind]
    }
    i = i+1
  }
  j = j+1
}



