#Step 0: prepare workstation
rm(list = ls())
source("code/model.R")
source("code/generate_figures.R")


#Step 1 : find <tcrit> and 'Itotal' (here equivalent to Atotal)

# Define model parameters
S0 = 0
M0 = 0.1
T_end = 64
kp = 1
kmax = 0.1
kr = 0.5
tcrit_av = 48;

S0 = 1 #gC
M0 = 10 #gC
T_end = 64
kp = 0.04 #gCgC-1d-1
kmax = 0.1 
kr = 0.02 #gCgC-1d-1

#Calculate deterministic model outputs for the parameters
ts = find_ts_SimpleModel1(M0,S0,T_end,tcrit_av,kp,kr,kmax,0,100,0.001);
out = simpleModel_1(M0,S0,T_end,tcrit_av,ts,kp,kr,kmax);

plot_biomass_storage(out$M, out$S)

ts = 1:(T_end-1)
out_sT <- sapply(ts, function(x){
  out = simpleModel_1(M0,S0,T_end,tcrit_av,x,kp,kr,kmax);
  return(out$S[length(out$S)])
})
out_mT <- sapply(ts, function(x){
  out = simpleModel_1(M0,S0,T_end,tcrit_av,x,kp,kr,kmax);
  return(out$M[length(out$M)])
})


plot(out_mT[out_sT>=0], type="l", xlab="t_s (day)", ylab="final biomass M(T)")
plot(out_sT[out_sT>=0], type="l", xlab="t_s (day)", ylab="final storage S(T)")



out_maxS <- simpleModel_1(M0,S0,T_end,tcrit_av,30,kp,kr,kmax)
plot_biomass_storage(out_maxS$M, out_maxS$S)
