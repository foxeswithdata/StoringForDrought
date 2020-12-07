rm(list=ls())

library(ggplot2)
library(wesanderson)
library(ggpubr)
library(RColorBrewer)


## Figure 1 - several examples of r

source("code/parameters.R")
source("code/models/water_model.R")

param <- param_list_drake_W0_1000()
param2 <- param_list_drake_W0_3000()


deltat = 0.1
# Initialise a vector of possible ts values and the length of this
# vector
t_seq <- seq(from=0, to=param$Tend, by=deltat)
simlength <- length(t_seq)

kf <- seq(from=0, to=1, by=0.01)

ts_out_1000 <- vector()
tcrit_out_1000 <- vector()
optim_out_1000 <- vector()
M_out_kf_1000 <- vector()
S_out_kf_1000 <- vector()

ts_out_3000 <- vector()
tcrit_out_3000 <- vector()
optim_out_3000 <- vector()
M_out_kf_3000 <- vector()
S_out_kf_3000 <- vector()

j = 1
for (f in kf){
  S_out <- vector(length=simlength)
  M_out <- vector(length=simlength)
  tcrit <- vector(length=simlength)
  S_out_3 <- vector(length=simlength)
  M_out_3 <- vector(length=simlength)
  tcrit_3 <- vector(length=simlength)
  
  k=1
  for (i in t_seq){
    # Find the output values of simulation for the specific ts value
    
    out <- simple_water_model_sim_time_breaks(param, i, deltat)
    
    # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
    S_out[k] <- out$S[simlength]
    M_out[k] <- out$M[simlength]
    tcrit[k] <- t_seq[which.min(out$W)]
    
    if(sum(out$S<0) > 0){
      S_out[k] <- NA
      M_out[k] <- NA
      tcrit[k] <- NA
    }
    
    out <- simple_water_model_sim_time_breaks(param2, i, deltat)
    
    # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
    S_out_3[k] <- out$S[simlength]
    M_out_3[k] <- out$M[simlength]
    tcrit_3[k] <- t_seq[which.min(out$W)]
    
    if(sum(out$S<0) > 0){
      S_out_3[k] <- NA
      M_out_3[k] <- NA
      tcrit_3[k] <- NA
    }
    
    k = k+1
  }
  
  ms_sum <- colSums(rbind (f * M_out, (1 - f) * S_out))
  ts_ind <- which.max(ms_sum)
  if(length(ts_ind) == 0){
    ts_out_1000[j] = NA
    tcrit_out_1000[j] = NA
    M_out_kf_1000[j] = NA
    optim_out_1000[j] = NA
    S_out_kf_1000[j] = NA
  }
  else{
    ts_out_1000[j] <- t_seq[ts_ind]
    tcrit_out_1000[j] = tcrit[ts_ind]
    M_out_kf_1000[j] = M_out[ts_ind]
    optim_out_1000[j] = ms_sum[ts_ind]
    S_out_kf_1000[j] = S_out[ts_ind]
  }
  
  ms_sum <- colSums(rbind (f * M_out_3, (1 - f) * S_out_3))
  ts_ind <- which.max(ms_sum)
  if(length(ts_ind) == 0){
    ts_out_3000[j] = NA
    tcrit_out_3000[j] = NA
    M_out_kf_3000[j] = NA
    optim_out_3000[j] = NA
    S_out_kf_3000[j] = NA
  }
  else{
    ts_out_3000[j] <- t_seq[ts_ind]
    tcrit_out_3000[j] = tcrit_3[ts_ind]
    M_out_kf_3000[j] = M_out_3[ts_ind]
    optim_out_3000[j] = ms_sum[ts_ind]
    S_out_kf_3000[j] = S_out_3[ts_ind]
  }
  
  j = j+1
}


df_1000 <- data.frame(kf = rep(kf, times = 4),
                 t = c(rep(ts_out_1000, times=2),rep(tcrit_out_1000, times = 2)),
                 t_type = rep(c("ts", "tcrit"), each=2*length(kf)),
                 C_pool = rep(c(M_out_kf_1000, S_out_kf_1000), times=4),
                 C_pool_type = rep(rep(c("Biomass", "Storage"), each=length(kf)), times=2) ,
                 optim_val = rep(optim_out_1000, times=4),
                 w0 = as.factor(rep(1000, times=4*length(kf)))
                 )
df_3000 <- data.frame(kf = rep(kf, times = 4),
                      t = c(rep(ts_out_3000, times=2),rep(tcrit_out_3000, times = 2)),
                      t_type = rep(c("ts", "tcrit"), each=2*length(kf)),
                      C_pool = rep(c(M_out_kf_3000, S_out_kf_3000), times=4),
                      C_pool_type = rep(rep(c("Biomass", "Storage"), each=length(kf)), times=2) ,
                      optim_val = rep(optim_out_3000, times=4),
                      w0 = as.factor(rep(3000, times=4*length(kf)))
)

df <- rbind(df_1000, df_3000)

save(df, param, file="opt_ts_kf_w0_fig5_dat.RData")


p <- ggplot(df, aes(x = kf, y = t)) +
  geom_line(aes(linetype=t_type, colour=w0)) +
  scale_y_continuous("Time (day)") +
  scale_x_continuous(expression(paste("Fitness Proxy Parameter ", k[f]))) +
  scale_color_discrete(expression(paste("Initial Water Content kg", H[2], "O")))+
  scale_linetype_manual ("Timepoint (day)", labels=c(expression(paste(t[crit])),expression(paste( t[s]))), 
                         values = c(2,1)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p
