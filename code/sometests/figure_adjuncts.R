rm(list=ls())

library(ggplot2)
library(wesanderson)
library(ggpubr)
library(gridExtra)


## Figure 1 - several examples of r

source("code/parameters.R")
source("code/models/water_model.R")

param <- param_list_drake_W0_3000()
param$W0 <- 4000
deltat = 0.1
# Initialise a vector of possible ts values and the length of this
# vector
t_seq <- seq(from=0, to=param$Tend, by=deltat)
simlength <- length(t_seq)

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

kf = 0
ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
ts_ind <- which.max(ms_sum)
ts_out = t_seq[ts_ind]
ts_out_maxS <- ts_out

out <- simple_water_model_sim_time_breaks(param, ts_out_maxS, deltat)

tcrit_out_maxS <- t_seq[which.min(out$W)]

mu = c(0.01)
out <- simple_water_model_adjuncts_sim_time_breaks(param, c(0,1,0), 30, 135, m, deltat)

adjunct_df <- data.frame(lambda = c(out$lamM, out$lamS, out$lamW),
                         t = rep(t_seq, times=3),
                         lambda_type = rep(c("M", "S","W"), each=length(t_seq)))


p <- ggplot(adjunct_df, aes(t, lambda, col=lambda_type))+
  geom_line(linetype="dashed") + 
  geom_vline(aes(xintercept=30), color="red", size = 0.25,linetype="dotted") +
  geom_vline(aes(xintercept=135), color="blue", size = 0.25,linetype="dotted") +
  # geom_text(aes(x=tcrit - 15, label=paste("tcrit = ", tcrit), y=6000), 
  #           mean_tcrit, colour="blue", angle=90, text=element_text(size=8), vjust = 1,) +
  scale_y_continuous(expression(paste("Adjunct Value ", lambda)))+
  scale_x_continuous("Time (day)", breaks=c(30, 135, 150), 
                     labels = c(expression(t[s]), expression(t[crit]), "T")) +
  scale_color_hue("Equivalent pool", guide=guide_legend(order=1)) +
  theme_bw() + 
  theme(axis.text.y = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p  
  
  
  