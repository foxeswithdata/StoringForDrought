rm(list=ls())

library(ggplot2)
library(wesanderson)
library(ggpubr)
library(lemon)
library(cowplot)
library(plyr)
library(reshape2)

source("code/parameters.R")
source("code/models/water_model.R")

param <- param_list_drake_W0_1000()
param_2 <- param_list_drake_W0_3000()

deltat = 0.1
t_seq <- seq(from=0, to=param$Tend, by=deltat)
simlength <- length(t_seq)

kf_1000 <- seq(from=0, to=1, by=0.20)
ts_1000 <- seq(from=0, to=30, by=0.5)

kf_3000 <- seq(from=0, to=1, by=0.20)
ts_3000 <- seq(from=0, to=100, by=1)

C_out_1000 <- matrix(nrow=length(kf_1000), ncol=length(ts_1000))
rownames(C_out_1000) <- kf_1000
colnames(C_out_1000) <- ts_1000


C_out_3000 <- matrix(nrow=length(kf_3000), ncol=length(ts_3000))
rownames(C_out_3000) <- kf_3000
colnames(C_out_3000) <- ts_3000

i = 1
for(t in ts_1000){
  out <- simple_water_model_sim_time_breaks(param, t, deltat)
  j = 1
  for(f in kf_1000){
    C_out_1000[j,i] = out$M[simlength]*f + out$S[simlength]*(1-f)
    j=j+1
  }
  i=i+1
}

i = 1
for(t in ts_3000){
  out <- simple_water_model_sim_time_breaks(param_2, t, deltat)
  j = 1
  for(f in kf_3000){
    C_out_3000[j,i] = out$M[simlength]*f + out$S[simlength]*(1-f)
    j=j+1
  }
  i=i+1
}


tcrit_1000 <- vector()
tcrit_3000 <- vector()
i = 1
for(t in seq(from=0, to=param$Tend)){
  out_1 <- simple_water_model_sim_time_breaks(param, t, deltat)
  out_2 <- simple_water_model_sim_time_breaks(param_2, t, deltat)
  tcrit_1000[i] = t_seq[which.min(out_1$W)]
  tcrit_3000[i] = t_seq[which.min(out_2$W)]
  i = i+1
}




C_out_1000_df <- melt(C_out_1000, value.name = "C_out")
names(C_out_1000_df) <- c("kf", "ts", "C_out")
C_out_1000_df$kf <- as.factor(C_out_1000_df$kf)
C_out_1000_df$W0 <- as.factor(rep(param$W0, times=length(C_out_1000_df$kf)))
ts_max_1000 <- ts_1000[which.min(abs(C_out_1000_df$C_out[C_out_1000_df$kf==0]))]

C_out_3000_df <- melt(C_out_3000, value.name = "C_out")
names(C_out_3000_df) <- c("kf", "ts", "C_out")
C_out_3000_df$kf <- as.factor(C_out_3000_df$kf)
C_out_3000_df$W0 <- as.factor(rep(param_2$W0, times=length(C_out_3000_df$kf)))
ts_max_3000 <- ts_3000[which.min(abs(C_out_3000_df$C_out[C_out_3000_df$kf==0]))]

C_out_df <- rbind(C_out_1000_df, C_out_3000_df)

tcrit_df <- data.frame(ts = c(0:150, 0:150), tcrit = c(tcrit_1000, tcrit_3000),
                       W0=as.factor(c(rep(param$W0, times=length(tcrit_1000)),
                            rep(param_2$W0,times=length(tcrit_3000)))))





p_1 <- ggplot(C_out_1000_df, aes(x = ts, y = C_out, colour = kf, linetype=kf))+
  geom_line() + 
  scale_y_continuous(name = "Terminating Pool Size (gC)") +
  scale_x_continuous(name = expression(paste("Time of switch ",t[s]," (d)"))) +
  scale_color_hue(name = expression(paste("Strategy ",k[f]))) +
  scale_linetype_manual(values=c("solid", "twodash","twodash","twodash","twodash","solid"))+
  geom_hline(yintercept=0, color="#000000", size = 0.25,linetype="dashed") + 
  geom_vline(xintercept=ts_max_1000, color="#000000", size = 0.25,linetype="dotted") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position="none")
p_1

p_2 <- ggplot(C_out_3000_df, aes(x = ts, y = C_out, colour = kf, linetype=kf))+
  geom_line() + 
  scale_y_continuous(name = "Terminating Pool Size (gC)") +
  scale_x_continuous(name = expression(paste("Time of switch ",t[s]," (d)"))) +
  scale_color_hue(name = expression(paste("Strategy ",k[f]))) +
  scale_linetype_manual(values=c("solid", "twodash","twodash","twodash","twodash","solid"))+
  geom_hline(yintercept=0, color="#000000", size = 0.25,linetype="dashed") + 
  geom_vline(xintercept=ts_max_3000, color="#000000", size = 0.25,linetype="dotted") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p_2

p_3 <- ggplot(tcrit_df, aes(x=ts, y=tcrit, linetype=W0))+
  geom_line() + 
  scale_y_continuous(name = expression(paste("Time of water depletion ", t[crit], " (d)"))) + 
  scale_x_continuous(name = expression(paste("Time of switch ",t[s]," (d)"))) +
  scale_linetype(name=expression(W[0])) + 
  geom_hline(yintercept=0, color="#000000", size = 0.25,linetype="dashed") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p_3

p <- ggarrange(ggarrange(p_1,p_2, labels=c("A","B")),p_3, nrow=2, labels=c("","C"))
p


