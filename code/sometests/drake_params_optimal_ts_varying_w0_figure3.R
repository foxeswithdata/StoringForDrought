rm(list=ls())

library(ggplot2)
library(wesanderson)
library(ggpubr)
library(RColorBrewer)


## Figure 1 - several examples of r

source("code/parameters.R")
source("code/models/water_model.R")

param <- param_list_drake_W0_1000()


W0 = seq(from=500, to=7000, by=100);

deltat = 0.1
# Initialise a vector of possible ts values and the length of this
# vector
t_seq <- seq(from=0, to=param$Tend, by=deltat)
simlength <- length(t_seq)

ts_maxm <- vector();
ts_maxs <- vector();
ts_maxms <- vector();

tcrit_maxm <- vector();
tcrit_maxs <- vector();
tcrit_maxms <- vector();

MT_maxm <- vector();
MT_maxs <- vector();
MT_maxms <- vector();

ST_maxm <- vector();
ST_maxs <- vector();
ST_maxms <- vector();

MST_maxm <- vector();
MST_maxs <- vector();
MST_maxms <- vector();

WT_maxm <- vector();
WT_maxs <- vector();
WT_maxms <- vector();



j = 1
for (w0 in W0){
  print(w0)
  S_out <- vector(length=simlength)
  M_out <- vector(length=simlength)
  W_out <- vector(length=simlength)
  tcrit_out <- vector(length=simlength)
  
  k=1
  for (i in t_seq){
    # Find the output values of simulation for the specific ts value
    param_2 <- param;
    param_2$W0 <- w0;
    param_2$x0 <- c(param_2$M0, param_2$S0, param_2$W0)
    
    out <- simple_water_model_sim_time_breaks(param_2, i, deltat)
    
    # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
    S_out[k] <- out$S[simlength]
    M_out[k] <- out$M[simlength]
    W_out[k] <- out$W[simlength]
    tcrit_out[k] <- t_seq[which.min(out$W)]
    
    if(sum(out$S<0) > 0){
      S_out[k] <- NA
      M_out[k] <- NA
      W_out[k] <- NA
      tcrit_out[k] <- NA
    }
    
    k = k+1
  }
  
  kf = 1
  ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
  ts_ind <- which.max(ms_sum)
  if(length(ts_ind) == 0){
    ts_maxm[j] = NA
    tcrit_maxm[j] = NA
    MT_maxm[j] = NA
    MST_maxm[j] = NA
    ST_maxm[j] = NA
    WT_maxm[j] = NA
  }
  else{
    ts_maxm[j] <- t_seq[ts_ind]
    tcrit_maxm[j] = tcrit_out[ts_ind]
    MT_maxm[j] = M_out[ts_ind]
    MST_maxm [j] = M_out[ts_ind] + S_out[ts_ind]
    ST_maxm[j] = S_out[ts_ind]
    WT_maxm[j] = W_out[ts_ind]
  }
  
  kf = 0.5
  ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
  ts_ind <- which.max(ms_sum)
  if(length(ts_ind) == 0){
    ts_maxms[j] = NA
    tcrit_maxms[j] = NA
    MT_maxms[j] = NA
    MST_maxms[j] = NA
    ST_maxms[j] = NA
    WT_maxms[j] = NA
  }
  else{
    ts_maxms[j] <- t_seq[ts_ind]
    tcrit_maxms[j] = tcrit_out[ts_ind]
    MT_maxms[j] = M_out[ts_ind]
    MST_maxms[j] = M_out[ts_ind] + S_out[ts_ind]
    ST_maxms[j] = S_out[ts_ind]
    WT_maxms[j] = W_out[ts_ind]
  }
  
  kf = 0
  ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
  ts_ind <- which.max(ms_sum)
  if(length(ts_ind) == 0){
    ts_maxs[j] = NA
    tcrit_maxs[j] = NA
    MT_maxs[j] = NA
    MST_maxs[j] = NA
    ST_maxs[j] = NA
    WT_maxs[j] = NA
  }
  else{
    ts_maxs[j] <- t_seq[ts_ind]
    tcrit_maxs[j] = tcrit_out[ts_ind]
    MT_maxs[j] = M_out[ts_ind]
    MST_maxs[j] = M_out[ts_ind] + S_out[ts_ind]
    ST_maxs[j] = S_out[ts_ind]
    WT_maxs[j] = W_out[ts_ind]
  }
  
  j = j+1
}


W0_length <- length(W0)

resdf <- data.frame(t_val = rep(c(ts_maxm, ts_maxs, ts_maxms, tcrit_maxm, tcrit_maxs, tcrit_maxms),2),
                    type = rep(rep(c("ts", "tcrit"), each=W0_length*3),2),
                    goal = rep(rep(c("MaxM", "MaxS", "MaxMS"), each=W0_length), times=4),
                    w0 = rep(W0, 12),
                    carbon_pool_val = c(rep(c(MT_maxm, MT_maxs, MT_maxms),2), 
                                        rep(c(ST_maxm, ST_maxs, ST_maxms),2)),
                    carbon_pool_perc = c(rep(c(MT_maxm/(MT_maxm+ST_maxm), 
                                               MT_maxs/(MT_maxs+ST_maxs), 
                                               MT_maxms/(MT_maxms+ST_maxms)),2), 
                                         rep(c(ST_maxm/(MT_maxm+ST_maxm), 
                                               ST_maxs/(MT_maxs+ST_maxs), 
                                               ST_maxms/(MT_maxms+ST_maxms)),2)),
                    carbon_pool_type = as.factor(rep(c("Biomass", "Storage"), each=W0_length*6)),
                    water_pool = rep(c(WT_maxm, WT_maxs, WT_maxms), times = 4)
                    )

save(resdf, W0, t_seq, param, file="optimal_traj_data.RData")


rm(list=ls())
load("optimal_traj_data.RData")


resdf_tsonly <- subset(resdf, type=="ts")
resdf_tsonly <- subset(resdf_tsonly, goal%in%c("MaxM", "MaxS"))



p1 <- ggplot(data=resdf_tsonly)+
  geom_line(aes(x = w0, y = t_val, colour = goal)) +
  scale_y_continuous(expression(paste("Time (day)"))) +
  scale_x_continuous(expression(paste("Initial Water Content ", W[0] ," (kg", H[2], "O)"))) +
  scale_color_manual("Plant Goal", values = brewer.pal(n = 3, name = "Dark2"),labels = c("Maximising \n Biomass \n", "Maximising \n Storage \n"))+
  # scale_linetype_manual ("Timepoint (day)", labels=c(expression(paste(t[crit])),expression(paste( t[s]))),
                         # values = c(2,1)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p1


p2 <- ggplot(data=resdf)+
  geom_line(aes(x = w0, y = carbon_pool_val, colour = goal, linetype=carbon_pool_type)) +
  scale_y_continuous(expression(paste("Final Carbon Pool Size (gC)"))) +
  scale_x_continuous(expression(paste("Initial Water Content ", W[0] ," (kg", H[2], "O)"))) +
  scale_color_manual("Plant Goal", values = brewer.pal(n = 3, name = "Dark2"),labels = c("Maximising \n Biomass \n", "Maximising \n Biomass \n and Storage\n", "Maximising \n Storage \n"))+
  scale_linetype_manual ("Carbon Pool", labels=c("Biomass", "Storage"), 
                         values = c(1,2)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p2


ggarrange(p1,p2, labels = c("A", "B"), nrow = 2)

p3 <- ggplot(data=resdf)+
  geom_line(aes(x = w0, y = carbon_pool_perc, colour = goal, linetype=carbon_pool_type)) +
  scale_y_continuous(expression(paste("Final Carbon Pool Size (gC)"))) +
  scale_x_continuous(expression(paste("Initial Water Content ", W[0] ," (kg", H[2], "O)"))) +
  scale_color_manual("Plant Goal", values = brewer.pal(n = 3, name = "Dark2"),labels = c("Maximising \n Biomass \n", "Maximising \n Biomass \n and Storage\n", "Maximising \n Storage \n"))+
  scale_linetype_manual ("Carbon Pool", labels=c("Biomass", "Storage"), 
                         values = c(1,2)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p3










####### Find difference between optimal allocation for W0 4100 and only storage allocation to W0 4100
####### Specifically for the storage pool

results = unique(subset(resdf, w0==4100 & carbon_pool_type == "Storage" & goal == "MaxS"))
opt_store <- unique(results$carbon_pool_val)


param_2 <- param;
param_2$W0 <- 4100;
param_2$x0 <- c(param_2$M0, param_2$S0, param_2$W0)

out <- simple_water_model_sim_time_breaks(param_2, 0, deltat)
onlyS_store <- out$S[simlength]

(opt_store - onlyS_store)/opt_store
(opt_store - onlyS_store)/onlyS_store




