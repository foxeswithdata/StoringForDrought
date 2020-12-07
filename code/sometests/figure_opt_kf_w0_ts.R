rm(list=ls())

library(ggplot2)
library(wesanderson)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(lattice)


## Figure 1 - several examples of r

source("code/parameters.R")
source("code/models/water_model.R")

param <- param_list_drake_W0_1000()


W0 = seq(from=500, to=7000, by=100);
kf = seq(from=0, to=1, by=0.01)

deltat = 0.1
# Initialise a vector of possible ts values and the length of this
# vector
t_seq <- seq(from=0, to=param$Tend, by=deltat)
simlength <- length(t_seq)

ts_out <- matrix(nrow = length(kf), ncol=length(W0))
ST_out <- matrix(nrow = length(kf), ncol=length(W0))
MT_out <- matrix(nrow = length(kf), ncol=length(W0))
tcrit_out <- matrix(nrow = length(kf), ncol=length(W0))
WT_out <- matrix(nrow = length(kf), ncol=length(W0))

j = 1
for (w0 in W0){
  print(w0)
  S_out <- vector(length=simlength)
  M_out <- vector(length=simlength)
  W_out <- vector(length=simlength)
  tcrit <- vector(length=simlength)

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
    tcrit[k] <- t_seq[which.min(out$W)]

    if(sum(out$S<0) > 0){
      S_out[k] <- NA
      M_out[k] <- NA
      W_out[k] <- NA
      tcrit[k] <- NA
    }
    k = k+1
  }

  k=1
  for (i in kf){
    ms_sum <- colSums(rbind (i * M_out, (1 - i) * S_out))
    ts_ind <- which.max(ms_sum)
    if(length(ts_ind) == 0){
      ts_out[k,j] = NA
      ST_out[k,j] = NA
      MT_out[k,j] = NA
      tcrit_out[k,j] = NA
      WT_out[k,j] = NA
    }
    else{
      ts_out[k,j] <- t_seq[ts_ind]
      ST_out[k,j] = S_out[ts_ind]
      MT_out[k,j] = M_out[ts_ind]
      tcrit_out[k,j] = tcrit[ts_ind]
      WT_out[k,j] = W_out[ts_ind]
    }
    k = k + 1
  }
  j = j+1
}


save(ts_out, kf, W0, t_seq, param, file="optimal_ts_kf_data.RData")

load("optimal_ts_kf_data.RData") 

xycoords <- expand.grid(x=kf,y=W0)
rownames(ts_out) <- kf
colnames(ts_out) <- W0
z <- melt(ts_out)


kf_W0_df <- data.frame(ts = z$value, 
                 ST = melt(ST_out)$value,
                 MT = melt(MT_out)$value, 
                 tcrit = melt(tcrit_out)$value, 
                 WT = melt(WT_out)$value,
                 kf = z$Var1, W0 = z$Var2)

save(kf_W0_df, t_seq, param, file="optimal_ts_full_kf_data.RData")

load(file="optimal_ts_full_kf_data.RData")

p1 <- ggplot(data=kf_W0_df) +
  geom_tile(aes(x = W0, y = kf, fill = ts)) + 
  labs(x = expression(paste("Initial Water Content ", W[0], " (kg", H[2], "O)")),
       y = expression(paste("Proportion of Biomass in Goal Function ", k[f]))) +
  scale_fill_distiller(expression(paste("Time of", " switch ", t[s], " (day)")), palette = "Greens") + 
  theme_bw() +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
 
p2 <- ggplot(data=kf_W0_df) +
  geom_tile(aes(x = W0, y = kf, fill = tcrit)) + 
  labs(x = expression(paste("Initial Water Content ", W[0], " (kg", H[2], "O)")),
       y = expression(paste("Proportion of Biomass in Goal Function ", k[f]))) +
  scale_fill_distiller(expression(paste("Time of", " water depletion ", t[crit], " (day)")), palette = "Blues") + 
  theme_bw() +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p3 <- ggplot(data=kf_W0_df) +
  geom_tile(aes(x = W0, y = kf, fill = ST)) + 
  labs(x = expression(paste("Initial Water Content ", W[0], " (kg", H[2], "O)")),
       y = expression(paste("Proportion of Biomass in Goal Function ", k[f]))) +
  scale_fill_distiller(expression(paste("Final Storage \n Size kgC")), palette = "Reds") + 
  theme_bw() +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p4 <- ggplot(data=kf_W0_df) +
  geom_tile(aes(x = W0, y = kf, fill = MT)) + 
  labs(x = expression(paste("Initial Water Content ", W[0], " (kg", H[2], "O)")),
       y = expression(paste("Proportion of Biomass in Goal Function ", k[f]))) +
  scale_fill_distiller(expression(paste("Final Biomass \n Size kgC")), palette = "Oranges") + 
  theme_bw() +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


ggarrange(p1,p2,p3,p4, ncol=2, nrow=2)



# From the lattice package
wireframe(ts ~ W0*kf, data = kf_W0_df,
          xlab = expression(paste(W[0])),
          ylab = expression(paste(k[f])),
          zlab = expression(paste(t[s])),
          drape = TRUE,
          colorkey = TRUE
          # screen = list(z = -60, x = -60)
)



W0_list <- unique(kf_W0_df$W0)
kf_out_list <- vector(length=length(W0_list))

j = 1
for (w0 in W0_list) {
  finder <- subset(kf_W0_df, W0==w0, select=c(kf,ts))
  ind <- which.max(diff(finder$ts))
  if(length(ind) == 0) {
    kf_out_list[j] <- NA
  }
  else{
    kf_out_list[j] <- finder$kf[ind]
  }
  
  j = j+1
}

df_extra <- data.frame(W0=W0_list, kf=kf_out_list)

p <- ggplot(df_extra, aes(x = W0, y = kf)) + 
  geom_line() + 
  scale_x_continuous(expression(paste("Initial Water Availability ", W[0], " kg", H[2], "O"))) +
  scale_y_continuous(expression(paste("Fitness Proxy Parameter ", k[f]))) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p



