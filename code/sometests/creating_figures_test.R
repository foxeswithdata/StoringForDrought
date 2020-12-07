library(ggplot2)
library(wesanderson)
library(ggpubr)


## Figure 1 - several examples of r

source("code/parameters.R")
source("code/models/water_model.R")

param <- param_list_water_short2()

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

kf = 1
ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
ts_ind <- which.max(ms_sum)
ts_out = t_seq[ts_ind]
ts_out_maxM <- ts_out

out <- simple_water_model_sim_time_breaks(param, ts_out_maxM, deltat)

out_maxM <- data.frame(S= out$S, M=out$M, W=out$W, t=t_seq)
tcrit_out_maxM <- t_seq[which.min(out$W)]

kf = 0
ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
ts_ind <- which.max(ms_sum)
ts_out = t_seq[ts_ind]
ts_out_maxS <- ts_out


out <- simple_water_model_sim_time_breaks(param, ts_out_maxS, deltat)
tcrit_out_maxS <- t_seq[which.min(out$W)]


out_maxS <- data.frame(S= out$S, M=out$M, W=out$W, t=t_seq)

coeff <- 1.5
colors <- wes_palette("Darjeeling1", n = 3)

# Plotting MaxS

p1 <- ggplot(out_maxS, aes(x=t)) +
  
  geom_line( aes(y=S, colour="Storage"), color=colors[3]) + 
  geom_line( aes(y=M, colour="Biomass"), color=colors[1]) +
  geom_line( aes(y=W / coeff), color=colors[2], size=1.5, alpha=0.5) + # Divide by 10 to get the same range than the temperature
  xlab("time (day)") +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Pool Size (gC)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name=expression(paste("Water Availability (kg", H[2], "O)")))
  ) +
  geom_vline(xintercept=ts_out_maxS, colour="grey") +
  geom_text(aes(x=ts_out_maxS, label="allocation switch", y=50), colour="grey", angle=90, text=element_text(size=11), vjust = 1,) +
  geom_vline(xintercept=tcrit_out_maxS, colour="grey") +
  geom_text(aes(x=tcrit_out_maxS, label="Water lost", y=50), colour="grey", angle=90, text=element_text(size=11), vjust = 1,) +
  scale_colour_discrete(name  ="Pool",
                          breaks=c("Storage", "Biomass"),
                          labels=c("Storage", "Biomass")) +
    theme_minimal()
p

# Plotting MaxM

p2 <- ggplot(out_maxM, aes(x=t)) +
  
  geom_line( aes(y=S, colour="Storage"), color=colors[3]) + 
  geom_line( aes(y=M, colour="Biomass"), color=colors[1]) +
  geom_line( aes(y=W / coeff), color=colors[2], size=1.5, alpha=0.5) + # Divide by 10 to get the same range than the temperature
  xlab("time (day)") +
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Pool Size (gC)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name=expression(paste("Water Availability (kg", H[2], "O)")))
  ) +
  geom_vline(xintercept=ts_out_maxM, colour="grey") +
  geom_text(aes(x=ts_out_maxM, label="allocation switch", y=50), colour="grey", angle=90, text=element_text(size=11), vjust = 1,) +
  geom_vline(xintercept=tcrit_out_maxM, colour="grey") +
  geom_text(aes(x=tcrit_out_maxM, label="Water lost", y=50), colour="grey", angle=90, text=element_text(size=11), vjust = 1,) +
  scale_colour_discrete(name  ="Pool",
                        breaks=c("Storage", "Biomass"),
                        labels=c("Storage", "Biomass")) +
  theme_minimal()
p












param$W0 = 150
deltat = 0.1
# Initialise a vector of possible ts values and the length of this
# vector
t_seq <- seq(from=0, to=param$Tend, by=deltat)
simlength <- length(t_seq)

S_out_2 <- vector(length=simlength)
M_out_2 <- vector(length=simlength)

k=1
for (i in t_seq){
  # Find the output values of simulation for the specific ts value
  
  out <- simple_water_model_sim_time_breaks(param, i, deltat)
  
  # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
  S_out_2[k] <- out$S[simlength]
  M_out_2[k] <- out$M[simlength]
  
  if(sum(out$S<0) > 0){
    S_out_2[k] <- NA
    M_out_2[k] <- NA
  }
  k = k+1
}

kf = 1
ms_sum <- colSums(rbind (kf * M_out_2, (1 - kf) * S_out_2))
ts_ind <- which.max(ms_sum)
ts_out = t_seq[ts_ind]
ts_out_maxM_2 <- ts_out

out <- simple_water_model_sim_time_breaks(param, ts_out_maxM_2, deltat)
tcrit_out_maxM_2 <- t_seq[which.min(out$W)]

out_maxM_2 <- data.frame(S= out$S, M=out$M, W=out$W, t=t_seq)

kf = 0
ms_sum <- colSums(rbind (kf * M_out_2, (1 - kf) * S_out_2))
ts_ind <- which.max(ms_sum)
ts_out = t_seq[ts_ind]
ts_out_maxS_2 <- ts_out

out <- simple_water_model_sim_time_breaks(param, ts_out_maxS_2, deltat)
tcrit_out_maxS_2 <- t_seq[which.min(out$W)]
out_maxS_2 <- data.frame(S= out$S, M=out$M, W=out$W, t=t_seq)





coeff <- 1.5
colors <- wes_palette("Darjeeling1", n = 3)

# Plotting MaxS

p3 <- ggplot(out_maxS_2, aes(x=t)) +
  
  geom_line( aes(y=S, colour="Storage"), color=colors[3]) + 
  geom_line( aes(y=M, colour="Biomass"), color=colors[1]) +
  geom_line( aes(y=W / coeff), color=colors[2], size=1.5, alpha=0.5) + # Divide by 10 to get the same range than the temperature
  xlab("time (day)") +
  scale_y_continuous(
    
    # Features of the first axis
    name = "Pool Size (gC)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name=expression(paste("Water Availability (kg", H[2], "O)")))
  ) +
  geom_vline(xintercept=ts_out_maxS_2, colour="grey") +
  geom_text(aes(x=ts_out_maxS_2, label="allocation switch", y=50), colour="grey", angle=90, text=element_text(size=11), vjust = 1,) +
  geom_vline(xintercept=tcrit_out_maxS_2, colour="grey") +
  geom_text(aes(x=tcrit_out_maxS_2, label="Water lost", y=50), colour="grey", angle=90, text=element_text(size=11), vjust = 1,) +
  scale_colour_discrete(name  ="Pool",
                        breaks=c("Storage", "Biomass"),
                        labels=c("Storage", "Biomass")) +
  theme_minimal()
p

# Plotting MaxM

p4 <- ggplot(out_maxM_2, aes(x=t)) +
  
  geom_line( aes(y=S, colour="Storage"), color=colors[3]) + 
  geom_line( aes(y=M, colour="Biomass"), color=colors[1]) +
  geom_line( aes(y=W / coeff), color=colors[2], size=1.5, alpha=0.5) + # Divide by 10 to get the same range than the temperature
  xlab("time (day)") +
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Pool Size (gC)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name=expression(paste("Water Availability (kg", H[2], "O)")), limits=c(0,220))
  ) +
  geom_vline(xintercept=ts_out_maxM_2, colour="grey") +
  geom_text(aes(x=ts_out_maxM_2, label="allocation switch", y=50), colour="grey", angle=90, text=element_text(size=11), vjust = 1,) +
  geom_vline(xintercept=tcrit_out_maxM_2, colour="grey") +
  geom_text(aes(x=tcrit_out_maxM_2, label="Water lost", y=50), colour="grey", angle=90, text=element_text(size=11), vjust = 1,) +
  scale_colour_discrete(name  ="Pool",
                        breaks=c("Storage", "Biomass"),
                        labels=c("Storage", "Biomass")) +
  theme_minimal()
p

ggarrange(p1, p2, p3, p4, 
          labels = c("A", "B", "C", "D"),
          ncol = 4, nrow = 1)





















# Plotting MaxM - changing axis

p2 <- ggplot(out_maxM, aes(x=t)) +
  
  geom_line( aes(y=S, colour="Storage"), color=colors[3]) + 
  geom_line( aes(y=M, colour="Biomass"), color=colors[1]) +
  geom_line( aes(y=W / coeff), color=colors[2], size=1.5, alpha=0.5) + # Divide by 10 to get the same range than the temperature
  xlab("time (day)") +
  geom_vline(xintercept=ts_out_maxM, colour="grey", alpha=0.8, linetype=2) +
  geom_vline(xintercept=tcrit_out_maxM, colour="grey", alpha=0.8, linetype=2) +
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Pool Size (gC)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name=expression(paste("Water Availability (kg", H[2], "O)")))
  ) +
  scale_colour_discrete(name  ="Pool",
                        breaks=c("Storage", "Biomass"),
                        labels=c("Storage", "Biomass")) +
  theme_minimal()
p2 + scale_x_continuous(breaks=c(ts_out_maxM, tcrit_out_maxM, param$Tend),
                      labels=c("ts","tcrit", "T"))















# Plotting MaxM - changing background

rect <- data.frame(xmin = c(0, ts_out_maxM, tcrit_out_maxM), xmax = c(ts_out_maxM, tcrit_out_maxM, param$Tend), col = c("growth", "storage", "stress"))

p2 <- ggplot() +
  geom_rect(data = rect, aes(xmin = xmin, xmax = xmax, fill = col, 
                ymin = -Inf, ymax = Inf), alpha = 0.4) +
  geom_line( data = out_maxM, aes(x = t, y=S, colour="Storage"), color=colors[3]) + 
  geom_line( data = out_maxM, aes(x = t, y=M, colour="Biomass"), color=colors[1]) +
  geom_line( data = out_maxM, aes(x = t, y=W / coeff), color=colors[2], size=1.5, alpha=0.5) + # Divide by 10 to get the same range than the temperature
  xlab("time (day)") +
  geom_vline(xintercept=ts_out_maxM, colour="grey", alpha=0.8, linetype=2) +
  geom_vline(xintercept=tcrit_out_maxM, colour="grey", alpha=0.8, linetype=2) +
  
  scale_y_continuous(
    
    # Features of the first axis
    name = "Pool Size (gC)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name=expression(paste("Water Availability (kg", H[2], "O)")))
  ) +
  scale_colour_discrete(name  ="Pool",
                        breaks=c("Storage", "Biomass"),
                        labels=c("Storage", "Biomass")) +
  theme_minimal()
# p2 + scale_x_continuous(breaks=c(ts_out_maxM, tcrit_out_maxM, param$Tend),
                        # labels=c("ts","tcrit", "T"))

p2 




