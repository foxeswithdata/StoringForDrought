rm(list=ls())

library(ggplot2)
library(wesanderson)
library(ggpubr)
library(gridExtra)


## Figure 1 - several examples of r

source("code/parameters.R")
source("code/models/water_model.R")

param <- param_list_drake_W0_1000()
param_2 <- param_list_drake_W0_3000()
deltat = 0.1
# Initialise a vector of possible ts values and the length of this
# vector
t_seq <- seq(from=0, to=param$Tend, by=deltat)
simlength <- length(t_seq)

S_out <- vector(length=simlength)
M_out <- vector(length=simlength)
S_out_2 <- vector(length=simlength)
M_out_2 <- vector(length=simlength)

k=1
for (i in t_seq){
  # Find the output values of simulation for the specific ts value
  
  out <- simple_water_model_sim_time_breaks(param, i, deltat)
  out_2 <- simple_water_model_sim_time_breaks(param_2, i, deltat)
  
  # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
  S_out[k] <- out$S[simlength]
  M_out[k] <- out$M[simlength]
  
  S_out_2[k] <- out_2$S[simlength]
  M_out_2[k] <- out_2$M[simlength]
  
  if(sum(out$S<0) > 0){
    S_out[k] <- NA
    M_out[k] <- NA
  }
  
  if(sum(out_2$S<0) > 0){
    S_out_2[k] <- NA
    M_out_2[k] <- NA
  }
  
  k = k+1
}

kf = 1
ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
ts_ind <- which.max(ms_sum)
ts_out = t_seq[ts_ind]
ts_out_maxM <- ts_out

out <- simple_water_model_sim_time_breaks(param, ts_out_maxM, deltat)

tcrit_out_maxM <- t_seq[which.min(out$W)]

out_maxM_1000 <- data.frame(carbon_pool_type = as.factor(rep(c("Storage", "Biomass"), each = simlength)),
                  carbon_pool_val = c(out$S, out$M),
                  water_pool_val = rep(out$W, times = 2),
                  water_pool_type = as.factor(rep("Soil Water Available", times = simlength * 2)),
                  initial_water = rep(param$W0, simlength * 2),
                  ts = (rep(ts_out_maxM, times = simlength * 2)),
                  tcrit = (rep(tcrit_out_maxM, times = simlength * 2)),
                  goal = rep("Biomass", simlength * 2),
                  time = rep(t_seq, times = 2))

kf = 0
ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
ts_ind <- which.max(ms_sum)
ts_out = t_seq[ts_ind]
ts_out_maxS <- ts_out

out <- simple_water_model_sim_time_breaks(param, ts_out_maxS, deltat)

tcrit_out_maxS <- t_seq[which.min(out$W)]

out_maxS_1000 <- data.frame(carbon_pool_type = as.factor(rep(c("Storage", "Biomass"), each = simlength)),
                            carbon_pool_val = c(out$S, out$M),
                            water_pool_val = rep(out$W, times = 2),
                            water_pool_type = as.factor(rep("Soil Water Available", times = simlength * 2)),
                            initial_water = rep(param$W0, simlength * 2),
                            ts = (rep(ts_out_maxS, times = simlength * 2)),
                            tcrit = (rep(tcrit_out_maxS, times = simlength * 2)),
                            goal = rep("Storage", simlength * 2),
                            time = rep(t_seq, times = 2))

kf = 1
ms_sum <- colSums(rbind (kf * M_out_2, (1 - kf) * S_out_2))
ts_ind <- which.max(ms_sum)
ts_out = t_seq[ts_ind]
ts_out_maxM <- ts_out

out <- simple_water_model_sim_time_breaks(param_2, ts_out_maxM, deltat)

tcrit_out_maxM <- t_seq[which.min(out$W)]

out_maxM_3000 <- data.frame(carbon_pool_type = as.factor(rep(c("Storage", "Biomass"), each = simlength)),
                            carbon_pool_val = c(out$S, out$M),
                            water_pool_val = rep(out$W, times = 2),
                            water_pool_type = as.factor(rep("Soil Water Available", times = simlength * 2)),
                            initial_water = rep(param_2$W0, simlength * 2),
                            ts = (rep(ts_out_maxM, times = simlength * 2)),
                            tcrit = (rep(tcrit_out_maxM, times = simlength * 2)),
                            goal = rep("Biomass", simlength * 2),
                            time = rep(t_seq, times = 2))

kf = 0
ms_sum <- colSums(rbind (kf * M_out_2, (1 - kf) * S_out_2))
ts_ind <- which.max(ms_sum)
ts_out = t_seq[ts_ind]
ts_out_maxS <- ts_out

out <- simple_water_model_sim_time_breaks(param_2, ts_out_maxS, deltat)

tcrit_out_maxS <- t_seq[which.min(out$W)]

out_maxS_3000 <- data.frame(carbon_pool_type = as.factor(rep(c("Storage", "Biomass"), each = simlength)),
                            carbon_pool_val = c(out$S, out$M),
                            water_pool_val = rep(out$W, times = 2),
                            water_pool_type = as.factor(rep("Soil Water Available", times = simlength * 2)),
                            initial_water = rep(param_2$W0, simlength * 2),
                            ts = (rep(ts_out_maxS, times = simlength * 2)),
                            tcrit = (rep(tcrit_out_maxS, times = simlength * 2)),
                            goal = rep("Storage", simlength * 2),
                            time = rep(t_seq, times = 2))

allres <- rbind(out_maxM_1000, out_maxS_1000, out_maxM_3000, out_maxS_3000)



coeff <- 0.4





mean_ts <- data.frame(
  ts = c(out_maxM_1000$ts[1], out_maxS_1000$ts[1], out_maxM_3000$ts[1], out_maxS_3000$ts[1]), 
  goal=c("Biomass", "Storage", "Biomass", "Storage"), 
  initial_water=c(param$W0, param$W0, param_2$W0, param_2$W0)
  )

mean_tcrit <- data.frame(
  tcrit = c(out_maxM_1000$tcrit[1], out_maxS_1000$tcrit[1], out_maxM_3000$tcrit[1], out_maxS_3000$tcrit[1]), 
  goal=c("Biomass", "Storage", "Biomass", "Storage"), 
  initial_water=c(param$W0, param$W0, param_2$W0, param_2$W0)
)





labels <- data.frame(
  lab = c("A 1000kgH20", "B 1000kgH20", "C 3000kgH20", "D 3000kgH20"),
  goal=c("Biomass", "Storage", "Biomass", "Storage"), 
  initial_water=c(param$W0, param$W0, param_2$W0, param_2$W0)
)

p <- ggplot(allres, aes(x=time, y=carbon_pool_val)) +
  geom_line(mapping=aes(colour=carbon_pool_type), show.legend=TRUE) + 
  geom_line(mapping=aes(y=water_pool_val/coeff, 
                        fill = water_pool_type), size=1.5, alpha=0.5, color = "blue") + 
  geom_hline(yintercept=0, color="#000000", size = 0.25,linetype="dashed") + 
  geom_vline(aes(xintercept=ts), mean_ts, color="red", size = 0.25,linetype="dashed") +
  geom_text(aes(x=15, label=lab, y=7750), 
            labels, text=element_text(size=8), vjust = 1,) +
  # geom_vline(aes(xintercept=tcrit), mean_tcrit, color="blue", size = 0.25,linetype="dashed") +
  # geom_text(aes(x=tcrit - 15, label=paste("tcrit = ", tcrit), y=6000), 
  #           mean_tcrit, colour="blue", angle=90, text=element_text(size=8), vjust = 1,) +
  scale_y_continuous(
    # limits = c(-625,4525),
    # Features of the first axis
    "Carbon Pool Size (gC)", 
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name = expression(paste("Available Water (kg", H[2], "O)")))
  ) +
  scale_color_hue("Carbon Pool", guide=guide_legend(order=1)) +
  scale_fill_manual(name  ="Water Pool", values = rep(1, 3),
                    guide=guide_legend(
                      override.aes = list(colour=c("blue")) , order=2),
                    labels=c("Available\nSoil Water")) +
  labs(x = "Time (day)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p 


p <- p + facet_grid(cols = vars(goal), rows=vars(initial_water), space= "fixed") +
#, label_both, am = label_value)
 #              labeller = label_bquote(cols = paste("Maximising ", .(goal)),
  #                                     rows = paste(W[0], " = ", .(initial_water), " kg", H[2], "O"))) +
  theme(strip.background = element_blank(), strip.text.y = element_blank())

p


t_df <- data.frame(
  t = c(mean_ts$ts, mean_tcrit$ts),
  t_type = rep(c("t[s]", "t[crit]"), each=4),
  goal = c(mean_ts$goal, mean_tcrit$goal),
  initial_water = c(mean_ts$initial_water, mean_tcrit$initial_water)
)

t_df <- mean_ts
t_df$tcrit <- mean_tcrit$tcrit
rownames(t_df) <- paste(t_df$goal, t_df$initial_water)


# Set theme to allow for plotmath expressions
tt <- ttheme_default(padding = unit(c(8, 8), "mm"),
                     colhead=list(fg_params=list(parse=TRUE)),
                     rowhead=list(fg_params=list(parse=TRUE)))
tbl <- tableGrob(subset(t_df, select=c(ts, tcrit)), 
                 rows = c("Phi[M]",
                          "Phi[S]",
                          "Phi[M]",
                          "Phi[S]"), 
                 cols = c("t[s]", "t[crit]"),
                 theme=tt)


grid.arrange(p, tbl,widths = 2:1)
             # as.table=TRUE)
