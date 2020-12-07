rm(list=ls())

library(ggplot2)
library(wesanderson)
library(ggpubr)
library(lemon)
library(cowplot)
library(plyr)

source("code/parameters.R")
source("code/models/water_model.R")

param <- param_list_drake_W0_1000()

deltat = 0.1
t_seq <- seq(from=0, to=param$Tend, by=deltat)
simlength <- length(t_seq)

ts_1 = 10;
ts_2 = 20;
ts_3 = 30;


out_1 <- simple_water_model_sim_time_breaks(param, ts_1, deltat)
out_2 <- simple_water_model_sim_time_breaks(param, ts_2, deltat)
out_3 <- simple_water_model_sim_time_breaks(param, ts_3, deltat)




out_df_1 <- data.frame(carbon_pool_type = as.factor(rep(c("Storage", "Biomass"), each = simlength)),
                            carbon_pool_val = c(out_1$S, out_1$M),
                            carbon_pool_perc = c(out_1$S, out_1$M)/rep(out_1$S+out_1$M, times=2),
                            water_pool_val = rep(out_1$W, times = 2),
                            water_pool_type = as.factor(rep("Soil Water Available", times = simlength * 2)),
                            ts = (rep(ts_1, times = simlength * 2)),
                            time = rep(t_seq, times = 2))
out_df_2 <- data.frame(carbon_pool_type = as.factor(rep(c("Storage", "Biomass"), each = simlength)),
                            carbon_pool_val = c(out_2$S, out_2$M),
                            carbon_pool_perc = c(out_2$S, out_2$M)/rep(out_2$S+out_2$M, times=2),
                            water_pool_val = rep(out_2$W, times = 2),
                            water_pool_type = as.factor(rep("Soil Water Available", times = simlength * 2)),
                            ts = (rep(ts_2, times = simlength * 2)),
                            time = rep(t_seq, times = 2))
out_df_3 <- data.frame(carbon_pool_type = as.factor(rep(c("Storage", "Biomass"), each = simlength)),
                            carbon_pool_val = c(out_3$S, out_3$M),
                            carbon_pool_perc = c(out_3$S, out_3$M)/rep(out_3$S+out_3$M, times=2),
                            water_pool_val = rep(out_3$W, times = 2),
                            water_pool_type = as.factor(rep("Soil Water Available", times = simlength * 2)),
                            ts = (rep(ts_3, times = simlength * 2)),
                            time = rep(t_seq, times = 2))


t_val = t_seq[which(out_1$S < 0)]
out_df_1$neg[out_df_1$time %in% t_val] = 
  out_df_1$carbon_pool_val[out_df_1$time %in% t_val]
t_val = t_seq[which(out_2$S < 0)]
out_df_2$neg[out_df_2$time %in% t_val] = 
  out_df_2$carbon_pool_val[out_df_2$time %in% t_val]
t_val = t_seq[which(out_3$S < 0)]
out_df_3$neg[out_df_3$time %in% t_val] = 
  out_df_3$carbon_pool_val[out_df_3$time %in% t_val]
out_df_3$carbon_pool_val[out_df_3$time %in% t_val] = 
  NA




allres <- rbind.fill(out_df_1,out_df_2,out_df_3)
#allres$ts <- as.factor(allres$ts)


coeff <- 0.4
colors <- wes_palette("Darjeeling1", n = 3)






p1 <- ggplot(out_df_1, aes(x=time, y=carbon_pool_val)) +
  geom_line(mapping=aes(colour=carbon_pool_type), show.legend=TRUE) + 
  geom_line(mapping=aes(y=water_pool_val/coeff, 
                        fill = water_pool_type), size=1.5, alpha=0.5, color = "blue") + 
  
  geom_hline(yintercept=0, color="#000000", size = 0.25,linetype="dashed") + 
  scale_y_continuous(
    limits = c(-625,4525),
    # Features of the first axis
    name = "Pool Size (gC)",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff,labels = NULL)
  ) +
  scale_color_hue("Carbon Pool", guide=guide_legend(order=1)) +
  scale_fill_manual(name  ="Water Pool", values = rep(1, 3),
                    guide=guide_legend(
                      override.aes = list(colour=c("blue")) , order=2),
                    labels=c("Available Soil Water")) +
  
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


p1


p2 <- ggplot(out_df_2, aes(x=time, y=carbon_pool_val)) +
  geom_line(mapping=aes(colour=carbon_pool_type), show.legend=TRUE) + 
  geom_line(mapping=aes(y=water_pool_val/coeff, 
                        fill = water_pool_type), size=1.5, alpha=0.5, color = "blue") + 
  geom_hline(yintercept=0, color="#000000", size = 0.25,linetype="dashed") + 
  scale_y_continuous(
    limits = c(-625,4525),
    # Features of the first axis
    name = "", labels = NULL,
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff,labels = NULL)
  ) +
  scale_color_hue("Carbon Pool", guide=guide_legend(order=1)) +
  scale_fill_manual(name  ="Water Pool", values = rep(1, 3),
                    guide=guide_legend(
                      override.aes = list(colour=c("blue")) , order=2),
                    labels=c("Available Soil Water")) +
  
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


p2

p3 <- ggplot(out_df_3, aes(x=time, y=carbon_pool_val)) +
  geom_line(mapping=aes(colour=carbon_pool_type), show.legend=TRUE) + 
  geom_line(mapping=aes(y=water_pool_val/coeff, 
                        fill = water_pool_type), size=1.5, alpha=0.5, color = "blue") + 
  geom_hline(yintercept=0, color="#000000", size = 0.25,linetype="dashed") + 
  scale_y_continuous(
    limits = c(-625,4525),
    # Features of the first axis
    name = "", labels = NULL,
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name = expression(paste("Available Water (kg", H[2], "O)")))
  ) +
  scale_color_hue("Carbon Pool", guide=guide_legend(order=1)) +
  scale_fill_manual(name  ="Water Pool", values = rep(1, 3),
                    guide=guide_legend(
                      override.aes = list(colour=c("blue")) , order=2),
                    labels=c("Available Soil Water")) +
  
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


p3

p3 + geom_line(aes(y=neg, colour = neg_carbon_pool),na.rm = TRUE, size=4)

ggarrange(p1, p2, p3, ncol=3, nrow=1, 
          labels = c("A", "B", "C"), common.legend = TRUE, legend = "bottom")

plot_grid(p1, p2, p3, labels=c("A", "B", "C"), ncol = 3, nrow = 1)




#FIND NEGATIVE VALUES IN OUT 3

allres$neg[allres$carbon_pool_type == "Storage" & allres$carbon_pool_val < 0] = allres$carbon_pool_val[allres$carbon_pool_type == "Storage" & allres$carbon_pool_val < 0]




wes_palette_mine <- function(n){
  return(wes_palette(n, "Darjeeling1"))
}

#### Plotting percentage of total biomass


p <- ggplot(allres, aes(x=time, y=carbon_pool_perc)) +
  geom_line(mapping=aes(colour=carbon_pool_type), show.legend=TRUE) + 
  # geom_line(aes(y=neg, group = carbon_pool_type),na.rm = TRUE, linetype="dotted") +
  geom_hline(yintercept=0, color="#000000", size = 0.25,linetype="dashed") + 
  # geom_vline(xintercept=ts, color="#000000", size = 0.25,linetype="dashed") + 
  scale_y_continuous(
    # limits = c(-625,4525),
    # Features of the first axis
    "Proportion of carbon pool per total plant size", 
    # Add a second axis and specify its features
    # sec.axis = sec_axis(~.*coeff, name = expression(paste("Available Water (kg", H[2], "O)")))
  ) +
  scale_color_hue("Carbon Pool", guide=guide_legend(order=1)) +
  # scale_fill_manual(name  ="Water Pool", values = rep(1, 3),
                    # guide=guide_legend(
                      # override.aes = list(colour=c("blue")) , order=2),
                    # labels=c("Available Soil Water")) +
  labs(x = "Time (day)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p 


p + facet_grid(cols = vars(ts), space= "fixed", 
               labeller = label_bquote(cols = paste(t[s], " = ", .(ts)))) +
  theme(strip.background = element_blank())






#### Plotting the original
 
p <- ggplot(allres, aes(x=time, y=carbon_pool_val)) +
  geom_line(mapping=aes(colour=carbon_pool_type), show.legend=TRUE) + 
  geom_line(mapping=aes(y=water_pool_val/coeff, 
                        fill = water_pool_type), size=1.5, alpha=0.5, color = "blue") + 
  geom_line(aes(y=neg, group = carbon_pool_type),na.rm = TRUE, linetype="dotted") +
  geom_hline(yintercept=0, color="#000000", size = 0.25,linetype="dashed") + 
  # geom_vline(xintercept=ts, color="#000000", size = 0.25,linetype="dashed") + 
  scale_y_continuous(
    limits = c(-625,4525),
    # Features of the first axis
    "Carbon Pool Size (gC)", 
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*coeff, name = expression(paste("Available Water (kg", H[2], "O)")))
  ) +
  scale_color_hue("Carbon Pool", guide=guide_legend(order=1)) +
  scale_fill_manual(name  ="Water Pool", values = rep(1, 3),
                    guide=guide_legend(
                      override.aes = list(colour=c("blue")) , order=2),
                    labels=c("Available Soil Water")) +
  labs(x = "Time (day)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p 


p + facet_grid(cols = vars(ts), space= "fixed", 
               labeller = label_bquote(cols = paste(t[s], " = ", .(ts)))) +
  theme(strip.background = element_blank())
