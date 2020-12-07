# set up ggplot
library(ggplot2)
library(ggpubr)
theme_set(
  theme_pubr() +
    theme(legend.position = "right") +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"),)
)

#import all necessary models and parameters
source('code/models/model_water_kp_rh/model_fb.R')
source('code/models/model_water_kp_rh/model_sim.R')
source('code/models/water_model.R')
source('code/parameters.R')
  
# import all the parameters we're going to work with
param <- param_list_water_short2()

#### Start with what the kp(W) function looks like
W <- seq(0, 100, by=0.01)
kp <- (param$kp * W)/(param$k + W)

smallk <- 1
bigk <- 10
kp_smallk <- (param$kp * W)/(smallk + W)
kp_bigk <- (param$kp * W)/(bigk + W)
 
kp_func <- data.frame(W = c(W,W,W), kp = c(kp,kp_smallk,kp_bigk), 
                      k_val = as.factor(c(rep(param$k, times = length(W)),
                                          rep(smallk, times = length(W)),
                                          rep(bigk, times = length(W)))))

p <- ggplot(kp_func, aes(x = W, y = kp, color=k_val)) +
  geom_line() +
  geom_hline(yintercept=param$kp, linetype="dashed")
  
p <- p + labs(title = "Photosynthetic Parameter Function") +
              xlab(expression(paste("Available Soil Water [kg", H[2], "0]"))) +
              ylab(expression(paste("Photosynthetic Parameter ", k[p], " []")))
p

# Let's see how the solution looks
out <- simulation_kp_rh(param, 40, deltat) 

# How does my simple simulation looks like
out1 <- simulation_kp_rh_simple(param, 10, deltat)
out2 <- simulation_kp_rh_simple(param, 20, deltat)
out3 <- simulation_kp_rh_simple(param, 30, deltat)
out4 <- simulation_kp_rh_simple(param, 40, deltat)
out5 <- simulation_kp_rh_simple(param, 50, deltat)
out6 <- simulation_kp_rh_simple(param, 60, deltat)

plot(out$t, out$S, type='l')
plot(out1$t, out1$S, type="l")
# plot(out2$t, out2$S, type="l", col="blue")
lines(out2$t, out2$S, col="blue")
lines(out3$t, out3$S, col="red")
lines(out4$t, out4$S, col="magenta")
lines(out5$t, out5$S, col="green")
lines(out6$t, out6$S, col="grey")


# look at all the results for different outputs (W, M, S etc)

plot(out3$t, out3$M, type="l", col="blue", 
     ylim = c(min(c(out3$M, out3$S, out3$W)), max(c(out3$M, out3$S, out3$W))))
lines(out3$t, out3$S, col="red")
lines(out3$t, out3$W, col="green")

# Let's compare a single result with ts = 30

out1 <- water_model_log(param, 30, deltat)
out2 <- simple_water_model_sim_time_breaks(param, 30, deltat)
out3 <- simulation_kp_rh_simple(param, 30, deltat)

plot(out1$t, out1$M, type="l", col="blue", 
     ylim = c(min(c(out1$M, out2$M, out1$M)), max(c(out1$M, out2$M, out1$M))))
lines(out2$t, out2$M, col="red")
lines(out3$t, out3$M, col="green")

plot(out1$t, out1$S, type="l", col="blue", 
     ylim = c(min(c(out1$S, out2$S, out1$S)), max(c(out1$S, out2$S, out1$S))))
lines(out2$t, out2$S, col="red")
lines(out3$t, out3$S, col="green")

plot(out1$t, out1$W, type="l", col="blue", 
     ylim = c(min(c(out1$W, out2$W, out1$W)), max(c(out1$W, out2$W, out1$W))))
lines(out2$t, out2$W, col="red")
lines(out3$t, out3$W, col="green")



 # Look at how the W0 value changes

W0_list =  seq(from=10, to=10000, by=10)
W0_list_testing = seq(from=100, to=500, by=10)

ts_maxs <- 0
mT_maxs <- 0
sT_maxs <- 0
wT_maxs <- 0
A_maxs <- 0

ts_maxs_log <- 0
mT_maxs_log <- 0
sT_maxs_log <- 0
wT_maxs_log <- 0
A_maxs_log <- 0

ts_maxs_rh_1 <- 0
mT_maxs_rh_1 <- 0
sT_maxs_rh_1 <- 0
wT_maxs_rh_1 <- 0
A_maxs_rh_1 <- 0

ts_maxs_rh_2 <- 0
mT_maxs_rh_2 <- 0
sT_maxs_rh_2 <- 0
wT_maxs_rh_2 <- 0
A_maxs_rh_2 <- 0

ts_maxs_rh_3 <- 0
mT_maxs_rh_3 <- 0
sT_maxs_rh_3 <- 0
wT_maxs_rh_3 <- 0
A_maxs_rh_3 <- 0

ts_maxs_rh_4 <- 0
mT_maxs_rh_4 <- 0
sT_maxs_rh_4 <- 0
wT_maxs_rh_4 <- 0
A_maxs_rh_4 <- 0

ts_maxs_rh_5 <- 0
mT_maxs_rh_5 <- 0
sT_maxs_rh_5 <- 0
wT_maxs_rh_5 <- 0
A_maxs_rh_5 <- 0

# Iterate through the list of examined Water Inputs W0_list
# Save j as the index of final outputs
j = 1

for(w in W0_list){
  
  # Set simulation time tolerance. The smaller this value the
  # more precise the simulation
  deltat = 0.1
  
  # Initialise a vector of possible ts values and the length of this
  # vector
  t_seq <- seq(from=0, to=param$Tend, by=deltat)
  simlength <- length(t_seq)
  
  # Initialise output vectors, these will contain the
  # Final outputs of simulations for each ts value
  S_out <- vector(length=simlength)
  S_out_log <- vector(length=simlength)
  S_out_rh_1 <- vector(length=simlength)
  S_out_rh_2 <- vector(length=simlength)
  S_out_rh_3 <- vector(length=simlength)
  S_out_rh_4 <- vector(length=simlength)
  S_out_rh_5 <- vector(length=simlength)
  
  # Iterate through the possible ts values
  k = 1
  for (i in t_seq){
    # Find the output values of simulation for the specific ts value
    param_2 <- param
    param_2$W0 <- w
    param_2$x0 <- c(param_2$M0, param_2$S0, param_2$W0)
    
    out <- simple_water_model_sim_time_breaks(param_2, i, deltat)
    out_log <- water_model_log(param_2, i, deltat)
    out_rh <- simulation_kp_rh(param_2, i, deltat)
    
    # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
    S_out[k] <- out$S[simlength]
    S_out_log[k] <- out_log$S[simlength]
    S_out_rh_1[k] <- out_rh$S[simlength]
    
    param_2$krh <- 0.5
    out_rh_2 <- simulation_kp_rh(param_2, i, deltat)
    param_2$krh <- 2
    out_rh_3 <- simulation_kp_rh(param_2, i, deltat)
    param_2$krh <- 7.5
    out_rh_4 <- simulation_kp_rh(param_2, i, deltat)
    param_2$krh <- 10
    out_rh_5 <- simulation_kp_rh(param_2, i, deltat)
    
    S_out_rh_2[k] <- out_rh_2$S[simlength]
    S_out_rh_3[k] <- out_rh_3$S[simlength]
    S_out_rh_4[k] <- out_rh_4$S[simlength]
    S_out_rh_5[k] <- out_rh_5$S[simlength]
    
    if(sum(out$S<0) > 0){
      S_out[k] <- NA
    }
    if(sum(out_log$S<0) > 0){
      S_out_log[k] <- NA
    }
    if(sum(out_rh$S < 0) > 0){
      S_out_rh_1[k] <- NA
    }
    if(sum(out_rh_2$S < 0) > 0){
      S_out_rh_2[k] <- NA
    }
    if(sum(out_rh_3$S < 0) > 0){
      S_out_rh_3[k] <- NA
    }
    if(sum(out_rh_4$S < 0) > 0){
      S_out_rh_4[k] <- NA
    }
    if(sum(out_rh_5$S < 0) > 0){
      S_out_rh_5[k] <- NA
    }
    k = k+1
  }
  
  # Find which ts maximises final storage
  # Because maximising storage will eliminate the
  # results that are unrealistic (unless only those exist)
  # there is no need to limit the search
  # switch model
  ts_ind <- which.max(S_out)
  if(length(ts_ind) == 0){
    ts_ind = 1
  }
  ts_maxs[j] <- t_seq[ts_ind]
  sT_maxs[j] <- S_out[ts_ind]
  
  # log model
  ts_ind <- which.max(S_out_log)
  if(length(ts_ind) == 0){
    ts_ind = 1
  }
  ts_maxs_log[j] <- t_seq[ts_ind]
  sT_maxs_log[j] <- S_out_log[ts_ind]
  
  # rh model
  
  ts_ind <- which.max(S_out_rh_1)
  if(length(ts_ind) == 0){
    ts_ind = 1
  }
  ts_maxs_rh_1[j] <- t_seq[ts_ind]
  sT_maxs_rh_1[j] <- S_out_rh_1[ts_ind]
  
  # rh model different krhs
  
  ts_ind <- which.max(S_out_rh_2)
  if(length(ts_ind) == 0){
    ts_ind = 1
  }
  ts_maxs_rh_2[j] <- t_seq[ts_ind]
  sT_maxs_rh_2[j] <- S_out_rh_2[ts_ind]
  
  ts_ind <- which.max(S_out_rh_3)
  if(length(ts_ind) == 0){
    ts_ind = 1
  }
  ts_maxs_rh_3[j] <- t_seq[ts_ind]
  sT_maxs_rh_3[j] <- S_out_rh_3[ts_ind]
  
  ts_ind <- which.max(S_out_rh_4)
  if(length(ts_ind) == 0){
    ts_ind = 1
  }
  ts_maxs_rh_4[j] <- t_seq[ts_ind]
  sT_maxs_rh_4[j] <- S_out_rh_4[ts_ind]
  
  ts_ind <- which.max(S_out_rh_5)
  if(length(ts_ind) == 0){
    ts_ind = 1
  }
  ts_maxs_rh_5[j] <- t_seq[ts_ind]
  sT_maxs_rh_5[j] <- S_out_rh_5[ts_ind]
  
  #update W0 index
  j = j + 1;
}

maxS_strategy_w0_var_rh <- data.frame(w0 = rep(W0_list, times=5),
                                      ts = c(ts_maxs_rh_2, ts_maxs_rh_3, ts_maxs_rh_1, ts_maxs_rh_4, ts_maxs_rh_5),
                                      sT = c(sT_maxs_rh_2, sT_maxs_rh_3, sT_maxs_rh_1, sT_maxs_rh_4, sT_maxs_rh_5),
                                      rh = as.factor(rep(c(0.5, 2, 5, 7.5, 10), each = length(ts_maxs_rh_1)))
                                      )


p <- ggplot(maxS_strategy_w0_var_rh, aes(w0, ts, col=rh)) +
  geom_line() +
  scale_y_continuous(name = "time of allocation switch ts [day]", breaks = 0:15, limits = c(0,15)) + 
  scale_x_continuous(name = expression(paste("Initial Water Content [kg", H[2], "O]")), limits = c(200,2000)) + 
  labs(title= expression(paste("Allocation Response To Initial Water Content")),
       subtitle= expression(paste("Varying RH parameter ", k[rh])),
       colour= expression(k[rh]))
p


maxS_strategy_w0 <- data.frame(w0 = W0_list_testing, 
                               mT = mT_maxs, 
                               ts = ts_maxs, 
                               sT = sT_maxs, 
                               a = A_maxs, 
                               wT = wT_maxs
                               )

maxS_strategy_w0_log <- data.frame(w0 = W0_list_testing, 
                               mT = mT_maxs_log, 
                               ts = ts_maxs_log, 
                               sT = sT_maxs_log, 
                               a = A_maxs_log, 
                               wT = wT_maxs_log
)

maxS_strategy_w0_rh <- data.frame(w0 = W0_list_testing, 
                               mT = mT_maxs_rh, 
                               ts = ts_maxs_rh, 
                               sT = sT_maxs_rh, 
                               a = A_maxs_rh, 
                               wT = wT_maxs_rh
)
p <- ggplot(maxS_strategy_w0, aes(x = w0, y = ts)) +
  geom_line() 
p

p <- ggplot(maxS_strategy_w0_log, aes(x = w0, y = ts)) +
  geom_line() 
p

p <- ggplot(maxS_strategy_w0_rh, aes(x = w0, y = ts)) +
  geom_line() 
p

p <- ggplot(maxS_strategy_w0, aes(x = w0, y = ts), xlim=c(0,15), ) +
  geom_line() +
  geom_line(data = maxS_strategy_w0_log, aes(x = w0, y = ts), col="blue") + 
  geom_line(data = maxS_strategy_w0_rh, aes(x = w0, y = ts), col="red") + 
  scale_y_continuous(name = "time of allocation switch ts [day]", breaks = 0:15, limits = c(0,15)) + 
  scale_x_continuous(name = expression(paste("Initial Water Content [kg", H[2], "O]"))) + 
  labs(title="Allocation Response To Initial Water Content")
p

#########################################
## TEST TS RESPONSE TO ONE VALUE OF W0 ##
#########################################

# using default krh value


# Set simulation time tolerance. The smaller this value the
# more precise the simulation
deltat = 0.01

# Initialise a vector of possible ts values and the length of this
# vector
t_seq <- seq(from=0, to=param$Tend, by=deltat)
simlength <- length(t_seq)

# Initialise output vectors, these will contain the
# Final outputs of simulations for each ts value
S_out <- vector(length=simlength)
S_out_log <- vector(length=simlength)
S_out_rh <- vector(length=simlength)

M_out <- vector(length=simlength)
M_out_log <- vector(length=simlength)
M_out_rh <- vector(length=simlength)

W_out <- vector(length=simlength)
W_out_log <- vector(length=simlength)
W_out_rh <- vector(length=simlength)

A_out <- vector(length=simlength)
A_out_log <- vector(length=simlength)
A_out_rh <- vector(length=simlength)

# Iterate through the possible ts values
k = 1
for (i in t_seq){
  # Find the output values of simulation for the specific ts value
  param_2 <- param
  param_2$W0 <- 400
  
  out <- simple_water_model_sim_time_breaks(param_2, i, deltat)
  out_log <- water_model_log(param_2, i, deltat)
  out_rh <- simulation_kp_rh(param_2, i, deltat)
  
  # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
  S_out[k] <- out$S[simlength]
  S_out_log[k] <- out_log$S[simlength]
  S_out_rh[k] <- out_rh$S[simlength]
  
  M_out[k] <- out$M[simlength]
  M_out_log[k] <- out_log$M[simlength]
  M_out_rh[k] <- out_rh$M[simlength]
  
  W_out[k] <- out$W[simlength]
  W_out_log[k] <- out_log$W[simlength]
  W_out_rh[k] <- out_rh$W[simlength]
  
  A_out[k] <- sum(out$A)
  A_out_log[k] <- sum(out_log$A)
  A_out_rh[k] <- out_rh$A[simlength]
  
  if(sum(out$S<0) > 0){
    S_out[k] <- NA
    M_out[k] <- NA
    W_out[k] <- NA
    A_out[k] <- NA
  }
  if(sum(out_log$S<0) > 0){
    S_out_log[k] <- NA
    M_out_log[k] <- NA
    W_out_log[k] <- NA
    A_out_log[k] <- NA
  }
  if(sum(out_rh$S < 0) > 0){
    S_out_rh[k] <- NA
    M_out_rh[k] <- NA
    W_out_rh[k] <- NA
    A_out_rh[k] <- NA
  }
  
  k = k+1
}

results_w0_400 <- data.frame(ts = rep(t_seq, times = 3),
                            STout = c(S_out, S_out_log, S_out_rh),
                            MTout = c(M_out, M_out_log, M_out_rh),
                            WTout = c(W_out, W_out_log, W_out_rh),
                            ATout = c(A_out, A_out_log, A_out_rh),
                            type = as.factor(rep(c("switch", "log", "rh"), each=simlength))
                            )

p <- ggplot(results_w0_400, aes(ts, STout, col=type)) +
  geom_line() +
  scale_x_continuous(name = "time of allocation switch ts [day]") + 
  scale_y_continuous(name = expression(paste("Final Storage Value [kgC]"))) + 
  geom_point(size = 2, 
             aes(x = t_seq[which.max(S_out)], 
                 y = max(S_out), 
                 label = paste("MaxST [",t_seq[which.max(S_out)], ",", max(S_out), "]"))
             ) +
  labs(title= expression(paste("Final Storage Output vs Time Switch ", t[s])),
       subtitle= expression(paste("Initial Water Content 400 [kg", H[2], "O]")),
       colour= expression(paste(k[p], " function type")))
p

t_seq[which.max(S_out)]

p <- ggplot(results_w0_400, aes(ts, MTout, col=type)) +
  geom_line() +
  scale_x_continuous(name = "time of allocation switch ts [day]") + 
  scale_y_continuous(name = expression(paste("Final Biomass Value [kgC]"))) + 
  labs(title= expression(paste("Final Biomass Output vs Time Switch ", t[s])),
       subtitle= expression(paste("Initial Water Content 400 [kg", H[2], "O]")),
       colour= expression(paste(k[p], " function type")))
p


