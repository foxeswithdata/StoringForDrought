source("./code/models/model.R")
source("./code/models/water_model.R")
source("./code/parameters.R")

param <- param_list_water_short();

dt = 0.0001;
ts = 20;

u <- c(rep(param$ks, times=ts/dt), rep(0, times = param$Tend/dt - ts/dt))
u2 <- c(rep(param$ks, times=ts), rep(0, times=param$Tend-ts))

out <- simple_water_model_sim_time_breaks(param$M0, param$S0, param$W0, param$Tend, param$kp, param$kr, param$ks, param$kw, ts, deltat=dt)
out2 <- simple_water_model_sim_time_breaks(param$M0, param$S0, param$W0, param$Tend, param$kp, param$kr, param$ks, param$kw, ts, deltat=1)


plot(out$t, out$S, type="l")
plot(out2$t, out2$S, type="l")


