rm(list=ls())

alpha_plus <- function(ks, kp, kr){
  a_p = (-ks + sqrt(ks + 4 * ks * (kp - kr)))/2
  return (a_p)
}

alpha_minus <- function(ks, kp, kr){
  a_m = (-ks - sqrt(ks + 4 * ks * (kp - kr)))/2
  return (a_m)
}

beta_plus <- function(S0, M0, ks, kp, kr){
  b_p = (S0*ks - M0 * alpha_minus(ks, kp, kr))/(alpha_plus(ks, kp, kr)- alpha_minus(ks, kp, kr))
  return(b_p)
}

beta_minus <- function(S0, M0, ks, kp, kr){
  b_m = M0 - (S0*ks - M0 * alpha_minus(ks, kp, kr))/(alpha_plus(ks, kp, kr)- alpha_minus(ks, kp, kr))
  return(b_m)
}

beta_plus_2 <- function(S0, M0, ks, kp, kr){
  b_p = M0 - (S0 * ks - M0 * alpha_plus(ks, kp, kr))/(alpha_minus(ks, kp, kr)- alpha_plus(ks, kp, kr))
  return(b_p)
}

#not working
ts_fun_orig <- function(kr, T_end, t_star, kp, kw, nu){
  ts = (1 + kr * T_end - t_star * (kp * (kw * nu - 1)))/(kw * nu * kp - kp + kr)
  return(ts)
}

ts_fun_fixed <- function(kr, T_end, t_star, kp, kw, nu){
  ts = (1 + kr * T_end + (nu * kw - 1) * kp * t_star)/(kr + (nu * kw - 1) * kp)
  return(ts)
}

#not working
tstar_fun_orig <- function(kr, T_end, ts, kp, kw, nu){
  tstar = (1 + kr * T_end + ts * (kp * (1 - kw * nu) - kr))/ (kp * (kw * nu - 1))
  return (tstar)
}

tstar_fun_fixed <- function(kr, T_end, ts, kp, kw, nu){
  tstar = ((kr + (nu * kw - 1) * kp) * ts - kr * T_end - 1)/((nu * kw - 1)* kp)
  return(tstar)
}

tstar_fun_notebook <- function(kr, T_end, ts, kp, kw, nu){
  tstar = (ts * (kp * (kw * nu - 1) + kr) - 1 - kr * T_end)/((kw * nu -1) * kp)
  return(tstar)
}

ST_func <- function(ts, nu, S0, M0, ks, kp, kr, kw, T_end){
  ap <- alpha_plus(ks, kp, kr)
  am <- alpha_minus(ks, kp, kr)
  bp <- beta_plus(S0, M0, ks, kp, kr)
  bm <- beta_minus(S0, M0, ks, kp, kr)
  SI_t_s = 1/ks * (bp * ap * exp(ap * ts) + bm * am * exp(am * ts))
  MI_t_s = bp * exp(ap * ts) + (M0 - bm) * exp(am * ts)
  tstar <- tstar_fun_fixed(kr, T_end, ts, kp, kw, nu)
  SII_t_star = SI_t_s + (kp - kr) * MI_t_s * (tstar - ts)
  SIII_T_end = SII_t_star - kr * MI_t_s * (T_end + tstar)
  return(SIII_T_end)
}

Wt_star_func <- function(ts, nu, W0, S0, M0, ks, kp, kr, kw, T_end){
  ap <- alpha_plus(ks, kp, kr)
  am <- alpha_minus(ks, kp, kr)
  bp <- beta_plus(S0, M0, ks, kp, kr)
  bm <- beta_minus(S0, M0, ks, kp, kr)
  tstar <- tstar_fun_fixed(kr, T_end, ts, kp, kw, nu)
  
  WII_star = W0 - kp * kw * (exp(ap * ts) * (bp * (1 - ap * (tstar - ts)))/(ap * (tstar - ts)) 
                           + (bm * (1 - am* (tstar - ts)))/(am * (tstar - ts))) * (tstar - ts);
  
  return(WII_star);
}


source("code/parameters.R")




nu_ar <- seq(-10, 10, 0.05);
ts <- seq(0, param_1$T_end/2, 0.5)

Wtstar <- matrix(nrow = length(nu_ar), ncol=length(ts))
ST_res <- matrix(nrow = length(nu_ar), ncol=length(ts))

#currently only working for half of the results??

# for(i in nu_ar){
#   for(j in ts){
#     nu_ind <- which(nu_ar == i)
#     
#     ts_ind <- which(ts == j)
#     
#     print(ts_ind)
#     
#     Wtstar[i,j] = Wt_star_func(j,i, param_1$W0, param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr, param_1$kw, param_1$T_end)
#     ST_res[i,j] = ST_func(j, i, param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr, param_1$kw, param_1$T_end)
#   }
# }

nu = -1000

Wts <- vector(length=30);
STts <- vector(length=30);

# for (i in 1:30){
#   Wts[i] = Wt_star_func(i,nu, param_1$W0, param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr, param_1$kw, param_1$T_end)
#   STts[i] = ST_func(i,nu, param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr, param_1$kw, param_1$T_end)
# }


plot(1:30, Wts)
plot(1:30, STts)

print(STts>0)

rownames(Wtstar) <- nu_ar
colnames(Wtstar) <- ts

library(reshape2)
melted_Wtstar <- melt(Wtstar)
head(melted_Wtstar)

melted_ST <- melt(ST_res)
head(melted_ST)



library(ggplot2)
ggplot(data = melted_Wtstar, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()

ggplot(data = melted_ST, aes(x=Var2, y=Var1, fill=value)) + 
  geom_tile()



source("code/models/water_model.R")

#test for ts = 15 and nu = 0.75

ts_test <- 15
nu_test <- -1000
u <- c(rep(param_1$ks, times=ts_test), rep(0, times = param_1$T_end - ts_test))

res <- simple_water_model_sim(param_1$M0, param_1$S0, param_1$W0, param_1$T_end, param_1$kp, param_1$kr, u, param_1$kw)
res$S[param_1$T_end]
ST_func(ts_test,nu_test, param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr, param_1$kw, param_1$T_end)

# What if I now change nu, is that the issue? let's set nu to -1.5

nu_test <- -1.5
ST_func(ts_test,nu_test, param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr, param_1$kw, param_1$T_end)

# tests for nu to be 0.00005, 1200 and negatives

nu_test <- 0.00005
ST_func(ts_test,nu_test, param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr, param_1$kw, param_1$T_end)

nu_test <- 1200000
ST_func(ts_test,nu_test, param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr, param_1$kw, param_1$T_end)

nu_test <- -0.00005
ST_func(ts_test,nu_test, param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr, param_1$kw, param_1$T_end)

nu_test <- -1200000
ST_func(ts_test,nu_test, param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr, param_1$kw, param_1$T_end)

nu_test <- 0
ST_func(ts_test,nu_test, param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr, param_1$kw, param_1$T_end)

c# Is beta plus 1 and 2 equivalent

beta_plus(param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr)
beta_plus_2(param_1$S0, param_1$M0, param_1$ks, param_1$kp, param_1$kr)

#yes they are

nu_test <- 0
tstar_fun_fixed(param_1$kr, param_1$T_end, ts_test, param_1$kp, param_1$kw, nu_test)
tstar_fun_orig(param_1$kr, param_1$T_end, ts_test, param_1$kp, param_1$kw, nu_test)

ts_test <- 20
tstar_fun_fixed(param_1$kr, param_1$T_end, ts_test, param_1$kp, param_1$kw, nu_test)
tstar_fun_orig(param_1$kr, param_1$T_end, ts_test, param_1$kp, param_1$kw, nu_test)

#this gives me 67 so that's really cool. I'll also test tstar and ts

tstar_test <- 67
ts_fun_orig(param_1$kr, param_1$T_end, tstar_test, param_1$kp, param_1$kw, nu_test)
ts_fun_fixed(param_1$kr, param_1$T_end, tstar_test, param_1$kp, param_1$kw, nu_test)

#this works both ways. Good! The fixed version that is.

# Now let's try the other version from my notebook (was not any different)
ts_test <- 20
tstar_fun_fixed(param_1$kr, param_1$T_end, ts_test, param_1$kp, param_1$kw, nu_test)
tstar_fun_orig(param_1$kr, param_1$T_end, ts_test, param_1$kp, param_1$kw, nu_test)
tstar_fun_notebook(param_1$kr, param_1$T_end, ts_test, param_1$kp, param_1$kw, nu_test)




### BELOW IS NOT FINISHED YET 
# 
# #Generate a vector of outputs
# x <- seq(min,max,increment)
# #Find where ST is negative or positive
# y <- sign(ST_func(x))
# #Find where the sign crosses from (+) to (-), i.e. crosses 0. 
# #We want to find that value of mu that gives S_T=0
# plus.to.minus <- which(diff(y)<0)   
# 
# #Find the roots of the function, ie S_T=0
# root <- uniroot(ST_func, interval=c(min,max), tol=.Machine$double.eps)
# mu <- root$root
# ts <- floor(ts_func(mu))
# 
# 
# 
# 
# 
# 
# 
# 
