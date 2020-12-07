rm(list=ls())

source("code/parameters.R")
source("code/models/water_model.R")
params <- param_list_drake_W0_1000()

r1 <- function(params, ut){
   return((-ut+sqrt(ut^2 + 4 * ut * (params$kp-params$kr)))/2)
}

r2 <- function(params, ut){
   return((-ut-sqrt(ut^2 - 4 * ut * (params$kp-params$kr)))/2)
}

C1 <- function(params, ut){
   r1_ = r1(params, ut)
   r2_ = r2(params, ut)
   return((params$M0 - (ut * params$S0 - params$M0 * r1_)/(r1_ - r2_)))
}

C2 <- function(params, ut){
   r1_ = r1(params, ut)
   r2_ = r2(params, ut)
   return((ut * params$S0 - params$M0 * r1_)/(r1_ - r2_))
}

Mt_1 <- function(params, t){
   C1_ = C1(params, params$ks)
   C2_ = C2(params, params$ks)
   r1_ = r1(params, params$ks)
   r2_ = r2(params, params$ks)
   
   out = C1_ * exp(r1_ * t) + C2_ * exp(r2_ * t)
   
   return(out)
}

St_1 <- function(params, t){
   C1_ = C1(params, params$ks)
   C2_ = C2(params, params$ks)
   r1_ = r1(params, params$ks)
   r2_ = r2(params, params$ks)
   Mt_ = Mt_1(params, t)
   
   out = (C1_ * r1_ * exp(r1_ * t) + C2_ * r2_ * exp(r2_ * t))/params$ks
   
   return(out)
}

St_2 <- function(params, ts, t){
   Mt_ = Mt_1(params, ts)
   St1_ts = St_1(params, ts)
   
   
   out = St1_ts + (t-ts) * (params$kp - params$kr) * Mt_
   
   return(out)
}

St_3 <- function(params, tcrit, ts, t){
   Mt_ = Mt_1(params, ts)
   St2_ts = St_2(params, ts, tcrit)
   
   out = St2_ts + (t-tcrit) * (- params$kr) * Mt_
   
   return(out)
}

Wt_1 <- function(params, t){
   C1_ = C1(params, params$ks)
   C2_ = C2(params, params$ks)
   r1_ = r1(params, params$ks)
   r2_ = r2(params, params$ks)
   
   out = params$W0 + params$kp * params$kw * (C1_/r1_ + C2_/r2_ - 
                                                 C1_/r1_ * exp(r1_ * t) - 
                                                 C2_/r2_ * exp(r2_ * t) )
   
   return(out)
}

Wt_2 <- function(params, ts, t){
   Mt_ = Mt_1(params, ts)
   Wt1_ts = Wt_1(params, ts)
   
   out = Wt1_ts = Wt_1(params, ts) - (t-ts) * Mt_ * params$kp * params$kw
   
   return(out)
}

tcrit_given_ts <- function(params, ts){

   f <- function(x){
      Wt_2(params, ts, x)
   }
   
   #Find the roots of the function, ie S_T=0
   root <- uniroot(f, interval=c(0,params$Tend), tol=.Machine$double.eps)
   
   return(root)
}








ts = 20
out_sim <- simple_water_model_sim_time_breaks(params, ts, deltat=0.1)

t_test <- which(out_sim$t == 0.5)

St_1_test <- St_1(params, 0.5)
Mt_1_test <- Mt_1(params, 0.5)
Wt_1_test <- Wt_1(params, 0.5)

out_sim$W[t_test]
out_sim$M[t_test]
out_sim$S[t_test]

St_1_test
Wt_1_test
Mt_1_test

out_sim$W[t_test] - Wt_1_test
out_sim$S[t_test] - St_1_test
out_sim$M[t_test] - Mt_1_test


tcrit_given_ts(params, 20)

plot(out_sim$t, out_sim$W, type="l")

