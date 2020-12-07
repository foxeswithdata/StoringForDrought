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

beta_minus_2 <- function(S0, M0, ks, kp, kr){
  b_m = (S0*ks - M0 * alpha_plus(ks, kp, kr))/(alpha_minus(ks, kp, kr)- alpha_plus(ks, kp, kr))
  return(b_m)
}

beta_plus_2 <- function(S0, M0, ks, kp, kr){
  b_p = M0 - (S0 * ks - M0 * alpha_plus(ks, kp, kr))/(alpha_minus(ks, kp, kr)- alpha_plus(ks, kp, kr))
  return(b_p)
}

ts_fun <- function(kr, T_end, t_star, kp, kw, mu){
  ts <- (1 + kr * T_end - (1 + mu * kw * (t_star - T_end)) * kp * t_star)/(kr - (1 + mu * kw * (t_star - T_end)) * kp)
  return(ts)
}

ST_func <- function(tstar, mu, S0, M0, ks, kp, kr, kw, T_end){
  ap <- alpha_plus(ks, kp, kr)
  am <- alpha_minus(ks, kp, kr)
  bp <- beta_plus(S0, M0, ks, kp, kr)
  bm <- beta_minus(S0, M0, ks, kp, kr)
  ts <- ts_fun(kr, T_end, tstar, kp, kw, mu)
  SI_t_s = 1/ks * (bp * ap * exp(ap * ts) + bm * am * exp(am * ts))
  MI_t_s = bp * exp(ap * ts) + bm * exp(am * ts)
  SII_t_star = SI_t_s + (kp - kr) * MI_t_s * (tstar - ts)
  SIII_T_end = SII_t_star + kr * MI_t_s * (tstar - T_end)
  return(SIII_T_end)
}

Wt_star_func <- function(tstar, mu, W0, S0, M0, ks, kp, kr, kw, T_end){
  ap <- alpha_plus(ks, kp, kr)
  am <- alpha_minus(ks, kp, kr)
  bp <- beta_plus(S0, M0, ks, kp, kr)
  bm <- beta_minus(S0, M0, ks, kp, kr)
  ts <- ts_fun(kr, T_end, tstar, kp, kw, mu)
  
  WII_star = W0 + kp * kw * ( bp/ap * (1-exp(ap * tstar)) - bp * exp(ap * ts) * (tstar - ts) + bm / am * (1 - exp(am * tstar)) - bm * exp(am * ts) * (tstar - ts) )

  return(WII_star);
}

source('code/parameters.R')
param <- param_list_water_short()

# mu must be positive
mu_ar <- seq(0, 10, 0.05);
tstar <- seq(0, param$Tend, 0.5)

Wtstar <- matrix(nrow = length(mu_ar), ncol=length(tstar))
ST_res <- matrix(nrow = length(mu_ar), ncol=length(tstar))

l = 1;
for(i in mu_ar){
  k = 1;
  for(j in tstar){
    Wtstar[l,k] = Wt_star_func(j,i, param$W0, param$S0, param$M0, param$ks, param$kp, param$kr, param$kw, param$Tend)
    ST_res[l,k] = ST_func(j, i, param$S0, param$M0, param$ks, param$kp, param$kr, param$kw, param$Tend)
    k=k+1
  }
  l = l+1;
}


(beta_plus(param$S0, param$M0, param$ks, param$kp, param$kr)/alpha_plus(param$ks, param$kp, param$kr) + 
  beta_minus(param$S0, param$M0, param$ks, param$kp, param$kr)/alpha_minus(param$ks, param$kp, param$kr)) * param$kp




