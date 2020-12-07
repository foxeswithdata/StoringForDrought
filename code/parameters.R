param_list_simple_early <- function(){
  param_1 <- list(S0 = 0,
                  M0 = 0.1,
                  Tend = 365,
                  kp = 1,
                  ks = 0.1,
                  kr = 0.5,
                  tcrit = 200)
  return(param_1)
}

param_list_simple_short <- function(){
  param_1 <- list(S0 = 10,
                  M0 = 100,
                  W0 = 200,
                  Tend = 64,
                  kp = 0.04,
                  ks = 0.1,
                  kr = 0.02,
                  kw = 1,
                  tcrit = 24)
  return(param_1)
}


param_list_water_short <- function(){
  param_2 <- list(S0 = 10,
                  M0 = 100,
                  W0 = 200,
                  Tend = 64,
                  kp = 0.04,
                  ks = 0.1,
                  kr = 0.02,
                  kw = 1)
  return(param_2)
}


param_list_water_short2 <- function(){
  param_2 <- list(S0 = 10,
                  M0 = 100,
                  W0 = 300,
                  Tend = 64,
                  kp = 0.04,
                  ks = 0.1,
                  kr = 0.02,
                  kw = 1,
                  tcrit = 45)
  return(param_2)
}

param_list_water_short2 <- function(){
  param <- list(S0 = 10,
                  M0 = 100,
                  W0 = 200,
                  Tend = 64,
                  kp = 0.04,
                  ks = 0.1,
                  kr = 0.02,
                  kw = 1,
                  krh = 5,
                  klog = 10,
                  w50 = 2,
                  tcrit = 45)
  param$x0 <- c(param$M0, param$S0, param$W0)
  return(param)
}

param_list_drake_W0_1000 <- function(){
  param <- list(S0 = 550,
                M0 = 3300,
                W0 = 1000,
                Tend = 150,
                kp = 0.014,
                ks = 0.06,
                kr = 0.004,
                kw = 0.4,
                krh = 5,
                klog = 10,
                w50 = 2,
                tcrit = 75)
  param$x0 <- c(param$M0, param$S0, param$W0)
  return(param)
}

param_list_drake_W0_3000 <- function(){
  param <- list(S0 = 550,
                M0 = 3300,
                W0 = 3000,
                Tend = 150,
                kp = 0.014,
                ks = 0.06,
                kr = 0.004,
                kw = 0.4,
                krh = 5,
                klog = 10,
                w50 = 2,
                tcrit = 75)
  param$x0 <- c(param$M0, param$S0, param$W0)
  return(param)
}


