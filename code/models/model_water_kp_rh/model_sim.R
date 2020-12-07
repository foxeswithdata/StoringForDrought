#' Find the trajectory of allocation to storage and biomass 
#' given a switch point
#' model version R-H-water-model
#'  
#' @description 
#' This function produces a simulation of a model run given a time point 
#' ts at which the allocation switches from growth to storage.
#' 
#' This simulation uses runge-kutta expansion to approximate the continuous 
#' functions. 
#' 
#' Full model description can be found under:
#' /documentation/models/r-h-water-model.html
#' 
#' @param param a list of model parameters
#' @param ts switch time point
#' @param deltat spacing between (time) nodes
#' 
#' @return list of vectors containing: time (sequence of time), S, M, W (sequences of states M, S, W), 
#'         A and WU (sequences of processes photosynthesis A and water use WU),
#'         u (sequence of control results)
simulation_kp_rh <- function(param, ts, deltat=0.1){
  
  # INITIALISE Timeline settings
  t = seq(from=0, to=param$Tend, by=deltat)
  simlength = length(t)
  # delta t = h
  h = param$Tend/(simlength-1)
  h2 = h/2
  
  # INITIALISE OUTPUT VECTORS
  M <- vector(length=simlength)
  S <- vector(length=simlength)
  WU <- vector(length=simlength)
  W <- vector(length=simlength)
  A <- vector(length=simlength)
  u <- c(rep(param$ks, times = sum(t<=ts)), rep(0, times = sum(t>ts)))
  
  # INITIALISE FIRST VALUE
  M[1] = param$M0
  A[1] = 0
  WU[1] = 0
  W[1] = param$W0
  S[1] = param$S0
  
  # RUN MODEL FOR REST OF THE PERIOD
  for(i in 2:simlength){
    k1p = param$kp * W[i-1] /(param$krh + W[i-1])
    if(W[i-1] <= 0){
      k1p = 0
    }
    k1a = k1p  * M[i-1]
    k1wu = k1p * M[i-1] * param$kw
    k1m = u[i-1] * S[i-1]
    k1s = ((k1p - param$kr) * M[i-1] - u[i-1]*S[i-1])
    k1w = - k1p * M[i-1] * param$kw
    
    k2p = param$kp * (W[i-1] + h2 * k1w)/ (param$krh + (W[i-1] + h2 * k1w))
    if((W[i - 1] + h/6 * k1w) <= 0){
      k2p = 0
    }
    k2a = (k2p * (M[i-1] + h2*k1m))
    k2wu = (k2p * (M[i-1] + h2*k1m) * param$kw)
    k2m = ((0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k1s))
    k2s = ((k2p - param$kr) * (M[i-1] + h2*k1m) - (0.5*(u[i-1] + u[i])) * (S[i-1] + h2*k1s)) # might need to play around a bit with u
    k2w = - (k2p * (M[i-1] + h2*k1m) * param$kw)
    
    k3p = param$kp * (W[i-1] + h2 * k2w)/ (param$krh + (W[i-1] + h2 * k2w))
    if((W[i - 1] + h/6 * (k1w + 2 * k2w)) <= 0){
      k3p = 0
    }
    k3a = (k3p * (M[i-1] + h2*k2m))
    k3wu = (k3p * (M[i-1] + h2*k2m) * param$kw)
    k3m = ((0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k2s))
    k3s = ((k3p - param$kr) * (M[i-1] + h2*k2m) - (0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k2s)) # might need to play around a bit with u
    k3w = - (k3p * (M[i-1] + h2*k2m) * param$kw)
    
    k4p = param$kp * (W[i-1] + h * k3w)/ (param$krh + (W[i-1] + h * k3w))
    if((W[i - 1] + h/6 * (k1w + 2 * k2w + 2 * k3w)) <= 0){
      k4p = 0
    }
    k4a = (k4p * (M[i-1] + h*k3m))
    k4wu = (k4p * (M[i-1] + h*k3m) * param$kw)
    k4m = (u[i] * (S[i-1] + h*k3s))
    k4s = ((k4p - param$kr) * (M[i-1] +h*k3m) - u[i] * (S[i-1] + h*k3s)) # might need to play around a bit with u
    k4w = -(k4p *(M[i-1] + h*k3m) * param$kw)
    
    A[i] = A[i-1] + h/6 * (k1a + 2 * k2a + 2 * k3a + k4a) 
    WU[i] = WU[i-1] + h/6 * (k1wu + 2 * k2wu + 2 * k3wu + k4wu) 
    M[i] = M[i-1] + h/6 * (k1m + 2 * k2m + 2 * k3m + k4m)
    S[i] = S[i-1] + h/6 * (k1s + 2 * k2s + 2 * k3s + k4s)
    W[i] = W[i-1] + h/6 * (k1w + 2 * k2w + 2 * k3w + k4w)
    
  }
  return(list(S=S,M=M,A=A,WU=WU,W=W,t=t, u=u))
}

simulation_kp_rh_simple <- function(param, ts, deltat=0.1){
  
  # INITIALISE Timeline settings
  t = seq(from=0, to=param$Tend, by=deltat)
  simlength = length(t)
  
  # INITIALISE OUTPUT VECTORS
  M <- vector(length=simlength)
  S <- vector(length=simlength)
  WU <- vector(length=simlength)
  W <- vector(length=simlength)
  A <- vector(length=simlength)
  u <- c(rep(param$ks, times = sum(t<=ts)), rep(0, times = sum(t>ts)))
  
  # INITIALISE FIRST VALUE
  M[1] = param$M0
  A[1] = 0
  WU[1] = 0
  W[1] = param$W0
  S[1] = param$S0
  
  # RUN MODEL FOR REST OF THE PERIOD
  for(i in 2:simlength){
    kp = param$kp * W[i-1] /(param$krh + W[i-1])
    if(W[i-1] < 0){
      kp = 0
    }
    A[i] = kp  * M[i-1] * deltat
    WU[i] = kp * M[i-1] * param$kw * deltat
    M[i] = M[i-1] + u[i-1] * S[i-1] * deltat
    S[i] = S[i-1] + ((kp - param$kr) * M[i-1] - u[i-1]*S[i-1]) * deltat
    W[i] = W[i-1] - kp * M[i-1] * param$kw * deltat
  }
  return(list(S=S,M=M,A=A,WU=WU,W=W,t=t, u=u))
}
