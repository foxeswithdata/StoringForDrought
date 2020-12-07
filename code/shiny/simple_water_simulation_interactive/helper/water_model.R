simple_water_model_sim <- function(M0, S0, W0, Tend, kp, kr, u, wue){
  # INITIALISE OUTPUT VECTORS
  # browser() 
  M <- rep(M0, times=Tend)
  S <- vector(length=Tend)
  WU <- vector(length=Tend)
  W <- vector(length=Tend)
  A <- vector(length=Tend)
  
  # INITIALISE FIRST VALUE
  M[1] = M0 + S0*u[1]
  A[1] = kp * M0
  WU[1] = A[1] * wue;
  if(WU[1] > W0){
    A[1] = 0;
    WU[1] = 0;
    kp = 0;
  }
  W[1] = W0 - WU[1]
  S[1] = S0 + A[1] - M0 * kr - S0*u[1]
  
  # RUN MODEL FOR REST OF THE PERIOD
  for(i in 2:Tend){
    A[i] = kp * M[i-1]
    WU[i] = A[i] * wue;
    if(WU[i] > W[i-1]){
      WU[i] = W[i-1];
      A[i] = W[i]/wue;
      kp = 0;
    }
    W[i] = W[i-1] - WU[i]
    M[i] = M[i-1] + S[i-1]*u[i]
    S[i] = S[i-1] + A[i] - M[i-1] * kr - S[i-1]*u[i]
    if(S[i] < 0){
      break;
    }
  }
  return(list(S=c(S0,S),M=c(M0, M),A=c(0,A), WU = c(0,WU), W=c(W0,W), t=0:Tend))
}

simple_water_model_sim_time_breaks <- function(M0, S0, W0, Tend, kp, kr, ks, kw, ts, deltat){
  
  # INITIALISE Timeline settings
  t = seq(from=0, to=Tend, by=deltat)
  simlength = length(t)
  # delta t = h
  h = Tend/(simlength-1)
  h2 = h/2
  
  # INITIALISE OUTPUT VECTORS
  M <- vector(length=simlength)
  S <- vector(length=simlength)
  WU <- vector(length=simlength)
  W <- vector(length=simlength)
  A <- vector(length=simlength)
  u <- c(rep(ks, times = sum(t<=ts)), rep(0, times = sum(t>ts)))
  
  # INITIALISE FIRST VALUE
  M[1] = M0
  A[1] = 0
  WU[1] = 0
  W[1] = W0
  S[1] = S0
  
  # RUN MODEL FOR REST OF THE PERIOD
  for(i in 2:simlength){
    if(W[i-1] <=0){
      kp = 0
    }
    k1a = h * kp * M[i-1]
    k1wu = h * kp * M[i-1] * kw
    k1m = u[i-1] * S[i-1]
    k1s = ((kp - kr) * M[i-1] - u[i-1]*S[i-1])
    k1w = - kp * M[i-1] * kw
    
    k2a = (kp * (M[i-1] + h2*k1m))
    k2wu = (kp * (M[i-1] + h2*k1m) * kw)
    k2m = ((0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k1s))
    k2s = ((kp - kr) * (M[i-1] + h2*k1m) - (0.5*(u[i-1] + u[i])) * (S[i-1] + h2*k1s)) # might need to play around a bit with u
    k2w = - (kp * (M[i-1] + h2*k1m) * kw)
    
    k3a = (kp * (M[i-1] + h2*k2m))
    k3wu = (kp * (M[i-1] + h2*k2m) * kw)
    k3m = ((0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k2s))
    k3s = ((kp - kr) * (M[i-1] + h2*k2m) - (0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k2s)) # might need to play around a bit with u
    k3w = - (kp * (M[i-1] + h2*k2m) * kw)
    
    k4a = (kp * (M[i-1] + h*k3m))
    k4wu = (kp * (M[i-1] + h*k3m) * kw)
    k4m = (u[i] * (S[i-1] + h*k3s))
    k4s = ((kp - kr) * (M[i-1] +h*k3m) - u[i] * (S[i-1] + h*k3s)) # might need to play around a bit with u
    k4w = -(kp *(M[i-1] + h*k3m))
    
    A[i] = A[i-1] + h/6 * (k1a + 2 * k2a + 2 * k3a + k4a) 
    WU[i] = WU[i-1] + h/6 * (k1wu + 2 * k2wu + 2 * k3wu + k4wu) 
    M[i] = M[i-1] + h/6 * (k1m + 2 * k2m + 2 * k3m + k4m)
    S[i] = S[i-1] + h/6 * (k1s + 2 * k2s + 2 * k3s + k4s)
    W[i] = W[i-1] + h/6 * (k1w + 2 * k2w + 2 * k3w + k4w)
    
  }
  return(list(S=S,M=M,A=A,WU=WU,W=W,t=t))
}

water_model_log <- function(M0, S0, W0, Tend, kp, kr, ks, kw, ts, deltat, w_mid, growth_rate){
  
  # INITIALISE Timeline settings
  t = seq(from=0, to=Tend, by=deltat)
  simlength = length(t)
  # delta t = h
  h = Tend/(simlength-1)
  h2 = h/2
  
  # INITIALISE OUTPUT VECTORS
  M <- vector(length=simlength)
  S <- vector(length=simlength)
  WU <- vector(length=simlength)
  W <- vector(length=simlength)
  A <- vector(length=simlength)
  u <- c(rep(ks, times = sum(t<=ts)), rep(0, times = sum(t>ts)))
  
  # INITIALISE FIRST VALUE
  M[1] = M0
  A[1] = 0
  WU[1] = 0
  W[1] = W0
  S[1] = S0
  
  # RUN MODEL FOR REST OF THE PERIOD
  for(i in 2:simlength){
    k1a = h * (kp / 1 + exp(-growth_rate * (W[i-1] - w_mid)))  * M[i-1]
    k1wu = h * (kp / 1 + exp(-growth_rate * (W[i-1] - w_mid))) * M[i-1] * kw
    k1m = u[i-1] * S[i-1]
    k1s = (((kp / 1 + exp(-growth_rate * (W[i-1] - w_mid))) - kr) * M[i-1] - u[i-1]*S[i-1])
    k1w = - (kp / 1 + exp(-growth_rate * (W[i-1] - w_mid))) * M[i-1] * kw
    
    k2a = ((kp / 1 + exp(-growth_rate * ((W[i-1] + h2 * k1w) - w_mid))) * (M[i-1] + h2*k1m))
    k2wu = ((kp / 1 + exp(-growth_rate * ((W[i-1] + h2 * k1w) - w_mid))) * (M[i-1] + h2*k1m) * kw)
    k2m = ((0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k1s))
    k2s = (((kp / 1 + exp(-growth_rate * ((W[i-1] + h2 * k1w) - w_mid))) - kr) * (M[i-1] + h2*k1m) - (0.5*(u[i-1] + u[i])) * (S[i-1] + h2*k1s)) # might need to play around a bit with u
    k2w = - ((kp / 1 + exp(-growth_rate * ((W[i-1] + h2 * k1w) - w_mid))) * (M[i-1] + h2*k1m) * kw)
    
    k3a = ((kp / 1 + exp(-growth_rate * ((W[i-1] + h2 * k2w) - w_mid))) * (M[i-1] + h2*k2m))
    k3wu = ((kp / 1 + exp(-growth_rate * ((W[i-1] + h2 * k2w) - w_mid))) * (M[i-1] + h2*k2m) * kw)
    k3m = ((0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k2s))
    k3s = (((kp / 1 + exp(-growth_rate * ((W[i-1] + h2 * k2w) - w_mid))) - kr) * (M[i-1] + h2*k2m) - (0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k2s)) # might need to play around a bit with u
    k3w = - ((kp / 1 + exp(-growth_rate * ((W[i-1] + h2 * k2w) - w_mid))) * (M[i-1] + h2*k2m) * kw)
    
    k4a = ((kp / 1 + exp(-growth_rate * ((W[i-1] + h * k3w) - w_mid))) * (M[i-1] + h*k3m))
    k4wu = ((kp / 1 + exp(-growth_rate * ((W[i-1] + h * k3w) - w_mid))) * (M[i-1] + h*k3m) * kw)
    k4m = (u[i] * (S[i-1] + h*k3s))
    k4s = (((kp / 1 + exp(-growth_rate * ((W[i-1] + h * k3w) - w_mid))) - kr) * (M[i-1] +h*k3m) - u[i] * (S[i-1] + h*k3s)) # might need to play around a bit with u
    k4w = -((kp / 1 + exp(-growth_rate * ((W[i-1] + h * k3w) - w_mid))) *(M[i-1] + h*k3m))
    
    A[i] = A[i-1] + h/6 * (k1a + 2 * k2a + 2 * k3a + k4a) 
    WU[i] = WU[i-1] + h/6 * (k1wu + 2 * k2wu + 2 * k3wu + k4wu) 
    M[i] = M[i-1] + h/6 * (k1m + 2 * k2m + 2 * k3m + k4m)
    S[i] = S[i-1] + h/6 * (k1s + 2 * k2s + 2 * k3s + k4s)
    W[i] = W[i-1] + h/6 * (k1w + 2 * k2w + 2 * k3w + k4w)
    
  }
  return(list(S=S,M=M,A=A,WU=WU,W=W,t=t))
}




