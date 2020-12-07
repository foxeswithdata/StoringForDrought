simple_model_sim_time_breaks <- function(M0, S0, Tend, kp_star, kr, ks, ts, tcrit, deltat){
  
  # INITIALISE Timeline settings
  t = seq(from=0, to=Tend, by=deltat)
  simlength = length(t)
  # delta t = h
  h = Tend/(simlength-1)
  h2 = h/2
  
  # INITIALISE OUTPUT VECTORS
  M <- vector(length=simlength)
  S <- vector(length=simlength)
  A <- vector(length=simlength)
  u <- c(rep(ks, times = sum(t<=ts)), rep(0, times = sum(t>ts)))
  kp <- c(rep(kp_star, times = sum(t<=tcrit)), rep(0, times = sum(t>tcrit)))
  
  # INITIALISE FIRST VALUE
  M[1] = M0
  A[1] = 0
  S[1] = S0
  
  # RUN MODEL FOR REST OF THE PERIOD
  for(i in 2:simlength){
    k1a = kp[i-1] * M[i-1]
    k1m = u[i-1] * S[i-1]
    k1s = ((kp[i-1] - kr) * M[i-1] - u[i-1]*S[i-1])
    
    k2a = (0.5 * (kp[i-1] + kp[i]) * (M[i-1] + h2*k1m))
    k2m = ((0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k1s))
    k2s = ((0.5 * (kp[i-1] + kp[i]) - kr) * (M[i-1] + h2*k1m) - (0.5*(u[i-1] + u[i])) * (S[i-1] + h2*k1s)) # might need to play around a bit with u
    
    k3a = (0.5 * (kp[i-1] + kp[i]) * (M[i-1] + h2*k2m))
    k3m = ((0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k2s))
    k3s = ((0.5 * (kp[i-1] + kp[i]) - kr) * (M[i-1] + h2*k2m) - (0.5 * (u[i-1] + u[i])) * (S[i-1] + h2*k2s)) # might need to play around a bit with u
    
    k4a = (kp[i] * (M[i-1] + h*k3m))
    k4m = (u[i] * (S[i-1] + h*k3s))
    k4s = ((kp[i] - kr) * (M[i-1] +h*k3m) - u[i] * (S[i-1] + h*k3s)) # might need to play around a bit with u
    
    A[i] = A[i-1] + h/6 * (k1a + 2 * k2a + 2 * k3a + k4a) 
    M[i] = M[i-1] + h/6 * (k1m + 2 * k2m + 2 * k3m + k4m)
    S[i] = S[i-1] + h/6 * (k1s + 2 * k2s + 2 * k3s + k4s)
    
  }
  return(list(S=S,M=M,A=A,t=t,u=u,kp=kp))
}

simplemodel_sim <- function(M0, S0, Tend, kp, kr, u){
  # INITIALISE OUTPUT VECTORS
  M <- vector(length=Tend)
  S <- vector(length=Tend)
  A <- vector(length=Tend)
  
  # INITIALISE FIRST VALUE
  M[1] = M0 + S0*u[1]
  A[1] = kp[1] * M0
  S[1] = S0 + A[1] - M0 * kr - S0*u[1]
  
  # RUN MODEL FOR REST OF THE PERIOD
  for(i in 2:Tend){
    M[i] = M[i-1] + S[i-1]*u[i]
    A[i] = kp[i] * M[i-1]
    S[i] = S[i-1] + A[i] - M[i-1] * kr - S[i-1]*u[i]
  }
  return(list(S=S,M=M,A=A))
}

simpleModel_1 <- function(M0,S0,Tend,tcrit,ts,kp,kr,ks){
  
  #Define output vectors: Mass M, Photosynthesis A, Storage S
  M <- vector(length=Tend)
  S <- vector(length=Tend)
  A <- vector(length=Tend)
  
  #First period of time 0 < t <= ts
  #Period of plenty, allocate to growth
  for(i in 1:ts){
    if(i == 1){
      S[i]=S0 + (kp-kr)*M0-ks*S0
      M[i]=M0 + (ks*S0)
      A[i]=kp*M0;
    }
    else{
      S[i]=S[i-1] + (kp-kr)*M[i-1] - ks* S[i-1];
      M[i]=M[i-1] + (ks* S[i-1]);
      A[i]=kp*M[i-1];
    }
  }
  #Second period ts < t <= tcrit
  #Period of plenty, allocate to storage
  for(i in (ts+1):tcrit){
    S[i]=S[i-1] + (kp-kr)*M[i-1];
    M[i]=M[i-1];
    A[i]=kp*M[i-1];
  }
  #Third period tcrit < t <= T
  #Period of drought, no photosynthesis
  for(i in (tcrit+1):T_end){
    S[i]=S[i-1] - kr*M[i-1]
    M[i]=M[i-1]
  }
  #Return all three components
  return(list(S=S,M=M,A=A))
}

#Find value of biomass in the first stage of the problem
M_t_1_func <- function(t, beta_p, beta_m, alpha_p, alpha_m){
  return(beta_p*exp(alpha_p*t)+beta_m*exp(alpha_m*t))
}

find_ts_SimpleModel1 <- function(M0,S0,Tend,tcrit,kp,kr,ks,min,max,increment){
  alpha_p <- (-ks+sqrt(ks^2+4*ks*(kp-kr)))/2
  alpha_m <- (-ks-sqrt(ks^2+4*ks*(kp-kr)))/2
  beta_m <- (S0*ks-M0*alpha_p)/(alpha_m-alpha_p)
  beta_p <- M0-beta_m
  
  #Find S_ts
  S_t_1_func <- function(t){
    return((beta_p*alpha_p*exp(alpha_p*t)+beta_m*alpha_m*exp(alpha_m*t))/ks)
  }
  
  #Find ts
  ts_func <- function(nu){
    return((nu + kr * nu * T_end - kp * tcrit * (1+nu))/(-kp-(kp-kr)*nu))
  }
  
  #Find S_T 
  ST_func <- function(mu){
    ts <- ts_func(mu)
    Mts <- M_t_1_func(ts, beta_p, beta_m, alpha_p, alpha_m)
    Sts <- S_t_1_func(ts)
    return(kr*Mts*(tcrit-T_end)+(kp-kr)*Mts*(tcrit-ts)+Sts)
  }
  
  #Generate a vector of outputs
  x <- seq(min,max,increment)
  #Find where ST is negative or positive
  y <- sign(ST_func(x))
  #Find where the sign crosses from (+) to (-), i.e. crosses 0. 
  #We want to find that value of mu that gives S_T=0
  plus.to.minus <- which(diff(y)<0)   
  
  #Find the roots of the function, ie S_T=0
  root <- uniroot(ST_func, interval=c(min,max), tol=.Machine$double.eps)
  mu <- root$root
  ts <- floor(ts_func(mu))
  return(ts)
}

#Find the total photosynthesis for a given model
Atotal <- function(M0,S0,tcrit,ts,kp,kr,ks){
  alpha_p <- (-ks+sqrt(ks^2+4*ks*(kp-kr)))/2
  alpha_m <- (-ks-sqrt(ks^2+4*ks*(kp-kr)))/2
  beta_m <- (S0*ks-M0*alpha_p)/(alpha_m-alpha_p)
  beta_p <- M0-beta_m
  A = 0
  for (i in 1:ts){
    A = A + M_t_1_func(i, beta_p, beta_m, alpha_p, alpha_m)*kp;
  }
  A = A + M_t_1_func(ts, beta_p, beta_m, alpha_p, alpha_m)*kp*(tcrit-ts)
}


