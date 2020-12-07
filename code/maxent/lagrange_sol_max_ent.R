find_lambda_err <- function(M0,S0,Tend,kp,kr,kmax,lambda,Itot){
  
  Z = NA;
  x = NA;
  
  for (i in 1:Tend) {
    u = c(rep(kmax,times=i), rep(0, times=Tend-i))
    out <- simplemodel_sim(M0, S0, Tend, kp, kr, u);
    E = sum(out$A)
    e_val = exp(lambda*E)
    if(i == 1){
      x = lambda*E + log(E)
      Z = lambda*E
    }
    else{
      x = x + log(1 + exp(lambda*E + log(E) - x))
      Z = Z + log(1 + exp(lambda*E - Z))
    }
    
  }
  
  result <- exp(x-Z) - Itot
  return(result)
}
pk_res <- function(M0, S0, Tend, kp, kr, kmax, ts, lambda){
  
  Z = NA;
  
  for (i in 1:Tend) {
    u = c(rep(kmax,times=i), rep(0, times=Tend-i))
    out <- simplemodel_sim(M0, S0, Tend, kp, kr, u);
    E = sum(out$A)
    if(i == 1){
      Z = lambda*E
    }
    else{
      Z = Z + log(1 + exp(lambda*E - Z))
    }
    
  }
  
  u = c(rep(kmax,times=ts), rep(0, times=Tend-ts))
  out <- simplemodel_sim(M0, S0, Tend, kp, kr, u);
  E = sum(out$A)
  
  pk_log <- lambda * E - Z
  pk <- exp(pk_log)
  
  return(pk)
  
}

find_lambda_err_alive <- function(M0,S0,Tend,kp,kr,kmax,max_ts,lambda,Itot){
  
  Z = NA;
  x = NA;
  
  for (i in 1:max_ts) {
    u = c(rep(kmax,times=i), rep(0, times=Tend-i))
    out <- simplemodel_sim(M0, S0, Tend, kp, kr, u);
    E = sum(out$A)
    e_val = exp(lambda*E)
    if(i == 1){
      x = lambda*E + log(E)
      Z = lambda*E
    }
    else{
      x = x + log(1 + exp(lambda*E + log(E) - x))
      Z = Z + log(1 + exp(lambda*E - Z))
    }
    
  }
  
  # print(exp(x-Z))
  # print(Itot)
  result <- exp(x-Z) - Itot
  # print(result)
  return(result)
}

pk_res_alive <- function(M0, S0, Tend, kp, kr, kmax, maxts, ts, lambda){
  
  Z = NA;
  
  for (i in 1:maxts) {
    u = c(rep(kmax,times=i), rep(0, times=Tend-i))
    out <- simplemodel_sim(M0, S0, Tend, kp, kr, u);
    E = sum(out$A)
    if(i == 1){
      Z = lambda*E
    }
    else{
      Z = Z + log(1 + exp(lambda*E - Z))
    }
    
  }
  
  u = c(rep(kmax,times=ts), rep(0, times=Tend-ts))
  out <- simplemodel_sim(M0, S0, Tend, kp, kr, u);
  E = sum(out$A)
  # E = sum(out$A)
  
  pk_log <- lambda * E - Z
  pk <- exp(pk_log)
  
  return(pk)
  
}
