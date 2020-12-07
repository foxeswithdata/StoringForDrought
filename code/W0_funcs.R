W0min <- function(param){
  return(param$kw * (param$kr * param$M0 * param$Tend - param$S0))
}

W0max <- function(param, water_func, W0range, deltaW0, deltat, kf){

  W0 = seq(from=W0range[1], to=W0range[2], by=deltaW0);
  # Initialise a vector of possible ts values and the length of this
  # vector
  t_seq <- seq(from=0, to=param$Tend, by=deltat)
  simlength <- length(t_seq)
  
  ts <- vector();
  res <- vector();
 
  j = 1
  for (w0 in W0){
    S_out <- vector(length=simlength)
    M_out <- vector(length=simlength)
    
    k=1
    for (i in t_seq){
      # Find the output values of simulation for the specific ts value
      param_2 <- param;
      param_2$W0 <- w0;
      param_2$x0 <- c(param_2$M0, param_2$S0, param_2$W0)
      
      out <- water_func(param_2, i, deltat)
      
      # Save the values of S[Tend], M[Tend], W[Tend] and total photosynthesis A
      S_out[k] <- out$S[simlength]
      M_out[k] <- out$M[simlength]
      
      if(sum(out$S<0) > 0){
        S_out[k] <- NA
        M_out[k] <- NA
      }
      k = k+1
    }
    
    ms_sum <- colSums(rbind (kf * M_out, (1 - kf) * S_out))
    ts_ind <- which.max(ms_sum)
    if(length(ts_ind) == 0){
      ts[j] = NA
      res[j] = NA
    }
    else{
      ts[j] <- t_seq[ts_ind]
      res[j] <- ms_sum[ts_ind]
    }
    j = j+1
  }
  
  W0_max <- W0[which.max(res)]
  return(W0_max)
}

W0maxM <- function(param){
  W0MaxM = param$kw * param$kp *((beta_p(param)/alpha_p(param))*exp(alpha_p(param)*param$Tend) + (beta_m(param)/alpha_m(param))*exp(alpha_m(param)*param$Tend))
  return(W0MaxM)
}

alpha_p <- function(param){
  return((-param$ks + sqrt(param$ks^2 + 4 * param$ks * (param$kp - param$kr)))/2)
}

alpha_m <- function(param){
  return((-param$ks - sqrt(param$ks^2 + 4 * param$ks * (param$kp - param$kr)))/2)
}

beta_p <- function(param){
  return(param$M0-(param$S0 * param$ks - param$M0 * alpha_p(param))/(alpha_m(param)-alpha_p(param)))
}

beta_m <- function(param){
  return((param$S0 * param$ks - param$M0 * alpha_p(param))/(alpha_m(param)-alpha_p(param)))
}


