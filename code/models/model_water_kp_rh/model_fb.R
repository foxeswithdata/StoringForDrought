#' Find the optimal trajectory of allocation to storage and biomass 
#' uses forward-backward algorithm 
#' model version R-H-water-model
#'  
#' @description 
#' This model uses the forward-backward algorithm and runge-kutta 
#' expansion to solve an optimal control problem. 
#' 
#' The forward-backward algorithm operates by simulating the 
#' model (states) forwards and corresponding adjuncts backwards.
#' The model then converges by comparing previous solutions with 
#' new ones. 
#' 
#' Full algorithm description can be found in:
#' [insert oct book name]
#' 
#' Full model description can be found under:
#' /documentation/models/r-h-water-model.html
#'0
#' @param param list of all the parameters for the model
#' @param x0 vector containing initial state values:
#'           x0[1] - M0, x0[2] - S0, x[3] - W0
#' @param delta_t spacing between nodes
#' @param mu_w lagrange multiplier responsible for the water constraint Wt>0
#' @param maxiter maximum number of iterations of the forward-backward loop
#' 
#' @return list containing vectors time (sequence of time), x (sequence of states M, S, W), 
#'         lambda (sequence of lambda values equivalent to each state), u (sequence of control results) 
forwardback_kp_rh <- function(param,
                                 x0,
                                 delta_t,
                                 mu_w,
                                 maxiter = 100000000) {
  iter = 0;
  test = -1; #this is to find convergence
      
  delta = 0.001; # convergence value
  t = seq(from=0,to=param$Tend,by=delta_t); # length of simulation 
  N = length(t)-1 #length of simulation
  h = param$Tend/N; # spacing between nodes
  h2 = h/2; # used for runge-kutta subroutine
  
  
  u = rep(sum(param$ubound)/2,times=N+1); # this is essentially ks
  
  x <- matrix(0, 3, N+1) # initialise guess for state variables
  
  x[,1] = x0; # substitute the first element for the initial state values
  
  lambda = matrix(1,3,N+1); # initialise guess for the adjoint
  
  #The following are from the analytical solution
  lambda[1,] = lambda[1,] * 0;
  lambda[2,] = lambda[2,] * (-1);
  lambda[3,] = lambda[3,] * 0;
  
  while(test < 0){
    iter = iter+1;
    
    oldu = u;
    # oldkp = kp;
    oldx = x;
    oldlambda = lambda;
    
    
    # The following is runge kutta as an approximation of the
    # continuous time - need to read up more about it
    
    for(i in 1:N){
      # CALCULATE XM (x(1,:)) XS (x(2,:)) and XW (x(3,:))  
      
      kp1 = param$kp * x[3,i] /(param$k + x[3,i])
      if(x[3,i] <= 0){
        kp1 = 0
      }
      km1 = x[2,i] * u[i];
      ks1 = (kp1 - param$kr) * x[1,i] - x[2,i] * u[i];
      kw1 = -param$kw * kp1 * x[1,i];
      
      kp2 = param$kp * (x[3,i] + h2 * kw1)/ (param$k + (x[3,i] + h2 * kw1))
      if(x[3,i] <= 0){
        kp2 = 0
      }
      km2 = (x[2,i] + h2*ks1) * 0.5*(u[i] + u[i+1]);
      ks2 = (kp2 - param$kr) * (x[1,i] + h2*km1) - (x[2,i] + h2*ks1) * 0.5*(u[i] + u[i+1]);
      kw2 = -kw * kp2 * (x[1,i] + h2 * km1);
      
      kp3 = param$kp * (x[3,i] + h2 * kw2)/ (param$k + (x[3,i] + h2 * kw2))
      if(x[3,i] <= 0){
        kp3 = 0
      }
      km3 = (x[2,i] + h2*ks2) * 0.5*(u[i] + u[i+1]);
      ks3 = (kp3 - param$kr) * (x[1,i] + h2*km2) - (x[2,i] + h2*ks2) * 0.5*(u[i] + u[i+1]);
      kw3 = -param$kw * kp3 * (x[1,i] + h2*km2);
      
      kp4 = param$kp * (x[3,i] + h * kw3)/ (param$k + (x[3,i] + h * kw3))
      if(x[3,i] <= 0){
        kp4 = 0
      }
      km4 = (x[2,i] + h*ks3) * u[i+1];
      ks4 = (kp4 - param$kr) * (x[1,i] + h*km3) - (x[2,i] + h*ks3) * u[i+1];
      kw4 = -param$kw * kp4 * (x[1,i] + h * km3);
      
      x[1,i+1] = x[1,i] + (h/6)*(km1+2*km2+2*km3+km4);
      x[2,i+1] = x[2,i] + (h/6)*(ks1+2*ks2+2*ks3+ks4);
      x[3,i+1] = x[3,i] + (h/6)*(kw1+2*kw2+2*kw3+kw4);
    }
    
    
    # Runge-Kutta sweep backwards to find lambda
    for(i in 1:N){
      j = N + 2 - i;
      
      kp1 = param$kp * x[3,j]/ (param$k + x[3,j])
      if(x[3,j] <= 0){
        kp1 = 0
      }
      km1 = -lambda[2,j] * (kp1 - param$kr) + kp1 * param$kw  * lambda[3,j];
      ks1 = (lambda[2,j] - lambda[1,j]) * u[j];
      dkpdw1 = param$k * kp1 / (param$k + x[3,j])^2
      kw1 = - (dkpdw1 * (lambda[2,j] * x[1,j] - param$kw * x[1,j] * lambda[3,j])) - mu_w
      
      kp2 = param$kp * 0.5 * (x[3,j] + x[3,j-1])/ (param$k + 0.5 * (x[3,j] + x[3,j-1]))
      if(x[3,j] <= 0){
        kp2 = 0
      }
      km2 = -(lambda[2,j]-h2*ks1) * (kp2 - param$kr) - kp2 * param$kw * (lambda[3,j]-h2*kw1);
      ks2 = ((lambda[2,j]-h2*ks1) - (lambda[1,j]-h2*km1)) * 0.5 * (u[j]+u[j-1]);
      dkpdw2 = param$k * kp2 / (param$k + 0.5 * (x[3,j] + x[3,j-1]))^2
      kw2 = - (dkpdw2 * ((lambda[2,j]-h2*ks1)  * 0.5* (x[1,j]+x[1,j-1]) - param$kw * 0.5 * (x[1,j] + x[1,j-1]) * (lambda[3,j] - h2 * kw1))) - mu_w
      
      kp3 = param$kp * 0.5 * (x[3,j] + x[3,j-1])/ (param$k + 0.5 * (x[3,j] + x[3,j-1]))
      if(x[3,j] <= 0){
        kp3 = 0
      }
      km3 = -(lambda[2,j]-h2*ks2) * (kp3 - param$kr) - kp3 * param$kw * (lambda[3,j]-h2*kw2);
      ks3 = ((lambda[2,j]-h2*ks2) - (lambda[1,j]-h2*km2)) * 0.5 * (u[j]+u[j-1]);
      dkpdw3 = param$k * kp3 / (param$k + 0.5 * (x[3,j] + x[3,j-1]))^2
      kw3 = - (dkpdw3 * ((lambda[2,j]-h2*ks2)  * 0.5* (x[1,j]+x[1,j-1]) - param$kw * 0.5 * (x[1,j] + x[1,j-1]) * (lambda[3,j] - h2 * kw2))) - mu_w
      
      kp4 = param$kp * x[3,j-1]/ (param$k + x[3,j-1])
      if(x[3,j-1] <= 0){
        kp4 = 0
      }
      km4 = -(lambda[2,j]-h*ks3) * (kp4 - param$kr) - kp4 * param$kw * (lambda[3,j]-h*kw3);
      ks4 = ((lambda[2,j]-h*ks3) - (lambda[1,j]-h*km3)) * u[j-1];
      dkpdw3 = param$k * kp4 / (param$k + x[3,j-1])^2
      kw4 = - (dkpdw4 * ((lambda[2,j]-h*ks3)  * x[1,j-1] - param$kw * x[1,j-1] * (lambda[3,j] - h * kw3))) - mu_w
      
      lambda[1,j-1] = lambda[1,j] - (h/6)*(km1 + 2*km2 + 2*km3 + km4);
      lambda[2,j-1] = lambda[2,j] - (h/6)*(ks1 + 2*ks2 + 2*ks3 + ks4);
      lambda[3,j-1] = lambda[3,j] - (h/6)*(kw1 + 2*kw2 + 2*kw3 + kw4);
    }
    u1 = u
    for(i in 1:N+1){
      temp = (lambda[1,i] - lambda[2,i]);
      if(temp > 0){
        u1[i] = ubound[1];
      }
      else
        u1[i] = ubound[2];
    }
    
    
    u = 0.5*(u1+oldu);
    
    temp1 = delta*sum(abs(u))-sum(abs(oldu-u));
    temp2 = delta*sum(abs(x))-sum(abs(oldx-x));
    temp3 = delta*sum(abs(lambda))-sum(abs(oldlambda-lambda));
    test = min(temp1, min(temp2, temp3));
    if(iter == maxiter){
      break;
    }
  }
  
  return(list(time = t, x = x, lambda = lambda, u = u));
}