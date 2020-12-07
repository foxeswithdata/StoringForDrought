
forwardback1 <- function(kr,kp,kw,x0,N,ubound, nu = 0){
  iter = 0;
  
  test = -1; #this is to find convergence
  
  delta = 0.001; # convergence value
  t = seq(from=0,to=N,length=N+1); 
  h = 1/N; # spacing between nodes
  h2 = h/2; # used for runge-kutta subroutine

  
  u = rep(0,times=N+1); # this is essentially ks
  
  x <- matrix(0, 3, N+1) # initialise guess for state variables

  x[,1] = x0;
  
  lambda = matrix(1,3,N+1); # initialise guess for the adjoint
  
  #The following are from the analytical
  lambda[1,] = lambda[1,] * 0;
  lambda[2,] = lambda[2,] * 1;
  lambda[3,] = lambda[3,] * nu;
  
  while(test < 0){
    iter = iter+1;
    
    oldu = u;
    oldx = x;
    oldlambda = lambda;
    
    kp_lambda = rep(1,times=N+1);
    
    # The following is runge kutta as an approximation of the
    # continuous time - need to read up more about it
    
    for(i in 1:N){
      # CALCULATE XM (x(1,:)) XS (x(2,:)) and XW (x(3,:))  
      kp_i = kp;
      if(x[3,i] <= 0){
        kp_i = 0;
      }
      
      kp_lambda[i] = kp_i;
      
      x[1,i+1] = x[1,i] + u[i] * x[2,i];
      x[2,i+1] = x[2,i] + (kp_i - kr) * x[1,i] - u[i] * x[2,i];
      x[3,i+1] = x[3,i] - kw * kp_i * x[1,i];
    }
    
    
    # Runge-Kutta sweep backwards to find lambda
    for(i in 1:N){
      j = N + 2 - i;
      
      lambda[1,j-1] = lambda[1,j] + lambda[3,j] + lambda[2,j] * (kp_lambda[j] - kr);
      lambda[2,j-1] = lambda[2,j] - u[j] * (lambda[2,j] - lambda[1,j]);
      lambda[3,j-1] = lambda[3,j];
    }
    u1 = u
    for(i in 1:N+1){
      temp = (lambda[1,i] - lambda[2,i]);
      if(temp < 0){
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
    
    
  }
  return(list(time = t, x = x, lambda = lambda, u = u));
  
}



forwardback1RK <- function(kr,kp,kw,x0,N,ubound, nu = 0){
  iter = 0;
  test = -1; #this is to find convergence
  
  delta = 0.001; # convergence value
  
  t = seq(from=0,to=N,length=N+1); 
  
  h = N/N; # spacing between nodes
  h2 = h/2; # used for runge-kutta subroutine
  
  u = rep(0,times=N+1); # this is essentially ks
  
  x <- matrix(0, 3, N+1) # initialise guess for state variables
  
  x[,1] = x0;
  
  lambda = matrix(1,3,N+1); # initialise guess for the adjoint
  
  # The following are from the analytical
  # These will never be changed so it is CRUCIAL to get these numbers right
  lambda[1,N+1] = 0;
  lambda[2,N+1] = 1;
  lambda[3,N+1] = nu;
  
  kp_lambda = rep(1,times=N+1);
  
  while(test < 0){
    iter = iter+1;
    
    oldu = u;
    oldx = x;
    oldlambda = lambda;
    
    # The following is runge kutta as an approximation of the
    # continuous time - need to read up more about it
    
    for(i in 1:N){
      # CALCULATE XM (x(1,:)) XS (x(2,:)) and XW (x(3,:))  
      kp_i = kp;
      if(x[3,i] <= 0){
        kp_i = 0;
      }
      kp_lambda[i] = kp_i;
      
      k11 = x[2,i] * u[i];
      k12 = (kp_i - kr) * x[1,i] - x[2,i] * u[i];
      k13 = - kw * kp_i * x[1,i];
      
      k21 = (x[2,i] + h2*k12) * 0.5 * (u[i] + u[i+1]);
      k22 = (kp_i - kr) * (x[1,i] + h2*k11) - (x[2,i] + h2*k12) * 0.5*(u[i] + u[i+1]);
      k23 = -kw * kp_i * (x[1,i] + h2 * k11);
      
      k31 = (x[2,i] + h2*k22) * 0.5*(u[i] + u[i+1]);
      k32 = (kp_i - kr) * (x[1,i] + h2*k21) - (x[2,i] + h2*k22) * 0.5*(u[i] + u[i+1]);
      k33 = -kw * kp_i * (x[1,i] + h2*k21);
      
      k41 = (x[2,i] + h*k32) * u[i+1];
      k42 = (kp_i - kr) * (x[1,i] + h*k31) - (x[2,i] + h*k32) * u[i+1];
      k43 = -kw * kp_i * (x[1,i] + h * k31);
      
      x[1,i+1] = x[1,i] + (h/6)*(k11+2*k21+2*k31+k41);
      x[2,i+1] = x[2,i] + (h/6)*(k12+2*k22+2*k32+k42);
      x[3,i+1] = x[3,i] + (h/6)*(k13+2*k23+2*k33+k43);
      
      
    }
    
    
    
    # Runge-Kutta sweep backwards to find lambda
    for(i in 1:N){
      j = N + 2 - i;
      k11 = -lambda[2,j] * (kp_lambda[j] - kr) - kp_lambda[j] * kw * lambda[3,j];
      k21 = (lambda[2,j] - lambda[1,j]) * u[j];
      k31 = 0;
      
      k12 = -(lambda[2,j]-h2 * k21) * (0.5* (kp_lambda[j]+kp_lambda[j-1]) - kr) - 0.5* (kp_lambda[j]+kp_lambda[j-1]) * kw * (lambda[3,j]-h2 * k31);
      k22 = ((lambda[2,j]-h2*k21) - (lambda[1,j]-h2*k11)) * 0.5 * (u[j]+u[j-1]);
      k32 = 0;
      
      k13 = -(lambda[2,j]-h2 * k22) * (0.5* (kp_lambda[j]+kp_lambda[j-1]) - kr) - 0.5* (kp_lambda[j]+kp_lambda[j-1]) * kw * (lambda[3,j] - h2 * k32);
      k23 = ((lambda[2,j]-h2*k22) - (lambda[1,j]-h2*k12)) * 0.5 * (u[j]+u[j-1]);
      k33 = 0;
      
      k14 = -(lambda[2,j]-h * k23) * (kp_lambda[j-1] - kr) - kp_lambda[j-1] * kw * (lambda[3,j] - h * k33);
      k24 = ((lambda[2,j]-h*k23) - (lambda[1,j]-h*k13)) * u[j-1];
      k34 = 0;
      
      lambda[1,j-1] = lambda[1,j] - (h/6)*(k11 + 2*k12 + 2*k13 + k14);
      lambda[2,j-1] = lambda[2,j] - (h/6)*(k21 + 2*k22 + 2*k23 + k24);
      lambda[3,j-1] = lambda[3,j] - (h/6)*(k31 + 2*k32 + 2*k33 + k34);
    }
    
    
    # df <- data.frame(time=t, m = x[1,], s = x[2,], w = x[3,], u = u, kp = kp_lambda)
    # 
    # p <- ggplot(df, aes(x=t)) +
    #   geom_line(aes(y=m, colour="MT")) +
    #   geom_line(aes(y=s, colour="ST")) +
    #   geom_line(aes(y=w, colour="WT")) +
    #   xlab("time") + 
    #   ylab("pools")+ 
    #   scale_colour_manual("", 
    #                       values = c("MT"="#32A287", "ST"="#6C464E", "WT"="#A63A50"))
    # 
    # ggsave(paste0("xiter", iter, ".png"), plot = p, device = png(), path = "figures/for_back_test", dpi = 400)
    # 
    # p <- ggplot(df, aes(x=t)) +
    #   geom_line(aes(y=u)) +
    #   xlab("time") + 
    #   ylab("u")
    # 
    # ggsave(paste0("uiter", iter, ".png"), plot = p, device = png(), path = "figures/for_back_test", dpi = 400)
    
    u1 = u
    for(i in 1:N+1){
      temp = (lambda[1,i] - lambda[2,i]);
      if(temp < 0){
        u1[i] = ubound[1];
      }
      else
        u1[i] = ubound[2];
    }
    
    u = 0.5*(u1+oldu);
    
    # graphics.off()
    
    temp1 = delta*sum(abs(u))-sum(abs(oldu-u));
    temp2 = delta*sum(abs(x))-sum(abs(oldx-x));
    temp3 = delta*sum(abs(lambda))-sum(abs(oldlambda-lambda));
    test = min(temp1, min(temp2, temp3));
    
  }
  return(list(time = t, x = x, lambda = lambda, u = u, iter=iter, kp = kp_lambda));
  
}



forwardback2 <- function(kr,kw,x0,N,ubound, kpbound, nu = 0){
  iter = 0;
  test = -1; #this is to find convergence
  
  delta = 0.001; # convergence value
  t = seq(from=0,to=N,length=N+1); 
  h = 1/N; # spacing between nodes
  h2 = h/2; # used for runge-kutta subroutine
  
  
  u = rep(0,times=N+1); # this is essentially ks
  kp = rep(0,times=N+1); # this is essentially ks
  
  x <- matrix(0, 3, N+1) # initialise guess for state variables
  
  x[,1] = x0;
  
  lambda = matrix(1,3,N+1); # initialise guess for the adjoint
  
  #The following are from the analytical
  lambda[1,] = lambda[1,] * 0;
  lambda[2,] = lambda[2,] * 1;
  lambda[3,] = lambda[3,] * nu;
  
  while(test < 0){
    iter = iter+1;
    
    oldu = u;
    oldkp = kp;
    oldx = x;
    oldlambda = lambda;
    
    
    # The following is runge kutta as an approximation of the
    # continuous time - need to read up more about it
    
    for(i in 1:N){
      # CALCULATE XM (x(1,:)) XS (x(2,:)) and XW (x(3,:))  

      x[1,i+1] = x[1,i] + u[i] * x[2,i];
      x[2,i+1] = x[2,i] + (kp[i] - kr) * x[1,i] - u[i] * x[2,i];
      x[3,i+1] = x[3,i] - kw * kp[i] * x[1,i];
      
    }
    
    
    # Runge-Kutta sweep backwards to find lambda
    for(i in 1:N){
      j = N + 2 - i;
      
      lambda[1,j-1] = lambda[1,j] + lambda[3,j] + lambda[2,j] * (kp[j] - kr);
      lambda[2,j-1] = lambda[2,j] - u[j] * (lambda[2,j] - lambda[1,j]);
      lambda[3,j-1] = lambda[3,j];
    }
    
    u1 = u
    for(i in 1:N+1){
      temp = (lambda[1,i] - lambda[2,i]);
      if(temp < 0){
        u1[i] = ubound[1];
      }
      else
        u1[i] = ubound[2];
    }
    
    kp1 = kp;
    for(i in 1:N+1){
      temp = (lambda[3,i]*kw - lambda[2,i]);
      if(temp < 0){
        kp1[i] = ubound[2];
      }
      else
        kp1[i] = ubound[1];
    }
    
    
    u = 0.5*(u1+oldu);
    kp = 0.5*(kp1+oldkp);
    
    temp1 = delta*sum(abs(u))-sum(abs(oldu-u));
    temp2 = delta*sum(abs(kp))-sum(abs(oldkp-kp))
    temp3 = delta*sum(abs(x))-sum(abs(oldx-x));
    temp4 = delta*sum(abs(lambda))-sum(abs(oldlambda-lambda));
    test = min(temp1, min(temp2, min(temp3, temp4)));
    
    
  }
  
  
  
  return(list(time = t, x = x, lambda = lambda, u = u));
  
}



forwardback2RK <- function(kr,kw,x0,N,ubound, kpbound, nu = 0){
  iter = 0;
  
  test = -1; #this is to find convergence
  
  delta = 0.001; # convergence value
  t = seq(from=0,to=N,length=N+1); 
  h = N/N; # spacing between nodes
  h2 = h/2; # used for runge-kutta subroutine
  
  
  u = rep(sum(ubound)/2,times=N+1); # this is essentially ks
  kp = rep(sum(kpbound)/2,times=N+1); # this is essentially ks
  
  x <- matrix(0, 3, N+1) # initialise guess for state variables
  
  x[,1] = x0;
  
  lambda = matrix(1,3,N+1); # initialise guess for the adjoint
  
  #The following are from the analytical
  lambda[1,] = lambda[1,] * 0;
  lambda[2,] = lambda[2,] * 1;
  lambda[3,] = lambda[3,] * nu;
  
  while(test < 0){
    iter = iter+1;
    
    oldu = u;
    oldkp = kp;
    oldx = x;
    oldlambda = lambda;
    
    
    # The following is runge kutta as an approximation of the
    # continuous time - need to read up more about it
    
    for(i in 1:N){
      # CALCULATE XM (x(1,:)) XS (x(2,:)) and XW (x(3,:))  
      
      k11 = x[2,i] * u[i];
      k21 = (kp[i] - kr) * x[1,i] - x[2,i] * u[i];
      k31 = -kw * kp[i] * x[1,i];
      
      k12 = (x[2,i] + h2*k21) * 0.5*(u[i] + u[i+1]);
      k22 = (0.5 * (kp[i] + kp[i+1]) - kr) * (x[1,i] + h2*k11) - (x[2,i] + h2*k21) * 0.5*(u[i] + u[i+1]);
      k32 = -kw * 0.5 * (kp[i] + kp[i+1]) * (x[1,i] + h2 * k11);
      
      k13 = (x[2,i] + h2*k22) * 0.5*(u[i] + u[i+1]);
      k23 = (0.5 * (kp[i] + kp[i+1]) - kr) * (x[1,i] + h2*k12) - (x[2,i] + h2*k22) * 0.5*(u[i] + u[i+1]);
      k33 = -kw * 0.5 * (kp[i] + kp[i+1]) * (x[1,i] + h2*k12);
      
      k14 = (x[2,i] + h*k23) * u[i+1];
      k24 = (kp[i+1] - kr) * (x[1,i] + h*k13) - (x[2,i] + h*k23) * u[i+1];
      k34 = -kw * kp[i+1] * (x[1,i] + h * k13);
      
      x[1,i+1] = x[1,i] + (h/6)*(k11+2*k12+2*k13+k14);
      x[2,i+1] = x[2,i] + (h/6)*(k21+2*k22+2*k23+k24);
      x[3,i+1] = x[3,i] + (h/6)*(k31+2*k32+2*k33+k34);
      
    }
    
    
    # Runge-Kutta sweep backwards to find lambda
    for(i in 1:N){
      j = N + 2 - i;
      
      k11 = -lambda[2,j] * (kp[j] - kr) - kp[j] * kw  * lambda[3,j];
      k21 = (lambda[2,j] - lambda[1,j]) * u[j];
      k31 = 0
      
      k12 = -(lambda[2,j]-h2*k21) * (0.5*(kp[j]+kp[j-1]) - kr) - 0.5* (kp[j]+kp[j-1]) * kw * (lambda[3,j]-h2*k31);
      k22 = ((lambda[2,j]-h2*k21) - (lambda[1,j]-h2*k11)) * 0.5 * (u[j]+u[j-1]);
      k32 = 0
      
      k13 = -(lambda[2,j]-h2*k22) * (0.5*(kp[j]+kp[j-1]) - kr) - 0.5* (kp[j]+kp[j-1]) * kw * (lambda[3,j]-h2*k32);
      k23 = ((lambda[2,j]-h2*k22) - (lambda[1,j]-h2*k12)) * 0.5 * (u[j]+u[j-1]);
      k33 = 0
      
      k14 = -(lambda[2,j]-h * k23) * (kp[j-1] - kr) - kp[j-1] * kw * (lambda[3,j] - h * k33);
      k24 = ((lambda[2,j]-h*k23) - (lambda[1,j]-h*k13)) * u[j-1];
      k34 = 0
      
      lambda[1,j-1] = lambda[1,j] - (h/6)*(k11 + 2*k12 + 2*k13 + k14);
      lambda[2,j-1] = lambda[2,j] - (h/6)*(k21 + 2*k22 + 2*k23 + k24);
      lambda[3,j-1] = lambda[3,j] - (h/6)*(k21 + 2*k22 + 2*k23 + k24);
    }
    u1 = u
    for(i in 1:N+1){
      temp = (lambda[1,i] - lambda[2,i]);
      if(temp < 0){
        u1[i] = ubound[1];
      }
      else
        u1[i] = ubound[2];
    }
    
    kp1 = kp;
    for(i in 1:N+1){
      temp = (lambda[2,i] - lambda[3,i]*kw);
      if(temp < 0){
        kp1[i] = ubound[1];
      }
      else
        kp1[i] = ubound[2];
    }
    
    
    u = 0.5*(u1+oldu);
    kp = 0.5*(kp1+oldkp);
    
    temp1 = delta*sum(abs(u))-sum(abs(oldu-u));
    temp2 = delta*sum(abs(kp))-sum(abs(oldkp-kp))
    temp3 = delta*sum(abs(x))-sum(abs(oldx-x));
    temp4 = delta*sum(abs(lambda))-sum(abs(oldlambda-lambda));
    test = min(temp1, min(temp2, min(temp3, temp4)));
    
    
  }
  
  
  return(list(time = t, x = x, kp = kp, lambda = lambda, u = u));
  
}


forwardback3DT <- function(kr,kp,kw,x0,N,ubound, nu = 0){
  iter = 0;
  kpmax = kp;
  test = -1; #this is to find convergence
  
  delta = 0.001; # convergence value
  t = seq(from=0,to=N,length=N+1); 
  h = 1/N; # spacing between nodes
  h2 = h/2; # used for runge-kutta subroutine
  
  
  u = rep(0,times=N+1); # this is essentially ks
  
  x <- matrix(0, 3, N+1) # initialise guess for state variables
  
  x[,1] = x0;
  
  lambda = matrix(1,3,N+1); # initialise guess for the adjoint
  
  #The following are from the analytical
  lambda[1,] = lambda[1,] * 0;
  lambda[2,] = lambda[2,] * 1;
  lambda[3,] = lambda[3,] * nu;
  
  while(test < 0 & iter < 1000){
    print(iter)
    iter = iter+1;
    
    oldu = u;
    oldx = x;
    oldlambda = lambda;
    
    kp_lambda = rep(1,times=N+1);
    
    # The following is runge kutta as an approximation of the
    # continuous time - need to read up more about it
    
    for(i in 1:N){
      # CALCULATE XM (x(1,:)) XS (x(2,:)) and XW (x(3,:))  
      kp_i = kpmax / (1 + exp(-100 * (x[3,i] - 0.1)));
      
      x[1,i+1] = x[1,i] + u[i] * x[2,i];
      x[2,i+1] = x[2,i] + (kp_i - kr) * x[1,i] - u[i] * x[2,i];
      x[3,i+1] = x[3,i] - kw * kp_i * x[1,i];
    }
    
    
    # Runge-Kutta sweep backwards to find lambda
    for(i in 1:N){
      j = N + 2 - i;
      alpha = exp(-100 * (x[3,j] - 0.1))
      kp_i = kpmax / (1 + alpha);
      lambda[1,j-1] = lambda[1,j] + lambda[3,j] + lambda[2,j] * (kp_i - kr);
      lambda[2,j-1] = lambda[2,j] - u[j] * (lambda[2,j] - lambda[1,j]);
      lambda[3,j-1] = lambda[3,j] + (x[1,j] * lambda[2,j] + lambda[3,j] * kw) * (100 * kpmax * alpha / (1 + alpha)^2) 
    }
    u1 = u
    for(i in 1:N+1){
      temp = (lambda[1,i] - lambda[2,i]);
      if(temp < 0){
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
    
    
  }
  return(list(time = t, x = x, lambda = lambda, u = u));
  
}



forwardback3RK <- function(kr,kp,kw,x0,N,ubound, nu = 0){
  iter = 0;
  kpmax = kp;
  test = -1; #this is to find convergence
  
  delta = 0.001; # convergence value
  t = seq(from=0,to=N,length=N+1); 
  h = N/N; # spacing between nodes
  h2 = h/2; # used for runge-kutta subroutine
  
  
  u = rep(0,times=N+1); # this is essentially ks
  
  x <- matrix(0, 3, N+1) # initialise guess for state variables
  
  x[,1] = x0;
  
  lambda = matrix(1,3,N+1); # initialise guess for the adjoint
  
  #The following are from the analytical
  lambda[1,] = lambda[1,] * 0;
  lambda[2,] = lambda[2,] * 1;
  lambda[3,] = lambda[3,] * nu;
  
  while(test < 0){
    iter = iter+1;
    
    oldu = u;
    oldx = x;
    oldlambda = lambda;
    
    kp_lambda = rep(1,times=N+1);
    
    # The following is runge kutta as an approximation of the
    # continuous time - need to read up more about it
    
    for(i in 1:N){
      # CALCULATE XM (x(1,:)) XS (x(2,:)) and XW (x(3,:))  
      
      k11 = x[2,i] * u[i];
      k21 = (kp_logistic(kmax, 100, 0.1, x[3,i]) - kr) * x[1,i] -x[2,i] * u[i];
      k31 = -kw * kp_logistic(kmax, 100, 0.1, x[3,i]) * x[1,i];
      
      k12 = (x[2,i] + h2*k21) * 0.5*(u[i] + u[i+1]);
      k22 = (kp_logistic(kmax, 100, 0.1, x[3,i] + h2 * k31) - kr) * (x[1,i] + h2*k11) - (x[2,i] + h2*k21) * 0.5*(u[i] + u[i+1]);
      k32 = -kw * kp_logistic(kmax, 100, 0.1, x[3,i] + h2 * k31) * (x[1,i] + h2 * k11);
      
      k13 = (x[2,i] + h2*k22) * 0.5*(u[i] + u[i+1]);
      k23 = (kp_logistic(kmax, 100, 0.1, x[3,i] + h2 * k32) - kr) * (x[1,i] + h2*k12) - (x[2,i] + h2*k22) * 0.5*(u[i] + u[i+1]);
      k33 = -kw * kp_logistic(kmax, 100, 0.1, x[3,i] + h2 * k32) * (x[1,i] + h2*k12);
      
      k14 = (x[2,i] + h*k23) * u[i+1];
      k24 = (kp_logistic(kmax, 100, 0.1, x[3,i] + h * k33) - kr) * (x[1,i] + h*k13) - (x[2,i] + h*k23) * u[i+1];
      k34 = -kw * kp_logistic(kmax, 100, 0.1, x[3,i] + h * k33) * (x[1,i] + h * k13);
      
      x[1,i+1] = x[1,i] + (h/6)*(k11+2*k12+2*k13+k14);
      x[2,i+1] = x[2,i] + (h/6)*(k21+2*k22+2*k23+k24);
      x[3,i+1] = x[3,i] + (h/6)*(k31+2*k32+2*k33+k34);
      
    }
    
    
    # Runge-Kutta sweep backwards to find lambda
    for(i in 1:N){
      j = N + 2 - i;
      k11 = -lambda[2,j] * (kp_logistic(kmax, 100, 0.1, x[3,j]) - kr) - kp_logistic(kmax, 100, 0.1, x[3,j]) * kw * lambda[3,j];
      k21 = (lambda[2,j] - lambda[1,j]) * u[j];
      k31 = -(lambda[2,j] * x[1,j] + kw * lambda[3,j]) * (100 * kpmax * exp(-100 * (x[3,j]))) / (exp(-100 * (x[3,j]))+1)^2
      
      k12 = -(lambda[2,j]-h2*k21) * (0.5*(kp_logistic(kmax, 100, 0.1, x[3,j]) + kp_logistic(kmax, 100, 0.1, x[3,j-1])) - kr) - 0.5* (kp_logistic(kmax, 100, 0.1, x[3,j]) + kp_logistic(kmax, 100, 0.1, x[3,j-1])) * kw * (x[3,j] - h2 * k31);
      k22 = ((lambda[2,j]-h2*k21) - (lambda[1,j]-h2*k11)) * 0.5 * (u[j]+u[j-1]);
      k31 = -((lambda[2,j]-h2*k21) * 0.5*(x[1,j] + x[1,j-1]) + kw * (lambda[3,j]-h2*k31)) * (100 * kpmax * exp(-100 * (0.5 * (x[3,j]+x[3,j-1])))) / (exp(-100 * (0.5 * (x[3,j]+x[3,j-1])))+1)^2
      
      
      k13 = -(lambda[2,j]-h2 * k22) * (0.5* (kp_logistic(kmax, 100, 0.1, x[3,j]) + kp_logistic(kmax, 100, 0.1, x[3,j-1])) - kr) - 0.5* (kp_logistic(kmax, 100, 0.1, x[3,j]) + kp_logistic(kmax, 100, 0.1, x[3,j-1])) * kw * (x[3,j] - h2 * k32);
      k23 = ((lambda[2,j]-h2*k22) - (lambda[1,j]-h2*k12)) * 0.5 * (u[j]+u[j-1]);
      k31 = ((lambda[2,j]-h2*k22) * 0.5*(x[1,j] + x[1,j-1]) + kw * (lambda[3,j]-h2*k32)) * (100 * kpmax * exp(-100 * (0.5 * (x[3,j]+x[3,j-1])))) / (exp(-100 * (0.5 * (x[3,j]+x[3,j-1])))+1)^2
      
      k14 = -(lambda[2,j]-h * k23) * (kp_logistic(kmax, 100, 0.1, x[3,j-1]) - kr) - kp_logistic(kmax, 100, 0.1, x[3,j-1]) * kw * (lambda[3,j] - h * k33);
      k24 = ((lambda[2,j]-h*k23) - (lambda[1,j]-h*k13)) * u[j-1];
      k34 = ((lambda[2,j]-h*k23) * x[1,j-1] + kw * (lambda[3,j]-h*k33)) * (100 * kpmax * exp(-100 * x[3,j-1])) / (exp(-100 * x[3,j-1])+1)^2
      
      lambda[1,j-1] = lambda[1,j] - (h/6)*(k11 + 2*k12 + 2*k13 + k14);
      lambda[2,j-1] = lambda[2,j] - (h/6)*(k21 + 2*k22 + 2*k23 + k24);
      lambda[3,j-1] = lambda[3,j] - (h/6)*(k31 + 2*k32 + 2*k33 + k34);
    }
    u1 = u
    for(i in 1:N+1){
      temp = (lambda[1,i] - lambda[2,i]);
      if(temp < 0){
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
    
    
  }
  return(list(time = t, x = x, lambda = lambda, u = u));
  
}

kp_logistic <- function(L, k, x0, x){
  
  return (L /(1 + exp(- k * (x - x0))));
  
}



forwardback_kp_log_RK <- function(kp, kr,kw,x0,N,ubound, k, w50, maxiter = 100000000){
  iter = 0;
  
  test = -1; #this is to find convergence
  
  delta = 0.001; # convergence value
  t = seq(from=0,to=N,length=N+1); 
  h = N/N; # spacing between nodes
  h2 = h/2; # used for runge-kutta subroutine
  
  #logistic values
  # w0 = 2
  # k = 10
  
  u = rep(sum(ubound)/2,times=N+1); # this is essentially ks
  
  x <- matrix(0, 3, N+1) # initialise guess for state variables
  
  x[,1] = x0;
  
  lambda = matrix(1,3,N+1); # initialise guess for the adjoint
  
  #The following are from the analytical
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
      
      kp1 = kp /(1 + exp(- k * (x[3,i] - w50)))
      km1 = x[2,i] * u[i];
      ks1 = (kp1 - kr) * x[1,i] - x[2,i] * u[i];
      kw1 = -kw * kp1 * x[1,i];
      
      kp2 = kp / (1 + exp(-k * ((x[3,i] + h2*kw1) - w50)))
      km2 = (x[2,i] + h2*ks1) * 0.5*(u[i] + u[i+1]);
      ks2 = (kp2 - kr) * (x[1,i] + h2*km1) - (x[2,i] + h2*ks1) * 0.5*(u[i] + u[i+1]);
      kw2 = -kw * kp2 * (x[1,i] + h2 * km1);
      
      kp3 = kp / (1 + exp(-k * ((x[3,i] + h2*kw2) - w50)))
      km3 = (x[2,i] + h2*ks2) * 0.5*(u[i] + u[i+1]);
      ks3 = (kp3 - kr) * (x[1,i] + h2*km2) - (x[2,i] + h2*ks2) * 0.5*(u[i] + u[i+1]);
      kw3 = -kw * kp3 * (x[1,i] + h2*km2);
      
      kp4 = kp / (1 + exp(-k * ((x[3,i] + h * kw3) - w50)))
      km4 = (x[2,i] + h*ks3) * u[i+1];
      ks4 = (kp4 - kr) * (x[1,i] + h*km3) - (x[2,i] + h*ks3) * u[i+1];
      kw4 = -kw * kp4 * (x[1,i] + h * km3);
      
      x[1,i+1] = x[1,i] + (h/6)*(km1+2*km2+2*km3+km4);
      x[2,i+1] = x[2,i] + (h/6)*(ks1+2*ks2+2*ks3+ks4);
      x[3,i+1] = x[3,i] + (h/6)*(kw1+2*kw2+2*kw3+kw4);
      
    }
    
    
    # Runge-Kutta sweep backwards to find lambda
    for(i in 1:N){
      j = N + 2 - i;
      
      kp1 = kp /(1 + exp(- k * (x[3,j] - w50)))
      km1 = -lambda[2,j] * (kp1 - kr) + kp1 * kw  * lambda[3,j];
      ks1 = (lambda[2,j] - lambda[1,j]) * u[j];
      dkpdw1 = (k * kp * exp(-k * (x[3,j] - w50)))/(1 + exp(-k * (x[3,j] - w50)))^2
      kw1 = - (dkpdw1 * (lambda[2,j] * x[1,j] - kw * x[1,j] * lambda[3,j]))
      
      kp2 = kp /(1 + exp(- k * (0.5 * (x[3,j] + x[3,j-1]) - w50)))
      km2 = -(lambda[2,j]-h2*ks1) * (kp2 - kr) - kp2 * kw * (lambda[3,j]-h2*kw1);
      ks2 = ((lambda[2,j]-h2*ks1) - (lambda[1,j]-h2*km1)) * 0.5 * (u[j]+u[j-1]);
      dkpdw2 = (k * kp * exp(-k * (0.5 * (x[3,j] + x[3,j-1]) - w50)))/(1 + exp(-k * (0.5 * (x[3,j] + x[3,j-1]) - w50)))^2
      kw2 = - (dkpdw2 * ((lambda[2,j]-h2*ks1)  * 0.5* (x[1,j]+x[1,j-1]) - kw * 0.5 * (x[1,j] + x[1,j-1]) * (lambda[3,j] - h2 * kw1)))
      
      kp3 = kp /(1 + exp(- k * (0.5 * (x[3,j] + x[3,j-1]) - w50)))
      km3 = -(lambda[2,j]-h2*ks2) * (kp3 - kr) - kp3 * kw * (lambda[3,j]-h2*kw2);
      ks3 = ((lambda[2,j]-h2*ks2) - (lambda[1,j]-h2*km2)) * 0.5 * (u[j]+u[j-1]);
      dkpdw3 = (k * kp * exp(-k * (0.5 * (x[3,j] + x[3,j-1]) - w50)))/(1 + exp(-k * (0.5 * (x[3,j] + x[3,j-1]) - w50)))^2
      kw3 = - (dkpdw3 * ((lambda[2,j]-h2*ks2)  * 0.5* (x[1,j]+x[1,j-1]) - kw * 0.5 * (x[1,j] + x[1,j-1]) * (lambda[3,j] - h2 * kw2)))
      
      kp4 = kp /(1 + exp(- k * (x[3,j-1] - w50)))
      km4 = -(lambda[2,j]-h*ks3) * (kp4 - kr) - kp4 * kw * (lambda[3,j]-h*kw3);
      ks4 = ((lambda[2,j]-h*ks3) - (lambda[1,j]-h*km3)) * u[j-1];
      dkpdw4 = (k * kp * exp(-k * (x[3,j-1] - w50)))/(1 + exp(-k * (x[3,j-1] - w50)))^2
      kw4 = - (dkpdw4 * ((lambda[2,j]-h*ks3)  * x[1,j-1] - kw * x[1,j-1] * (lambda[3,j] - h * kw3)))
      
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




