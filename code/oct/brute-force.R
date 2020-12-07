#nonlinear run

install.packages("nloptr")

library('nloptr')
library('NlcOptim')

eval_f0 <- function(x){
  
  S0 = 10;
  M0 = 100;
  W0 = 200;
  x0 = c(M0 ,S0 ,W0);
  T_end = 64;
  kp = 0.04;
  kmax = 0.1;
  kr = 0.02;
  delta = 1;
  ubound = c(0, kmax);
  kpbound = c(0, kp);
  kw = 1;
  
  N = T_end;
  
  u <- x
  x <- matrix(nrow = 3, ncol = N+1)
  x[,1] = x0;
  
  
  for(i in 1:N){
    # CALCULATE XM (x(1,:)) XS (x(2,:)) and XW (x(3,:))  
    kp_i = kp;
    if(x[3,i] <= 0){
      kp_i = 0;
    }
    x[1,i+1] = x[1,i] + u[i] * x[2,i];
    x[2,i+1] = x[2,i] + (kp_i - kr) * x[1,i] - u[i] * x[2,i];
    x[3,i+1] = x[3,i] - kw * kp_i * x[1,i];
  }
  
  return(-x[1,N+1])
  
}


res1 <- nloptr( x0 = rep(0.05, times=65),
                eval_f = eval_f0,
                lb  = rep(0, times=65),
                ub = rep(0.1, times=65),
                opts = list("algorithm"="NLOPT_LN_COBYLA",
                            "xtol_rel"=1.0e-8))
res1

plot(1:65, res1$solution)

