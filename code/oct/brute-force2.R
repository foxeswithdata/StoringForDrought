# nonlinear optimisation package

library('NlcOptim')

source("code/parameters.R")

#Define optimisation inputs
X0 = matrix(rep(0, times = 4*param_1$T_end), nrow=4, ncol=param_1$T_end)
# X0[4,] = rep(param_1$W0, param_1$T_end) # <- might not be necessary

lb = matrix(rep(0, times = 4*param_1$T_end), nrow=4, ncol=param_1$T_end)
ub = matrix(rep(1, times = 4*param_1$T_end), nrow=4, ncol=param_1$T_end)
ub[1,] = ub[1,]*param_1$ks
ub[2:4,] = ub[2:4,]*999999

lb = c(rep(0, times=4*param_1$T_end))
ub = c(rep(param_1$ks, times=param_1$T_end), rep(999999, times=3*param_1$T_end))





# Optimisation functions
nonlcon <- function(X, param, delta=1){ #returns c and ceq
  
  ceq = matrix(rep(0, times=3*(param$T_end+1)),nrow = 3, ncol = param$T_end + 1)
  
  #be careful of the subscripts:
  #X = [ u(0), u(1), ... , u(N-1)
         #      xM(1), xM(2), ..., xM(N)
         #      xS(1), xS(2), ..., xS(N)
         #      xW(1), xW(2), ..., xW(N)]
  #kp = [kp(0), kp(1), ... kp(N-1)]
  for (i in 2:ncol(X)){
    xm = X[2,i-1] + X[1,i] * X[3,i-1];
    kp = 0
    if(X[2,i-1] > 0){
      kp = param$kp
    }
    xs = X[3,i-1] + delta * ((kp - param$kr) * X[2,i-1] - X[1,i] * X[3,i-1])
    xw = X[4,i-1] - param$kw * kp * X[2,i-1]
    ceq[,i] = X[2:4,i] - t(c(xm,xs,xw))
  }
  
  
  #Find the values for initial states
  xm = param$M0 + X[1,1] * param$S0;
  kp = 0
  if(X[2,1] > 0){
    kp = param$kp
  }
  xs = param$S0 + delta * ((kp - param$kr) * param$M0 - X[1,1] * param$S0);
  xw = param$W0 - param$kw * kp * param$M0
  ceq[,1] = X[2:4,1] - t(c(xm,xs,xw))
  
  

  return(list(ceq=ceq,c=NULL))
}

#The objective function
objective <- function(X, param){
  eval = X[3, ncol(X)]
  return(eval)
}


obj <- function(x){
  return(objective(x, param=param_1))
}

nlc <- function(x){
  return(nonlcon(x, param_1))
}


res = solnl(X = X0, objfun = obj, confun = nlc, lb = lb, ub = ub, tolX = 1e-04,
            tolFun = 1e-06, tolCon = 1e-06, maxnFun = 1e+07, maxIter = 3000)

