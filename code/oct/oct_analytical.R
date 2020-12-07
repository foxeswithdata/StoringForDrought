rm(list=ls())

library(ggplot2)

source('code/oct/forward_backward.R')

S0 = 10;
M0 = 100;
W0 = 250;
X0 = c(M0 ,S0 ,W0);
T_end = 64;
kp = 0.04;
kmax = 0.1;
kr = 0.02;
delta = 1;
ubound = c(0, kmax);
kpbound = c(0, kp);
kw = 1;

# VERSION 1

output1 = forwardback1(kr,kp,kw,X0,T_end,ubound);
output1RK = forwardback1RK(kr,kp,kw,X0,T_end,ubound, nu=2);

plot(output1RK$time, output1RK$x[1,], type='l', col="red", ylim=c(0,W0))
lines(output1RK$time, output1RK$x[2,], col="green")
lines(output1RK$time, output1RK$x[3,], col="blue")



plot(output1RK$time, output1RK$lambda[1,], type='l', col="red", ylim=c(min(output1RK$lambda[1,], min(output1RK$lambda[2,],output1RK$lambda[3,])),max(output1RK$lambda[1,], max(output1RK$lambda[2,],output1RK$lambda[3,]))))
lines(output1RK$time, output1RK$lambda[2,], col="green")
lines(output1RK$time, output1RK$lambda[3,], col="blue")

plot(output1RK$time, output1RK$u, type='l')

output1RK$iter

# NOW LET's WRAP IT to FIND THE CORRECT nu

flag = -1;

a = -1;
b = 5;

z = forwardback1RK(kr,kp,kw,X0,T_end,ubound, nu=a);
Va = z$x[3, T_end+1] - 0; # final W should be 0
z = forwardback1RK(kr,kp,kw,X0,T_end,ubound, nu=b);
Vb = z$x[3, T_end+1] - 0;

iter = 0;

res <- matrix(, ncol=1000, nrow=7);
rownames(res) <- c('a','b','Va','Vb','WT','MT','ST');
res[,iter+1] <- c(a,b,Va,Vb, z$x[3,T_end+1], z$x[1,T_end+1], z$x[2, T_end+1]);



while(flag < 0 && iter < 1000){
  iter = iter + 1;
  if(abs(Va) > abs(Vb)){
    k = a;
    a = b;
    b = k;
    k = Va;
    Va = Vb;
    Vb = k;
  }
  
  d = Va * (b-a)/(Vb-Va);
  b = a;
  Vb = Va;
  a = a - d;
  print(iter)
  print(a)
  print(d)
  z = forwardback1RK(kr,kp,kw,X0,T_end,ubound, nu=a);
  Va = z$x[3, T_end+1] - 0;
  res[,iter+1] <- c(a,b,Va,Vb,z$x[3,T_end+1], z$x[2,T_end+1], z$x[1, T_end+1]);
    
  if(abs(Va) < 10^(-8)){
    flag = 1;
  }
}


### SOO that does not work very well. Let's examine how nu impacts firstly WT and ST

nu = c(-50:-1,2:50);
WTlist = rep(0, times=length(nu))
STlist = rep(0, times=length(nu))

iter = 0
for(i in nu){
  print(i)
  iter = iter+1;
  z = forwardback1RK(kr,kp,kw,X0,T_end,ubound, nu=i);
  WTlist[iter] = z$x[3,T_end+1]
  STlist[iter] = z$x[2,T_end+1]
}

plot(nu, STlist, type="l", col="red")
lines(nu, WTlist, type="l", col="blue")




# VERSION 2





output2 = forwardback2(kr,kw,X0,T_end,ubound, kpbound);
output2RK = forwardback2RK(kr,kw,X0,T_end,ubound, kpbound, nu=-3)
output3DT =forwardback3DT(kr,kp,kw,X0,T_end,ubound);
output3RK =forwardback3RK(kr,kp,kw,X0,T_end,ubound);

output2$u == output2RK$u


plot(output3RK$x[2,] ~ output3DT$time)
 
kprun = kp_logistic(kp, 100, 0, output3DT$x[3,])
kprun

plot(output3DT$time, output3DT$x[1,], type='l', col="red", ylim=c(0,W0))
lines(output3DT$time, output3DT$x[2,], col="green")
lines(output3DT$time, output3DT$x[3,], col="blue")

plot(output3DT$time, output3DT$lambda[1,], type='l', col="red", ylim=c(min(output3DT$lambda[1,], min(output3DT$lambda[2,],output3DT$lambda[3,])),max(output3DT$lambda[1,], max(output3DT$lambda[2,],output3DT$lambda[3,]))))
lines(output3DT$time, output3DT$lambda[2,], col="green")
lines(output3DT$time, output3DT$lambda[3,], col="blue")



# VERSION 4

w50 = 2
k = 10

output <- forwardback_kp_log_RK(kp, kr, kw, X0, T_end, ubound, k = 10, w50 = 2)

plot(output$time, output$x[1,], type='l', col="red", ylim=c(-40,W0))
lines(output$time, output$x[2,], col="green")
lines(output$time, output$x[3,], col="blue")
