library('rSymPy')

ts <- Var('ts')
tstar <- Var('tstar')
tstar2 <- Var('tstar2')
kp <- Var('kp')
kw <- Var('kw')
nu <- Var('nu')
kr <- Var('kr')
t <- Var('t')

tstar2 <- sympy("1/(kp * kw * nu)")
tstar3 <- sympy("tstar2 / (kp * (1 - kw * nu) - kr)")
sympy("tstar = tstar / (kp * (1 - kw * nu) - kr)")



sympy("simplify((kp * (1- kw * nu) - kr + kr * t * ( kp * (1 - kw * nu) - kr) + (kr * t + 1) * (kp * (kw * nu - 1) + kr)) / (kp * kw * nu * (kp * (1 - kw * nu) - kr) - (kp * kw * nu - 1) * (kp * (kw * nu - 1) + kr)))")

sympy("simplify((kr * t + 1 - (kr * t + 1) * kp * kw * nu)/(kp * kw * nu * ( kp * (1 - kw * nu) - kr ) + kp * (kw * nu - 1) + kr))")

sympy("simplify((1 - ((1 + kr * t)/(kr - kp + kp * kw * nu)) * (kp * (kw * nu - 1) + kr) + kr * t) / (kp * kw * nu))")


kp_s <- Var("kp_s")
x0 <- Var("x0")
W <- Var("W")
k <- Var("k")
kp <- Var("kp")

sympy("diff(kp_s / (1 + exp(-k * (W - x0))), W, 1)")

W = seq(from = 0, to = 6, by = 0.05)
kp = 0.04 / (1 + exp(-10 * (W - 2)))
kp_diff = (10 * 0.04 * exp(-10 * (W - 2)))/ (1 + exp(-10 * (W - 2)))^2

plot(W, kp, type="l", col="red")
lines(W, kp_diff, col="blue")



##########################
## integrating W(t)_log ##
########################## 
kp <- Var("kp")
M0 <- Var("M0")
kw <- Var("kw")
klog <- Var("klog")
W50 <- Var("W50")
y <- Var("y")
t <- Var("t")
y <- sympy("kp * M0 * kw / (1 + exp(-klog * (y (t) - W50)))")


sympy("integrate(y,t)")

