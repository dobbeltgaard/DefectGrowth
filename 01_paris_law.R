

rm(list = ls())

d = read.csv("defect_trajectories.csv")

#DEFINE VARIABLE FOR EXTRACTION
essential_cols = c("ID","dSize", "t", "UIC")
reference_cols = c("Size1", "Size2","Date1", "Date2", "Rail_string")
maintenan_cols = c("Maintenance_date","Maintenance_indicator", "Grinding", "Milling", "Planing", "Passages","Removed_status")
structure_cols = c("In_straight_track", "In_curve", "In_trans_curve","Turnout_indicator","Weld","Aluminothermic_weld","Flash_butt_weld","Other_weld")
covariate_cols = c("Age", "Curve", "Rail_weight_1m", "Steel_hardness", "Line_speed", "MGT_max")


summary(lm(Size2~Size1 + Maintenance_indicator + load_cycle + Line_speed + Age, d))



table(d$UIC)



library(deSolve)

#THE LIKELIHOOD IS VERY SLOW AND UNSTABLE
#EITHER IMPLEMENT EULER SCHEME IN C++
#OR FIND ALGEBRAIC EXPRESSION

loglik_ode = function(y1, y2, t, pam){
  
  # y1 = d$Size1[idx]
  # y2 = d$Size2[idx]
  # t = d$load_cycle[idx]
  # pam = c(0.1,1, 1)
  
  n = length(y1)
  pars = c(a = pam[1], m = pam[2], X_max = 140)
  sigma = pam[3]
  
  Paris <- function(Time, State, Pars) {
    with(as.list(c(State, Pars)), {
      X <- min(X, X_max)  # Enforce upper bound
      dX <- ifelse(X >= X_max, 0, a * X^m)  # Growth stops, but X stays at X_max
      return(list(c(dX)))
    })
  }
  
  loglik = 0;
  for(i in 1:n){
    times = seq(0, t[i], by = t[i]/200)
    yini <- c(X = y1[i])
    out <- ode(yini, times, Paris, pars)
    xt = tail(out,1)[,2]
    
    loglik = loglik + -1/2*log(sigma^2) - 1/(2*sigma^2)*(xt - y2[i])^2
  }
  return(loglik)
}

idx = 1000:2000
loglik_ode(y1 = d$Size1[idx], y2 = d$Size2[idx], t = d$load_cycle[idx], pam = c(0.1,1, 1))







Paris <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    #X <- min(X, X_max)  # Enforce upper bound
    #dX <- ifelse(X >= X_max, 0, a * X^m)  # Growth stops, but X stays at X_max
    dX <-  a * X^m  # Growth stops, but X stays at X_max
    return(list(c(dX)))
  })
}

pars = c(a = 1, m = 2, X_max = 140)
times = seq(0, 1, by = 0.1)
yini <- c(X = 0.1)
out <- ode(yini, times, Paris, pars)
tail(out)
plot(out[,1], out[,2], type = "p")




Paris_exact = function(t, t0, x0, Pars, xmax = 140){
  a = Pars[1]
  m = Pars[2]
  if(abs(m-1)<1E-10){ 
    sol = x0*exp(a*(t - t0))
  } else {
    sol = (a*m*(t0-t) - a*(t0-t) + x0^(1-m))^(-1/(m-1))}
  if(!is.finite(sol)){ return(Inf) }
  if(sol > xmax){sol = xmax}
  return( sol )
}

Paris_exact(1,0,0.1, c(1.033895,2.810418))



library(nleqslv)

equations <- function(x, params) {
  t_hat <- x[1]
  c <- x[2]
  
  a <- params[1]
  m <- params[2]
  y_1 <- params[3]
  y_2 <- params[4]
  t_delta <- params[5]
  
  eq1 <- (-a * t_hat*m + a * t_hat + c)^(-1 / (m - 1)) - y_1
  eq2 <- (-a * (t_hat+t_delta)*m - a * (t_hat + t_delta) + c)^(-1 / (m - 1)) - y_2
  
  return(c(eq1, eq2))
}

params = c(1,2,12,22,3)

initial_guess <- c(0, 1)

solution <- nleqslv(initial_guess, equations, params = params)

nleqslv()



library("nleqslv")
Paris_pars = function(){
  
  
}

loglik_Paris = function(y1, y2, t, pam){
  n = length(y1)
  pam = exp(pam)
  sigma = pam[3]
  loglik = 0;
  for(i in 1:n){
    i = 10939
    pam = c(1,2); y1 = d$Size1[idx]; y2 = d$Size2[idx]; t = d$t[idx]
    
    a = 1; m=2; t0=0; x0 = y1[i]; t=t[i]
    (a*m*(t0-t) - a*(t0-t) + x0^(1-m))^(-1/(m-1))
    
    
    x = Paris_exact(t[i], 0, y1[i], pam[1:2])
    #if(is.infinite(x)){print(i)}
    loglik = loglik - (x - y2[i])^2 / (2 * sigma^2) 
  }
  loglik = loglik/(2*sigma^2) - n/2 * log(2 * pi * sigma^2) 
  return(-loglik)
}

idx = d$dSize>0
fit = optim(par = c(0.1,1, 1), fn = loglik_Paris, y1 = d$Size1[idx], y2 = d$Size2[idx], t = d$t[idx])
fit
exp(fit$par)


loglik_Paris(y1 = d$Size1[idx], y2 = d$Size2[idx], t = d$t[idx], pam = c(1,2, 1))



Paris_exact(1,0,1, c(1,2))
Paris_exact(100,0,1, c(1,2))

loglik_Paris(y1 = d$Size1[idx], y2 = d$Size2[idx], t = d$t[idx], pam = c(1,2, 1))


plot(d$t[idx], d$dSize[idx])





