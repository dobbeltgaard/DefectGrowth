

rm(list = ls())

d = read.csv("defect_trajectories.csv")

#DEFINE VARIABLE FOR EXTRACTION
essential_cols = c("ID","dSize", "t", "UIC")
reference_cols = c("Size1", "Size2","Date1", "Date2", "Rail_string")
maintenan_cols = c("Maintenance_date","Maintenance_indicator", "Grinding", "Milling", "Planing", "Passages","Removed_status")
structure_cols = c("In_straight_track", "In_curve", "In_trans_curve","Turnout_indicator","Weld","Aluminothermic_weld","Flash_butt_weld","Other_weld")
covariate_cols = c("Age", "Curve", "Rail_weight_1m", "Steel_hardness", "Line_speed", "MGT_max")

d$load_cycle = d$t*d$MGT_max 
d = d[!(is.na(d$load_cycle) | d$load_cycle == 0), ]




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



