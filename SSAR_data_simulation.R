# data simulation script

# This script simulates data on nesting beach counts of leatherbacks. Simulated data are analyzed
# using various models to see how imputations and Pareto - k statistics perform.

# parameters

TT <- 16  # years of time series data
n.states <- 2   # # nesting beaches
n.months <- 12
periods <- c(12, 6)

N0.mean <- c(8.1, 7.9)   # initial values for the abundance mean
N0_sd <- c(1.4, 1.2)       # SD of the initial N 

U.mean <- c(-0.12, -0.10)    # annual growth rates for the two population
U.sd <- c(0.15, 0.15)

N0 <- numeric(length = n.states)
N <- predN <- array(dim = c(n.states, TT))

for (k in 1:n.states){
  N0[k] <- rnorm(n = 1, mean = N0.mean[k], sd = N0_sd[k])
  predN[k, 1] <- N0[k] + rnorm(n = 1, mean = U.mean[k], sd = U.sd[k])
  N[k, 1] <- predN[k, 1]
}

C0_Q <- C0 <- c(0.46, 2.91)

p.beta.cos <- c(-1.27, -4.08)
p.beta.sin <- c(4.18, -2.81)

sigma.R <- c(0.96, 0.61)
df <- c(49.9, 39.8)

p.const <- p <- sigma.Q <- array(dim = c(n.states, n.months))
for (j in 1:n.states){
  for (k in 1:n.months){
    p.const[j, k] <-  2 * pi * k / period[j]
    p[j, k] <- (C0[j] + p.beta.cos[j] * cos(p.const[j,k]) + p.beta.sin[j] * sin(p.const[j,k]))/(n.months * C0[j] + p.beta.cos[j] * C_cos[j] + p.beta.sin[j] * C_sin[j])
    
    sigma.Q[j,k] <- (C0_Q[j] + p.beta.cos[j] * cos(p.const[j,k]) + p.beta.sin[j] * sin(p.const[j,k]))/(n.months * C0_Q[j] + p.beta.cos[j] * C_cos[j] + p.beta.sin[j] * C_sin[j])
    
  }
}

for (j in 1:n.states){
  for (t in 1:n.months){
    predX[1,t,j] <- log(p[j, t]) + N[j,1]
    X[1,t,j] ~ dnorm(predX[1,t,j], sigma.Q[j,t])
    
    # observation
    predY[1,t,j] <- inprod(Z[j, 1:n.states], X[1,t,])
    
    #y[1,t,j] ~  dnorm(predY[1,t,j], tau.R[j])
    # observation
    y[1,t,j] ~ dt(predY[1,t,j], tau.R[j], df[j])
    
    loglik[1,t,j] <- logdensity.t(y[1,t,j], X[1,t,j], tau.R[j], df[j])
  }
}