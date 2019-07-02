# data simulation script

# This script simulates data on nesting beach counts of leatherbacks. Simulated data are analyzed
# using various models to see how imputations and Pareto - k statistics perform.

rm(list=ls())

library(ggplot2)
library(tidyverse)

# create a log density function for t because R doesn't give the same definition as JAGS:
logdensity.t <- function(y, mu, tau, k){
  log.t <- log(gamma((k + 1)/2)) - log(gamma(k/2)) + (1/2) * log(tau/(k * pi)) - ((k + 1)/2) * log(1 + (tau * (y - mu)^2)/k)
  return(log.t)
}


set.seed(123)

# parameters
TT <- 17  # years of time series data
n.states <- 2   # # nesting beaches
n.months <- 12  # # months
periods <- c(12, 6)  # periods of nesting cycles

N0.mean <- c(8.1, 7.9)   # initial values for the abundance mean in log scale
N0_sd <- c(1.4, 1.2)       # SD of the initial N in log scale

# multiple U
# matrix with rows as n.timeseries and cols as n.states (pops)
Z <- matrix(0, n.states, n.states)   

# # add a row of NAs to keep jagsUI from converting single time series matrix into vector
# Z[n.timeseries+1, ] <- NA                  
# # add a col of NAs to keep jagsUI from converting single state matrix into vector
# Z[ , n.states+1] <- NA                     
for(i in 1:n.states) Z[i, i] <- 1

U.mean <- c(-0.12, -0.10)    # annual growth rates for the two population
U.sd <- c(0.15, 0.15)        # SD of the growth rates

sigma.N <- 0.6   # SD of N

N0 <- numeric(length = n.states)
N <- predN <- array(dim = c(n.states, TT))

# Initital states
for (j in 1:n.states){
  N0[j] <- rnorm(n = 1, mean = N0.mean[j], sd = N0_sd[j])

}

C0_Q <- C0 <- c(15, 15)

p.beta.cos <- c(-1.27, -4.08)
p.beta.sin <- c(4.18, -2.81)

C_cos <- c(sum(apply(matrix(1:12, nrow=1), 
                    MARGIN = 1, 
                    FUN = function(x) cos(2 * pi * x/periods[1]))), 
          sum(apply(matrix(1:12, nrow=1), 
                    MARGIN = 1, 
                    FUN = function(x) cos(2 * pi * x/periods[2]))))

C_sin <- c(sum(apply(matrix(1:12, nrow=1), 
                    MARGIN = 1, 
                    FUN = function(x) sin(2 * pi * x/periods[1]))),
          sum(apply(matrix(1:12, nrow=1), 
                    MARGIN = 1, 
                    FUN = function(x) sin(2 * pi * x/periods[2]))))


sigma.R <- c(0.96, 0.61)
df <- c(49.9, 39.8)
cv.Q <- c(0.08, 0.02)

p.const <- p <- sigma.Q <- array(dim = c(n.states, n.months))
for (j in 1:n.states){
  for (k in 1:n.months){
    p.const[j, k] <-  2 * pi * k / periods[j]
    p[j, k] <- (C0[j] + p.beta.cos[j] * cos(p.const[j,k]) + p.beta.sin[j] * sin(p.const[j,k]))/(n.months * C0[j] + p.beta.cos[j] * C_cos[j] + p.beta.sin[j] * C_sin[j])
    
    # This is for sigma-Q to be Fourier function... 
    #sigma.Q[j,k] <- (C0_Q[j] + p.beta.cos[j] * cos(p.const[j,k]) + p.beta.sin[j] * sin(p.const[j,k]))/(n.months * C0_Q[j] + p.beta.cos[j] * C_cos[j] + p.beta.sin[j] * C_sin[j])
    
    
  }
}

predX <- X <- predY <- y.norm <- y.t <- loglik.norm <- loglik.t <- sigma.Q <- array(dim = c(TT, n.months, n.states))
tt <- 2
j <- k <- 1
for (tt in 1:TT){
  for (j in 1:n.states){

    # state
    if (tt == 1){
      predN[j, tt] <- N0[j] + rnorm(n = 1, mean = U.mean[j], sd = U.sd[j])
    } else {
      predN[j, tt] <- rnorm(n = 1, mean = U.mean[j], sd = U.sd[j]) + N[j, tt-1]
    }
    
    # Total N for the year
    N[j, tt] <- rnorm(n = 1, mean = predN[j, tt], sd = sigma.N) 
    
    for (k in 1:n.months){
      predX[tt, k, j] <- log(p[j, k]) + N[j, tt]
      
      # This is for constant CV model
      sigma.Q[tt,k,j] <- predX[tt, k, j] * cv.Q[j]
      
      X[tt, k, j] <- rnorm(n = 1, mean = predX[tt, k, j], sd = sigma.Q[tt, k, j])
      
      #predY[tt, k, j] <- Z[j, 1:n.states] %*% X[1,k,]

      # observation - normal
      y.norm[tt, k, j] <- rnorm(n = 1, mean = X[tt, k, j], sd = sigma.R[j])
      
      # observation - t - note rt() does not have sigma.R component... 
      y.t[tt, k, j] <- rt(n = 1, ncp = X[tt, k, j], df = df[j])
      
      loglik.norm[tt, k, j] <- dnorm(y.norm[tt, k, j], 
                                     mean = X[tt, k, j], 
                                     sd = sigma.R[j], log = TRUE)
      
      loglik.t[tt, k, j] <- logdensity.t(y.t[tt, k, j], 
                                         mu = X[tt, k, j], 
                                         tau = 1/((sigma.R[j])^2), 
                                         k = df[j])
    }
  }  
}

# observed data for locations 1 and 2:
y.1.norm.melt <- reshape::melt(t(y.norm[,,1]))
y.2.norm.melt <- reshape::melt(t(y.norm[,,2]))

y.1.t.melt <- reshape::melt(t(y.t[,,1]))
y.2.t.melt <- reshape::melt(t(y.t[,,2]))

ys.df <- data.frame(year = y.1.norm.melt[,2],
                    month = y.1.norm.melt[,1],
                    y.1.norm = y.1.norm.melt[,3],
                    y.2.norm = y.2.norm.melt[,3],
                    y.1.t = y.1.t.melt[,3],
                    y.2.t = y.2.t.melt[,3])

params.list <- list(N0.mean = N0.mean,
                    N0.sd = N0_sd,
                    U.mean = U.mean,
                    U.sd = U.sd,
                    N.sd = sigma.N,
                    p.beta.cos = p.beta.cos,
                    p.beta.sin = p.beta.sin,
                    sigma.R = sigma.R,
                    df = df,
                    cv.Q = cv.Q)

out.list <- list(data = ys.df,
                 parameters = params.list)

saveRDS(out.list, "RData/sim_constCV_independentUQ_data_parameters.rds")
