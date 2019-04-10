#time series analysis


rm(list=ls())

tic <- Sys.time()
Sys <- Sys.info()
source('Dc_Indonesia_nesting_fcns.R')
#library(rjags)

library(jagsUI)
library(coda)
# library(ggplot2)
library(loo)

save.RData <- T
save.fig <- T

MCMC.n.chains <- 5
MCMC.n.samples <- 500000
MCMC.n.burnin <- 350000
MCMC.n.thin <- 50

MCMC.params <- list(n.chains = MCMC.n.chains,
                    n.samples = MCMC.n.samples,
                    n.burnin = MCMC.n.burnin,
                    n.thin = MCMC.n.thin)

year.begin <- 1999
year.end <- 2018
loc <- "JM"
data.jags <- data.extract(location = loc, 
                          year.begin = year.begin, 
                          year.end = year.end)

data.jags$jags.data$pi <- 3.141593
data.jags$jags.data$period <- 12

jags.params <- c('beta.cos', 'beta.sin',
                 'sigma.pro1', "sigma.obs",
                 "mu", "y", "X", "deviance", "loglik")

jags.out <- run.jagsUI(data.jags$jags.data, 
                       jags.params, 
                       model.file = 'models/model_SSAR1_logY_norm_norm_varM_theta_Four.txt', 
                       MCMC.params)

Xs.stats <- jags.out$Xs.stats

Xs.stats$time <- data.jags$data.1$Frac.Year
Xs.stats$obsY <- data.jags$data.1$Nests
Xs.stats$month <- data.jags$data.1$Month
Xs.stats$year <- data.jags$data.1$Year

ys.stats <- jags.out$ys.stats
ys.stats$time <- data.jags$data.1$Frac.Year
ys.stats$obsY <- data.jags$data.1$Nests
ys.stats$month <- data.jags$data.1$Month
ys.stats$year <- data.jags$data.1$Year


p.1 <- ggplot() +
  #geom_point(data = ys.stats,
  #           aes(x = time, y = mode_y), color = "blue") +
  #geom_line(data = Xs.stats,
  #          aes(x = time, y = mode_X), color = 'blue') +
  geom_line(data = Xs.stats,
            aes(x = time, y = exp(high_X)), color = "red",
            linetype = 2) +
  geom_point(data = Xs.stats,
             aes(x = time, y = exp(median_X)), color = "red",
             alpha = 0.5) +
  geom_line(data = Xs.stats,
            aes(x = time, y = exp(median_X)), color = "red",
            alpha = 0.5) +
  geom_line(data = Xs.stats,
            aes(x = time, y = exp(low_X)), color = "red",
            linetype = 2) +
  geom_point(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5)


