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

# JM
year.begin <- 1999
year.end <- 2018
loc <- "JM"
period <- 12

# W
# year.begin <- 2006
# year.end <- 2018
# loc <- "W"
# period <- 6

data.jags <- data.extract(location = loc, 
                          year.begin = year.begin, 
                          year.end = year.end)

data.jags$jags.data$pi <- 3.141593
data.jags$jags.data$period <- period

jags.params <- c('beta.cos', 'beta.sin',
                 'sigma.pro1', "sigma.obs",
                 "mu", "y", "X", "deviance", "loglik")

jags.out <- run.jagsUI(data.jags$jags.data, 
                       jags.params, 
                       model.file = 'models/model_SSAR1_logY_norm_norm_varM_theta_Four_trend.txt', 
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

bayesplot::mcmc_trace(jags.out$jm$samples, "beta.cos")
bayesplot::mcmc_trace(jags.out$jm$samples, "beta.sin")

pareto.k <- jags.out$loo.out$loo.out$diagnostics$pareto_k
data.y <- na.omit(jags.out$jags.data$y)

pareto.df <- data.frame(y = data.y,
                        khat = pareto.k,
                        datapoint = seq(from = 1, to = length(data.y)),
                        k0.7 = cut(pareto.k,
                                   breaks = c(0, 0.7, 1.5),
                                   labels = c("<=0.7", ">0.7")))
ggplot(data = pareto.df) +   
  geom_path(aes(x = datapoint, y = exp(y)), alpha = 0.5) +
  geom_point(aes(x = datapoint, y = exp(y), 
                 size = khat,
                 color = k0.7)) +
  scale_size_continuous(limits = c(0.0, 1.3),
                        range = c(1, 4))+ 
  scale_color_manual(values = c("<=0.7" = "black", 
                                ">0.7" = "red")) 
