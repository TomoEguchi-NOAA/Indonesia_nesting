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

jags.params <- c('N', 'theta', "p",
                 'sigma.pro', "sigma.obs",
                 "mu", "y", "X", "deviance", "loglik")

# when running with parallel=T, error returns...:
# Error in mcmc.list(x) : Different start, end or thin values in each chain

# jags.out <- run.jagsUI(data.jags$jags.data2, 
#                        jags.params, 
#                        model.file = 'models/model_SSAR1_logY_norm_norm_trend.txt', 
#                        MCMC.params)
# 
# Xs.stats <- jags.out$Xs.stats

# it runs fine with no-parallelized - slow but works. 
jm <- jags(data.jags$jags.data2,
           inits = NULL,
           parameters.to.save= jags.params,
           model.file = 'models/model_SSAR1_logY_norm_norm_trend.txt',
           n.chains = MCMC.params$n.chains,
           n.burnin = MCMC.params$n.burnin,
           n.thin = MCMC.params$n.thin,
           n.iter = MCMC.params$n.samples,
           DIC = T, parallel=F)

ys.stats <- data.frame(low_y = jm$q2.5$y,
                       median_y = jm$q50$y,
                       high_y = jm$q97.5$y)


# extract Xs - the state model
Xs.stats <- data.frame(low_X = jm$q2.5$X,
                       median_X = jm$q50$X,
                       high_X = jm$q97.5$X)

# stopped here.... 2019-04-09

# the following doesn't work 2019-04-09
loo.out <- pareto.k.diag(jm, MCMC.params, data.jags$jags.data2)

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
