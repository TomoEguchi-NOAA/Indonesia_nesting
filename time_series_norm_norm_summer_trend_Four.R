#time series analysis


rm(list=ls())

# tic <- Sys.time()
# Sys <- Sys.info()
source('Dc_Indonesia_nesting_fcns.R')
#library(rjags)

library(jagsUI)
library(coda)
library(tidyverse)
library(loo)

fill.color <-  "darkseagreen"
fill.color.N <- "blue4"
fill.color.summer <- "darksalmon"
fill.color.winter <- "gray65"
fill.alpha <-  0.65
line.color <-  "darkblue"
line.color.N <- "cadetblue"
line.color.summer <- "red3"
line.color.winter <- "greenyellow"
data.color <- "black"
data.size <- 1.5
obsd.color <- "red2"

save.data <- T
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
# year.begin <- 1999
# year.end <- 2018
# loc <- "JM"
# loc.name <- "Jamursba-Medi"
# period <- 12
# maxN <- 10000

# W
year.begin <- 2003
year.end <- 2018
loc <- "W"
loc.name <- "Wermon"
period <- 6
maxN <- 3000

data.jags <- data.extract(location = loc, 
                          year.begin = year.begin, 
                          year.end = year.end)

jags.params <- c('N', 'theta', "p", "p.beta.cos", "p.beta.sin",
                 "sigma.N", 'sigma.pro', "sigma.obs",
                 "mu", "y", "X", "deviance", "loglik")

jags.data <- data.jags$jags.data.summer
jags.data$C0 <- 10
jags.data$n.months <- 6

jags.data$C_cos <- sum(apply(matrix(1:6, nrow=1), 
                             MARGIN = 1, 
                             FUN = function(x) cos(2 * pi * x/period)))

jags.data$C_sin <- sum(apply(matrix(1:6, nrow=1), 
                             MARGIN = 1, 
                             FUN = function(x) sin(2 * pi * x/period)))

jags.data$pi <- pi
jags.data$period <- period

# when running with parallel=T, error returns...:
# Error in mcmc.list(x) : Different start, end or thin values in each chain
# Restarting the computer fixed that problem next day... strange...

# it runs fine with no-parallelized - slow but works. 
jm <- jags(jags.data,
           inits = NULL,
           parameters.to.save= jags.params,
           model.file = 'models/model_norm_norm_trend_Four.txt',
           n.chains = MCMC.params$n.chains,
           n.burnin = MCMC.params$n.burnin,
           n.thin = MCMC.params$n.thin,
           n.iter = MCMC.params$n.samples,
           DIC = T, parallel=F)

# extract ys - include estimated missing data
# these need to be arranged in vectors
ys.stats <- data.frame(low_y = as.vector(t(jm$q2.5$y)),
                       median_y = as.vector(t(jm$q50$y)),
                       high_y = as.vector(t(jm$q97.5$y)))


# extract Xs - the state model
Xs.stats <- data.frame(low_X = as.vector(t(jm$q2.5$X)),
                       median_X = as.vector(t(jm$q50$X)),
                       high_X = as.vector(t(jm$q97.5$X)))

loo.out <- pareto.k.diag.3D(jm, MCMC.params, jags.data)

Xs.stats$time <- data.jags$data.summer$Frac.Year
Xs.stats$obsY <- data.jags$data.summer$Nests
Xs.stats$month <- data.jags$data.summer$Month
Xs.stats$year <- data.jags$data.summmer$Year

ys.stats$time <- data.jags$data.summer$Frac.Year
ys.stats$obsY <- data.jags$data.summer$Nests
ys.stats$month <- data.jags$data.summer$Month
ys.stats$year <- data.jags$data.summer$Year

Ns.stats <- data.frame(time = year.begin:year.end,
                       low_N = as.vector(t(jm$q2.5$N)),
                       median_N = as.vector(t(jm$q50$N)),
                       high_N = as.vector(t(jm$q97.5$N)))

p.1 <- ggplot() +
  geom_ribbon(data = Xs.stats,
              aes(x = time, 
                  ymin = exp(low_X), 
                  ymax = exp(high_X)),
              fill = fill.color,
              alpha = fill.alpha) +
  geom_point(data = Xs.stats,
             aes(x = time, 
                 y = exp(median_X)), 
             color = line.color,
             alpha = 0.5) +
  geom_line(data = Xs.stats,
            aes(x = time, 
                y = exp(median_X)), 
            color = line.color,
            alpha = 0.5) +
  geom_point(data = ys.stats,
             aes(x = time, y = obsY), 
             color = obsd.color,
             alpha = 0.5) + 
  geom_ribbon(data = Ns.stats,
              aes(x = time, 
                  ymin = exp(low_N),
                  ymax = exp(high_N)),
              fill = fill.color.N,
              alpha = fill.alpha) + 
  geom_line(data = Ns.stats,
            aes(x = time, y = exp(median_N)),
            color = line.color.N,
            alpha = 0.5,
            size = 1.5) + 
  scale_x_continuous(breaks = seq(year.begin, year.end, 5),
                     limits = c(year.begin, year.end)) +
  scale_y_continuous(limits = c(0, maxN)) + 
  labs(x = '', y = '# nests', title = loc.name)  +
  theme(axis.text = element_text(size = 12),
        text = element_text(size = 12))


p.1a <- ggplot() +
  geom_ribbon(data = Xs.stats,
              aes(x = time, 
                  ymin = low_X, 
                  ymax = high_X),
              fill = fill.color,
              alpha = fill.alpha) +
  geom_point(data = Xs.stats,
             aes(x = time, 
                 y = median_X), 
             color = line.color,
             alpha = 0.5) +
  geom_line(data = Xs.stats,
            aes(x = time, 
                y = median_X), 
            color = line.color,
            alpha = 0.5) +
  geom_point(data = ys.stats,
             aes(x = time, 
                 y = log(obsY)), 
             color = obsd.color,
             alpha = 0.5) + 
  geom_ribbon(data = Ns.stats,
              aes(x = time, 
                  ymin = low_N,
                  ymax = high_N),
              fill = fill.color.N,
              alpha = fill.alpha) + 
  geom_line(data = Ns.stats,
            aes(x = time, y = median_N),
            color = line.color.N,
            alpha = 0.5,
            size = 1.5) + 
  scale_x_continuous(breaks = seq(year.begin, year.end, 5),
                     limits = c(year.begin, year.end)) +
  scale_y_continuous(limits = c(0, log(maxN))) + 
  labs(x = '', y = 'log(# nests)', title = loc.name)  +
  theme(axis.text = element_text(size = 12),
        text = element_text(size = 12))

#p.1a

bayesplot::mcmc_trace(jm$samples, "theta")

bayesplot::mcmc_trace(jm$samples, "p.beta.cos")
bayesplot::mcmc_trace(jm$samples, "p.beta.sin")

bayesplot::mcmc_dens(jm$samples, "theta")
bayesplot::mcmc_dens(jm$samples, "p.beta.cos")
bayesplot::mcmc_dens(jm$samples, "p.beta.sin")


pareto.k <- loo.out$loo.out$diagnostics$pareto_k
data.y <- na.omit(as.vector(t(jags.data$y))

pareto.df <- data.frame(y = data.y,
                        khat = pareto.k,
                        datapoint = seq(from = 1, to = length(data.y)),
                        k0.7 = cut(pareto.k,
                                   breaks = c(0, 0.7, 1.5),
                                   labels = c("<=0.7", ">0.7")))
p.2 <- ggplot(data = pareto.df) +   
  geom_path(aes(x = datapoint, y = exp(y)), alpha = 0.5) +
  geom_point(aes(x = datapoint, y = exp(y), 
                 size = khat,
                 color = k0.7)) +
  scale_size_continuous(limits = c(0.0, 1.3),
                        range = c(1, 4))+ 
  scale_color_manual(values = c("<=0.7" = "black", 
                                ">0.7" = "red")) 

results.all <- list(jm = jm,
                    Xs.stats = Xs.stats,
                    ys.stats = ys.stats,
                    Ns.stats = Ns.stats)
if (save.data)
  saveRDS(results.all,
          file = paste0('RData/SSAR1_norm_norm_summer_trend_Four_', loc, "_",
                        year.begin, "_", year.end, "_", Sys.Date(), '.rds'))

if (save.fig){
  ggsave(filename = paste0('figures/SSAR1_norm_norm_summer_trend_Four_', loc, "_",
                           year.begin, "_", year.end, "_", Sys.Date(), ".png"),
         plot = p.1,
         dpi = 600,
         device = "png")
  
  ggsave(filename = paste0('figures/SSAR1_norm_norm_log_summer_trend_Four_', loc, "_",
                           year.begin, "_", year.end, "_", Sys.Date(), ".png"),
         plot = p.1a,
         dpi = 600,
         device = "png")

    ggsave(filename = paste0('figures/SSAR1_norm_norm_summer_trend_Four_', loc, "_",
                           year.begin, "_", year.end, "_", Sys.Date(), "_pareto.png"),
         plot = p.2,
         dpi = 600,
         device = "png")
  
}
