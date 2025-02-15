#time series analysis


rm(list=ls())

source('Dc_Indonesia_nesting_fcns.R')
library(rjags)
library(bayesplot)

Sys <- Sys.info()
tic <- Sys.time()

save.RData <- T
save.fig <- T

MCMC.params <- list(n.chains = 3,
                    n.iter = 50000)

# get JM data first:
data.0.JM <- read.csv('data/NestCounts_JM_09Feb2018.csv')

# create time-duration filed (in yrs)
# define dates with begin and end dates:
data.0.JM %>% reshape2::melt(id.vars = "YEAR",
                             variable.name = "month",
                             value.name = "count") -> data.1.JM
data.1.JM$MONTH <- unlist(lapply(data.1.JM$month, FUN = mmm2month))

data.1.JM <- mutate(data.1.JM, f.month = as.factor(MONTH),
                    f.year = as.factor(YEAR))%>%
  mutate(Frac.Year = YEAR + (MONTH-0.5)/12)

data.1.JM.2005 <- filter(data.1.JM, YEAR > 2004)

bugs.data <- list(y = data.1.JM.2005$count,
                  T = 156)

inits.function <- function(){
  mu <- rnorm(1, 0, 10)
  theta <- rnorm(1, 0, 1)
  #sigma.pro <- runif(1, 0, 50)
  #sigma.obs <- runif(1, 0, 50)
  A <- list(mu = mu, theta = theta)
  #          sigma.pro = sigma.pro, sigma.obs = sigma.obs)
  return(A)
}

load.module('dic')
params <- c('theta', 'sigma.pro', 'sigma.obs')

jm <- jags.model(file = 'models/model_SSAR1.txt',
                 data = bugs.data,
                 #inits = inits.function,
                 n.chains = MCMC.params$n.chains,
                 n.adapt = MCMC.params$n.iter)

# check for convergence first.
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = MCMC.params$n.iter)
g.diag <- gelman.diag(zm)

# plot posterior densities using bayesplot functions:
mcmc_dens(zm, 'theta')
#mcmc_trace(zm, 'theta')
mcmc_dens(zm, 'sigma.pro')
#mcmc_trace(zm, 'sigma')
mcmc_dens(zm, 'sigma.obs')

# then sample y
params <- c(params, 'deviance', 'y', 'X')
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = MCMC.params$n.iter)

summary.zm <- summary(zm)

# extract ys
ys.stats <- data.frame(summary.zm$quantiles[grep(pattern = 'y[/[]',
                                                 row.names(summary.zm$quantiles)),
                                            c('2.5%', '50%', '97.5%')])
colnames(ys.stats) <- c('low_y', 'mode_y', 'high_y')
ys.stats$time <- data.1.JM.2005$Frac.Year
ys.stats$obsY <- data.1.JM.2005$count
ys.stats$month <- data.1.JM.2005$MONTH
ys.stats$year <- data.1.JM.2005$YEAR

# extract Xs - the state model
Xs.stats <- data.frame(summary.zm$quantiles[grep(pattern = 'X[/[]',
                                                 row.names(summary.zm$quantiles)),
                                            c('2.5%', '50%', '97.5%')])
colnames(Xs.stats) <- c('low_X', 'mode_X', 'high_X')
Xs.stats$time <- data.1.JM.2005$Frac.Year
Xs.stats$obsY <- data.1.JM.2005$count
Xs.stats$month <- data.1.JM.2005$MONTH
Xs.stats$year <- data.1.JM.2005$YEAR

Xs.year <- group_by(Xs.stats, year) %>% summarize(mode = sum(mode_X),
                                                  low = sum(low_X),
                                                  high = sum(high_X))

p.1 <- ggplot() +
  #geom_point(data = ys.stats,
  #           aes(x = time, y = mode_y), color = "blue") +
  #geom_line(data = Xs.stats,
  #          aes(x = time, y = mode_X), color = 'blue') +
  geom_line(data = Xs.stats,
            aes(x = time, y = high_X), color = "red",
            linetype = 2) +
  geom_point(data = Xs.stats,
             aes(x = time, y = mode_X), color = "red",
             alpha = 0.5) +
  geom_line(data = Xs.stats,
            aes(x = time, y = mode_X), color = "red",
            alpha = 0.5) +
  geom_point(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5)
p.1
toc <- Sys.time()
dif.time <- toc - tic

if (save.fig)
  ggsave(plot = p.1,
         filename = 'figures/predicted_counts_SSAR1_JM.png',
         dpi = 600)

if (save.RData)
  save(data.1.JM.2005, summary.zm, Xs.stats, Xs.year, ys.stats, zm,
       tic, toc, dif.time, Sys, MCMC.params, g.diag,
       file = paste0('RData/SSAR1_', Sys.Date(), '.RData'))
