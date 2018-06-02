#time series analysis


rm(list=ls())

tic <- Sys.time()
Sys <- Sys.info()
source('Dc_Indonesia_nesting_fcns.R')
library(rjags)
library(bayesplot)

save.RData <- T
save.fig <- F

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
  mutate(Frac.Year = YEAR + (MONTH-0.5)/12) %>%
  mutate(log.count = log(count)) %>%
  reshape::sort_df(.,vars = "Frac.Year")

data.1.JM.2005 <- filter(data.1.JM, YEAR > 2004)

bugs.data <- list(y = data.1.JM.2005$log.count,
                  m = data.1.JM.2005$MONTH,
                  T = 156)

inits.function <- function(){
  mu <- rnorm(1, 0, 10)
  theta <- rnorm(1, 0, 1)
  phi <- rnorm(1, 0, 1)
  #sigma.pro <- runif(1, 0, 50)
  #sigma.obs <- runif(1, 0, 50)
  A <- list(mu = mu, theta = theta)
  #          sigma.pro = sigma.pro, sigma.obs = sigma.obs)
  return(A)
}

load.module('dic')
params <- c('theta', 'sigma.pro1', 'sigma.pro2',
            'sigma.obs', 'mu')

jm <- jags.model(file = 'models/model_SSAR1_month.txt',
                 data = bugs.data,
                 inits = inits.function,
                 n.chains = MCMC.params$n.chains,
                 n.adapt = MCMC.params$n.iter)

# check for convergence first.
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = MCMC.params$n.iter)
g.diag <- gelman.diag(zm)

# plot posterior densities using bayesplot functions:
#mcmc_dens(zm, 'theta')
#mcmc_trace(zm, 'theta')
#mcmc_dens(zm, 'phi1')
#mcmc_dens(zm, 'phi2')
#mcmc_trace(zm, 'phi2')
#mcmc_dens(zm, 'sigma.pro1')
#mcmc_dens(zm, 'sigma.pro2')
#mcmc_trace(zm, 'sigma')
#mcmc_dens(zm, 'sigma.obs')

# then sample y
params <- c(params, 'y', 'X', 'deviance')
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
            aes(x = time, y = exp(high_X)), color = "red",
            linetype = 2) +
  geom_point(data = Xs.stats,
             aes(x = time, y = exp(mode_X)), color = "red",
             alpha = 0.5) +
  geom_line(data = Xs.stats,
            aes(x = time, y = exp(mode_X)), color = "red",
            alpha = 0.5) +
  geom_point(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5)

toc <- Sys.time()
dif.time <- toc - tic
results.JM_SSAR1_logN_month <- list(data.1 = data.1.JM.2005,
                                    summary.zm = summary.zm,
                                    Xs.stats = Xs.stats,
                                    Xs.year = Xs.year,
                                    ys.stats = ys.stats,
                                    zm = zm,
                                    tic = tic,
                                    toc = toc,
                                    dif.time = dif.time,
                                    Sys = Sys,
                                    MCMC.params = MCMC.params,
                                    g.diag = g.diag,
                                    jm = jm)
if (save.fig)
  ggsave(plot = p.1,
         filename = 'figures/predicted_logN_SSAR1month_JM.png',
         dpi = 600)

if (save.RData)
  save(results.JM_SSAR1_logN_month,
       file = paste0('RData/SSAR1_logN_month_', Sys.Date(), '.RData'))
