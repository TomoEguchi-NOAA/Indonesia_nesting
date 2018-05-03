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

data.0 <- read.csv("data/NestCounts_Warmon_27March2018.csv")

# create time-duration filed (in yrs)
# define dates with begin and end dates:
data.0 %>% reshape2::melt(id.vars = "YEAR",
                             variable.name = "month",
                             value.name = "count") -> data.1
data.1$MONTH <- unlist(lapply(data.1$month, FUN = mmm2month))

data.1 <- mutate(data.1, f.month = as.factor(MONTH),
                    f.year = as.factor(YEAR))%>%
  mutate(Frac.Year = YEAR + (MONTH-0.5)/12) %>%
  reshape::sort_df(.,vars = "Frac.Year")

bugs.data <- list(y = data.1$count,
                  m = data.1$MONTH,
                  T = 168)

# bugs.data <- list(y = data.1$count,
#                   T = 168)

inits.function <- function(){
  mu <- rnorm(1, 0, 10)
  theta1 <- rnorm(1, 0, 10)
  theta2 <- rnorm(1, 0, 10)
  #phi <- rnorm(1, 0, 1)
  #sigma.pro <- runif(1, 0, 50)
  #sigma.obs <- runif(1, 0, 50)
  A <- list(mu = mu, theta1 = theta1, theta2 = theta2)
  #          sigma.pro = sigma.pro, sigma.obs = sigma.obs)
  return(A)
}

load.module('dic')
params <- c('theta1', 'theta2', 'sigma.pro1', 'sigma.pro2',
            'sigma.obs', 'mu', 'predY')

jm <- jags.model(file = 'models/model_SSAR2_month_Warmon.txt',
                 data = bugs.data,
                 #inits = inits.function,
                 n.chains = MCMC.params$n.chains,
                 n.adapt = MCMC.params$n.iter)

# check for convergence first.
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = MCMC.params$n.iter)

# I think theta1 and theta2 are indistinguishable and the samples turned out
# to be non positive definite - can't compute gelman diagnostic for convergence.
# This probably means that AR(1) is a better pick.

g.diag <- gelman.diag(zm)

# plot posterior densities using bayesplot functions:
# mcmc_dens(zm, 'theta')
# mcmc_dens(zm, 'sigma.pro')
# mcmc_dens(zm, 'sigma.obs')

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
ys.stats$time <- data.1$Frac.Year
ys.stats$obsY <- data.1$count
ys.stats$month <- data.1$MONTH
ys.stats$year <- data.1$YEAR

# extract Xs - the state model
Xs.stats <- data.frame(summary.zm$quantiles[grep(pattern = 'X[/[]',
                                                 row.names(summary.zm$quantiles)),
                                            c('2.5%', '50%', '97.5%')])
colnames(Xs.stats) <- c('low_X', 'mode_X', 'high_X')
Xs.stats$time <- data.1$Frac.Year
Xs.stats$obsY <- data.1$count
Xs.stats$month <- data.1$MONTH
Xs.stats$year <- data.1$YEAR

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
             alpha = 0.5)+
  geom_line(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5)

toc <- Sys.time()
dif.time <- toc - tic

results.Warmon_SSAR1_month <- list(data.1 = data.1,
                                   summary.zm = summary.zm,
                                   Xs.stats = Xs.stats,
                                   Xs.year = Xs.year,
                                   ys.stats = ys.stats,
                                   jm = jm,
                                   zm = zm,
                                   tic = tic,
                                   toc = toc,
                                   dif.time = dif.time,
                                   Sys = Sys,
                                   MCMC.params = MCMC.params,
                                   g.diag = g.diag)
if (save.fig)
  ggsave(plot = p.1,
         filename = 'figures/predicted_counts_SSAR2_Warmon.png',
         dpi = 600)

if (save.RData)
  save(results.Warmon_SSAR1_month,
       file = paste0('RData/SSAR2_month_Warmon_', Sys.Date(), '.RData'))
