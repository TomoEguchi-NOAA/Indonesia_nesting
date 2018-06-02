#time series analysis


rm(list=ls())

source('Dc_Indonesia_nesting_fcns.R')
library(rjags)
library(bayesplot)

n.chains <- 3
n.iter <- 10000
n.update <- 10000

# get JM data first:
data.0.JM <- read.csv('data/NestCounts_JM_09Feb2018.csv')

# create time-duration filed (in yrs)
# define dates with begin and end dates:
data.0.JM %>% reshape2::melt(id.vars = "YEAR",
                             variable.name = "month",
                             value.name = "count") -> data.1.JM
data.1.JM$MONTH <- unlist(lapply(data.1.JM$month, FUN = mmm2month))

data.1.JM <- mutate(data.1.JM, f.month = as.factor(MONTH),
                    f.year = as.factor(YEAR))

data.1.JM.2005 <- filter(data.1.JM, YEAR > 2004)

bugs.data <- list(y = data.1.JM.2005$count,
                  T = 156)

inits.function <- function(){
  eps <- rnorm(1, 0, 100)
  theta <- rnorm(1, 0, 100)
  c <- rnorm(1, 0, 100)
  A <- list(eps = eps, theta = theta, c = c)
  return(A)
}

params <- c('c', 'theta', 'sigma')

jm <- jags.model(file = 'models/model_AR1.txt',
                 data = bugs.data,
                 inits = inits.function,
                 n.chains = n.chains,
                 n.adapt = n.iter)

# check for convergence first.
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = n.iter)
gelman.diag(zm)
plot(zm)

# then sample y
params <- c('c', 'theta', 'sigma', 'y')
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = n.iter)

summary.zm <- summary(zm)

# plot posterior densities using bayesplot functions:
mcmc_dens(zm, "sigma")


# extract ys
ys.stats <- data.frame(summary.zm$quantiles[grep(pattern = 'y[/[]',
                                                 row.names(summary.zm$quantiles)),
                                            c('2.5%', '50%', '97.5%')])

colnames(ys.stats) <- c('lowN', 'modeN', 'highN')
ys.stats$time <- data.1.JM.2005$Frac.Year

plot(ys.stats$time, ys.stats$X50., col = 'red')

points(data.1.JM.2005$Frac.Year, data.1.JM.2005$count)


# SSAR1:
inits.function <- function(){
  #eps <- rnorm(1, 0, 100)
  theta <- rnorm(1, 0, 100)
  c <- rnorm(1, 0, 100)
  A <- list(theta = theta, c = c)
  return(A)
}

params <- c('c', 'theta', 'sigma')

jm <- jags.model(file = 'models/model_SSAR1.txt',
                 data = bugs.data,
                 inits = inits.function,
                 n.chains = n.chains,
                 n.adapt = n.iter)

zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = n.iter)
gelman.diag(zm)

plot(zm)

# redo it with ys:
params <- c('c', 'theta', 'sigma', 'y')

zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = n.iter)


summary.zm <- summary(zm)



# ARMA(2, 1):  Doesn't work! 3/16/2018
# bugs.data <- list(y = data.1.JM.2005$count,
#                   T = 156)

# inits.function <- function(){
#   phi <- rnorm(1, 0, 100)
#   #theta <- rnorm(1, 0, 100)
#   c <- rnorm(1, 0, 100)
#   A <- list(phi = phi, c = c)
#   return(A)
# }
#
# #params <- c('c', 'theta', 'sigma', 'y')
#
# jm <- jags.model(file = 'models/model_ARMA21.txt',
#                  data = bugs.data,
#                  inits = inits.function,
#                  n.chains = n.chains,
#                  n.adapt = n.iter)
#
# zm <- coda.samples(jm,
#                    variable.names = params,
#                    n.iter = n.iter)
# gelman.diag(zm)
#
# plot(zm)
# summary.zm <- summary(zm)
#
#
#
