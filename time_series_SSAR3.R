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
                    f.year = as.factor(YEAR))%>%
  mutate(Frac.Year = YEAR + (MONTH-0.5)/12)

data.1.JM.2005 <- filter(data.1.JM, YEAR > 2004)

bugs.data <- list(y = data.1.JM.2005$count,
                  T = 156)

inits.function <- function(){
  mu <- rnorm(1, 0, 10)
  theta1 <- rnorm(1, 0, 1)
  theta2 <- rnorm(1, 0, 1)
  theta3 <- rnorm(1, 0, 1)
  #sigma.pro <- runif(1, 0, 50)
  #sigma.obs <- runif(1, 0, 50)
  A <- list(mu = mu, theta1 = theta1, theta2=theta2, theta3=theta3) 
  #          sigma.pro = sigma.pro, sigma.obs = sigma.obs)
  return(A)
}

params <- c('mu', 'theta1', 'theta2','theta3','sigma.pro', 'sigma.obs')

jm <- jags.model(file = 'models/model_SSAR3.txt',
                 data = bugs.data,
                 inits = inits.function,
                 n.chains = n.chains,
                 n.adapt = n.iter)

# check for convergence first.
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = n.iter)
gelman.diag(zm)

# plot posterior densities using bayesplot functions:
mcmc_dens(zm, 'theta1')
#mcmc_trace(zm, 'theta1')
mcmc_dens(zm, 'theta2')
#mcmc_trace(zm, 'theta2')
mcmc_dens(zm, 'theta3')
#mcmc_trace(zm, 'theta3')

mcmc_dens(zm, 'sigma.pro')
#mcmc_trace(zm, 'sigma')
mcmc_dens(zm, 'sigma.obs')

# then sample y
params <- c('theta1', 'theta2','theta3',
            'sigma.pro', 'sigma.obs', 'y', 'X', 'mu')
zm <- coda.samples(jm,
                   variable.names = params,
                   n.iter = n.iter)

summary.zm <- summary(zm)

# extract ys
ys.stats <- data.frame(summary.zm$quantiles[grep(pattern = 'y[/[]',
                                                 row.names(summary.zm$quantiles)),
                                            c('2.5%', '50%', '97.5%')])
colnames(ys.stats) <- c('low_y', 'mode_y', 'high_y')
ys.stats$time <- data.1.JM.2005$Frac.Year
ys.stats$obsY <- data.1.JM.2005$count

# extract Xs - the state model
Xs.stats <- data.frame(summary.zm$quantiles[grep(pattern = 'X[/[]',
                                                 row.names(summary.zm$quantiles)),
                                            c('2.5%', '50%', '97.5%')])
colnames(Xs.stats) <- c('low_X', 'mode_X', 'high_X')
Xs.stats$time <- data.1.JM.2005$Frac.Year
Xs.stats$obsY <- data.1.JM.2005$count

ggplot() +
  geom_point(data = Xs.stats,
             aes(x = time, y = mode_X), color = "blue") +
  geom_line(data = Xs.stats,
            aes(x = time, y = mode_X), color = 'blue') +
  geom_point(data = ys.stats,
             aes(x = time, y = mode_y), color = "red",
             alpha = 0.5) + 
  geom_line(data = ys.stats,
            aes(x = time, y = mode_y), color = "red",
            alpha = 0.5) + 
  geom_point(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5) 
