#time series analysis


rm(list=ls())

tic <- Sys.time()
Sys <- Sys.info()
source('Dc_Indonesia_nesting_fcns.R')
library(jagsUI)
library(coda)

save.RData <- T
save.fig <- T


MCMC.n.chains <- 5
MCMC.n.samples <- 50000
MCMC.n.burnin <- 25000
MCMC.n.thin <- 5

MCMC.params <- list(n.chains = MCMC.n.chains,
                    n.samples = MCMC.n.samples,
                    n.burnin = MCMC.n.burnin,
                    n.thin = MCMC.n.thin)

# get JM data first:
data.0.JM <- read.csv('data/NestCounts_JM_09Feb2018.csv')

#data.0 <- read.csv("data/NestCount_Warmon_27March2018.csv")


# create time-duration filed (in yrs)
# define dates with begin and end dates:
data.0.JM %>% reshape2::melt(id.vars = "YEAR",
                             variable.name = "month",
                             value.name = "count") -> data.1.JM
data.1.JM$MONTH <- unlist(lapply(data.1.JM$month, FUN = mmm2month))

data.1.JM %>% mutate(., f.month = as.factor(MONTH),
                     f.year = as.factor(YEAR))%>%
  #filter(YEAR > 1999) %>%
  filter(YEAR > 2000) %>%
  #filter(YEAR > 1998) %>%
  mutate(Frac.Year = YEAR + (MONTH-0.5)/12) %>%
  reshape::sort_df(.,vars = "Frac.Year") -> data.2

# first three data points are NAs - remove
#data.2 <- data.2[4:nrow(data.2),]

#data.1.JM.2005 <- filter(data.1.JM, YEAR > 2004)

jags.data <- list(y = data.2$count,
                  m = data.2$MONTH,
                  T = nrow(data.2))

jags.params <- c('theta.1', 'theta.2',
                 'sigma.pro1', 'sigma.pro2',
                 'sigma.obs', 'mu', 'y', 'X', 'deviance')

jm <- jags(jags.data,
           inits = NULL,
           parameters.to.save= jags.params,
           model.file = 'models/model_SSAR1_month_var_theta.txt',
           n.chains = MCMC.n.chains,
           n.burnin = MCMC.n.burnin,
           n.thin = MCMC.n.thin,
           n.iter = MCMC.n.samples,
           DIC = T, parallel=T)

#g.diag1 <- gelman.diag(jm$samples)
Rhat <- jm$Rhat

# extract ys
ys.stats <- data.frame(low_y = jm$q2.5$y,
                       median_y = jm$q50$y,
                       high_y = jm$q97.5$y,
                       time = data.2$Frac.Year,
                       obsY = data.2$count,
                       month = data.2$MONTH,
                       year = data.2$YEAR)


# extract Xs - the state model
Xs.stats <- data.frame(low_X = jm$q2.5$X,
                       median_X = jm$q50$X,
                       high_X = jm$q97.5$X,
                       time = data.2$Frac.Year,
                       obsY = data.2$count,
                       month = data.2$MONTH,
                       year = data.2$YEAR)


Xs.year <- group_by(Xs.stats, year) %>% summarize(median = sum(median_X),
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
             aes(x = time, y = median_X), color = "red",
             alpha = 0.5, size = 2) +
  geom_line(data = Xs.stats,
            aes(x = time, y = median_X), color = "red",
            alpha = 0.5) +
  geom_line(data = Xs.stats,
            aes(x = time, y = low_X), color = "red",
            linetype = 2) +
  geom_point(data = ys.stats,
             aes(x = time, y = obsY), color = "green",
             alpha = 0.5) +
  labs(title = '', x = '', y = "Nest counts")

toc <- Sys.time()
dif.time <- toc - tic

results.JM_SSAR1_month_var_theta <- list(data.1 = data.2,
                                         jags.data = jags.data,
                                         Xs.stats = Xs.stats,
                                         Xs.year = Xs.year,
                                         ys.stats = ys.stats,
                                         tic = tic,
                                         toc = toc,
                                         dif.time = dif.time,
                                         Sys = Sys,
                                         MCMC.params = MCMC.params,
                                         Rhat = Rhat,
                                         jm = jm)

if (save.fig)
  ggsave(plot = p.1,
         filename = paste0('figures/predicted_counts_JM_month_var_May_Aug_theta_',
                           min(data.2$YEAR), '.png'),
         dpi = 600)

if (save.RData)
  save(results.JM_SSAR1_month_var_theta,
       file = paste0('RData/SSAR1_month_JM_var_May_Aug_theta_',
                     min(data.2$YEAR), '_',
                     Sys.Date(), '.RData'))

base_theme <- ggplot2::theme_get()
library(bayesplot)

# set back to the base theme:
ggplot2::theme_set(base_theme)
mcmc_trace(jm$samples, c('theta.1', "theta.2", 
                         "sigma.pro1", "sigma.pro2", "sigma.obs"))
mcmc_dens(jm$samples, c('theta.1', "theta.2", 
                        "sigma.pro1", "sigma.pro2", "sigma.obs"))

