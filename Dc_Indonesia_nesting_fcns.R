#Dc_Indonesia_nesting_fcns

# these are functions that are used in this project.
# TomosFunctions.R is needed also.

#sysInfo <- Sys.info()
ifelse(Sys.info()[1] == 'Linux',
       source('~/Documents/R/tools/TomosFunctions.R'),
       source('~/R/tools/TomosFunctions.R'))

library(ggplot2)
library(tidyverse)
library(lubridate)
library(loo)
library(dplyr)

save.fig <- F

sum.posterior <- function(yr, months = c(1:12), Xs.stats, zm) {
  Xs.stats %>%
    mutate(var.name = rownames(Xs.stats)) %>%
    filter(year == yr) %>%
    filter(month %in% months) %>%
    select(var.name) -> Xs.name
  zm.yr <- apply(Xs.name,
                 MARGIN = 1,
                 FUN = extract.samples,
                 zm)

  return(list(samples = zm.yr, var.names = Xs.name))
}

data.extract <- function(location, year.begin, year.end){
  if (location == "JM"){
    data.0 <- read.csv("data/JM_nests_March2019.csv")
    data.0 %>% mutate(Nests = JM.1) -> data.0
  } else if (location == "W"){
    data.0 <- read.csv('data/W_nests_March2019.csv')
    data.0 %>% mutate(Nests = W.1) -> data.0
  }
  # create regularly spaced time series:
  data.2 <- data.frame(Year = rep(min(data.0$Year_begin,
                                      na.rm = T):max(data.0$Year_begin,
                                                     na.rm = T),
                                  each = 12),
                       Month_begin = rep(1:12,
                                         max(data.0$Year_begin,
                                             na.rm = T) -
                                           min(data.0$Year_begin,
                                               na.rm = T) + 1)) %>%
    mutate(begin_date = as.Date(paste(Year,
                                      Month_begin,
                                      '01', sep = "-"),
                                format = "%Y-%m-%d"),
           Frac.Year = Year + (Month_begin-0.5)/12) %>%
    select(Year, Month_begin, begin_date, Frac.Year)
  
  data.0 %>% mutate(begin_date = as.Date(paste(Year_begin,
                                               Month_begin,
                                               '01', sep = "-"),
                                         format = "%Y-%m-%d")) %>%
    mutate(Year = Year_begin,
           Month = Month_begin,
           f_month = as.factor(Month),
           f_year = as.factor(Year),
           Frac.Year = Year + (Month_begin-0.5)/12) %>%
    select(Year, Month, Frac.Year, begin_date, Nests) %>%
    na.omit() %>%
    right_join(.,data.2, by = "begin_date") %>%
    transmute(Year = Year.y,
              Month = Month_begin,
              Frac.Year = Frac.Year.y,
              Nests = Nests) %>%
    reshape::sort_df(.,vars = "Frac.Year") %>%
    filter(Year >= year.begin & Year <= year.end) -> data.1
  
  jags.data <- list(y = log(data.1$Nests),
                    m = data.1$Month,
                    T = nrow(data.1))
  
  out <- list(jags.data = jags.data,
              data.1 = data.1)
  return(out)
}

run.jagsUI <- function(jags.data, jags.params, model.file, MCMC.params = list(n.chains = 5,
                                                                              n.samples = 500000,
                                                                              n.burnin = 350000,
                                                                              n.thin = 50)){
  tic <- Sys.time()
  jm <- jags(jags.data,
             inits = NULL,
             parameters.to.save= jags.params,
             model.file = model.file,
             n.chains = MCMC.params$n.chains,
             n.burnin = MCMC.params$n.burnin,
             n.thin = MCMC.params$n.thin,
             n.iter = MCMC.params$n.samples,
             DIC = T, parallel=T)
  
  # extract ys
  ys.stats <- data.frame(low_y = jm$q2.5$y,
                         median_y = jm$q50$y,
                         high_y = jm$q97.5$y)
  
  
  # extract Xs - the state model
  Xs.stats <- data.frame(low_X = jm$q2.5$X,
                         median_X = jm$q50$X,
                         high_X = jm$q97.5$X)
  
  #Xs.year <- group_by(Xs.stats, year) %>% summarize(median = sum(median_X),
  #                                                  low = sum(low_X),
  #                                                  high = sum(high_X))
  
  loo.out <- pareto.k.diag(jm, MCMC.params, jags.data)
  
  toc <- Sys.time()
  dif.time <- toc - tic
  
  results <- list(jags.data = jags.data,
                  Xs.stats = Xs.stats,
                  ys.stats = ys.stats,
                  tic = tic,
                  toc = toc,
                  dif.time = dif.time,
                  Sys = Sys.info(),
                  MCMC.params = MCMC.params,
                  jm = jm,
                  loo.out = loo.out)
  return(results)
}


extract.posterior.jagsUI <- function(yr, Xs.stats, samples){
  Xs.stats %>%
    mutate(var.name = rownames(Xs.stats)) %>%
    filter(season == yr) %>%
    select(var.name, summer) -> Xs.name
  
  all.samples <- do.call(rbind, samples)
  Xnames <- apply(as.matrix(Xs.name$var.name), 
                  FUN = function(x) paste0("X[", x, "]"), 
                  MARGIN = 1)
  
  out.samples <- all.samples[, Xnames]
  
  return(list(samples = out.samples, var.names = Xnames, summer = Xs.name$summer))
  
}

pareto.k.diag <- function(jm, MCMC.params, jags.data){
  
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  loglik.obs <- jm$sims.list$loglik[, !is.na(jags.data$y)]
  # get rid of NA columns - even if data existed (for example the first value) - no likelihood
  # for the first data point
  loglik.obs <- loglik.obs[, colSums(is.na(loglik.obs)) == 0]
  
  #loglik.obs <- jm$sims.list$loglik[, 2:jags.data$T]
  # cores = 1 is needed in the relative_eff function if the number of cores was set to more than
  # 1 with options(mc.cores = parallel::detectCores()) or something similear. See also here:
  # https://discourse.mc-stan.org/t/error-in-loo-relative-eff-related-to-options-mc-cores/5610/2
  
  Reff <- relative_eff(exp(loglik.obs), 
                       chain_id = rep(1:MCMC.params$n.chains, 
                                      each = n.per.chain),
                       cores = 1)
  
  loo.out <- loo(loglik.obs, r_eff = Reff, cores = 1)
  return(list(loglik.obs = loglik.obs,
              Reff = Reff,
              loo.out = loo.out))
}

model.names <- function(){
  norm.norm.models <- c("model_SSAR1_logY_norm_norm.txt",
                        "model_SSAR1_logY_norm_norm_theta.txt",
                        "model_SSAR1_logY_norm_norm_var.txt",
                        "model_SSAR1_logY_norm_norm_var_theta.txt",
                        "model_SSAR1_logY_norm_norm_varM_thetaM.txt",
                        "model_SSAR1_logY_norm_norm_varM_theta.txt",
                        "model_SSAR1_logY_norm_norm_var_thetaM.txt",
                        "model_SSAR1_logY_norm_norm_thetaM.txt",
                        "model_SSAR1_logY_norm_norm_varM.txt")
  
  norm.t.models <- c("model_SSAR1_logY_norm_t.txt",
                     "model_SSAR1_logY_norm_t_theta.txt",
                     "model_SSAR1_logY_norm_t_var.txt",
                     "model_SSAR1_logY_norm_t_var_theta.txt",
                     "model_SSAR1_logY_norm_t_varM_thetaM.txt",
                     "model_SSAR1_logY_norm_t_varM_theta.txt",
                     "model_SSAR1_logY_norm_t_var_thetaM.txt",
                     "model_SSAR1_logY_norm_t_thetaM.txt",
                     "model_SSAR1_logY_norm_t_varM.txt")
  
  all.models <- c(norm.norm.models, norm.t.models)
  
  models <- data.frame(number = paste0("M", 1:length(all.models)),
                       names = all.models)
  
  return(models)
  
}

run.singleUQ <- function(dat, 
                         jags.params, 
                         model.loc, 
                         MCMC.params, 
                         seed){
  
  n.timeseries= nrow(dat)
  # Model-specific parameters
  # multiple time series -> single population process
  n.states <- 1
  whichPop <- rep(1, n.timeseries)
  
  Z <- matrix(0, n.timeseries+1, n.states+1)   # matrix with rows as n.timeseries and cols as n.states (pops)
  Z[n.timeseries+1, ] <- NA                  # add a row of NAs to keep jagsUI from converting single time series matrix into vector
  Z[ , n.states+1] <- NA                     # add a col of NAs to keep jagsUI from converting single state matrix into vector
  for(i in 1:n.timeseries) Z[i, whichPop[i]] <- 1
  
  # add a row (ie time series of data) to trick jagsUI into NOT converting single time series matrix to vector
  jags.data <- list(Y = rbind(dat, NA), 
                    n.yrs = ncol(dat),
                    n.timeseries= n.timeseries,
                    Z = Z,
                    a_mean = 0,
                    a_sd = 4,
                    u_mean = 0,
                    u_sd = 0.5,
                    q_alpha = 0.01,
                    q_beta = 0.01,
                    r_alpha = 0.01,
                    r_beta = 0.01,
                    x0_mean = dat[1,1],
                    x0_sd = 10)
  
  jags.model <- jags(jags.data, 
                     inits = NULL, 
                     parameters.to.save= jags.params, 
                     model.file=model.loc, 
                     n.chains = MCMC.params$n.chains, 
                     n.burnin = MCMC.params$n.burnin, 
                     n.thin = MCMC.params$n.thin, 
                     n.iter = MCMC.params$n.samples, 
                     DIC = T, 
                     parallel=T, 
                     seed=seed)
  
  # not sure how pareto statistic can be computed for multi-location data
  # with missing data... 
  #loo.out <- pareto.k.diag(jags.model, MCMC.params, jags.data)
  
  jm.out <- list(jags.data = jags.data,
                 params = jags.params,
                 out = jags.model)
  
  return(jm.out)
  
}

run.independentUQ <- function(dat, 
                              jags.params, 
                              model.loc, 
                              MCMC.params, 
                              seed){
  
  n.timeseries <- nrow(dat)
  x0_mean <- numeric(length=n.timeseries)
  x0_sd <- numeric(length=n.timeseries)
  for(i in 1:n.timeseries) {
    icol <- min(which(!is.na(dat[i,])))   # use the first non-NA data point in time series i as it's prior mean
    x0_mean[i] <- dat[i, icol]            
    x0_sd[i] <- 10     # we went with wide sd by testing for both JM and W for leathers; 10 works to make it super wide so as to not influence Wermon
  }
  
  
  # Model-specific parameters
  whichPop <- 1:n.timeseries         # multiple time series -> unique population processes
  n.states <- max(whichPop)
  
  Z <- matrix(0,n.timeseries+1,n.states+1)   # matrix with rows as n.timeseries and cols as n.states (pops)
  Z[n.timeseries+1, ] <- NA                  # add a row of NAs to keep jagsUI from converting single time series matrix into vector
  Z[ , n.states+1] <- NA                     # add a col of NAs to keep jagsUI from converting single state matrix into vector
  for(i in 1:length(whichPop)) Z[i,whichPop[i]] <- 1
  
  jags.data <- list(Y = rbind(dat, NA),
                    n.yrs = ncol(dat),
                    n.timeseries = n.timeseries,
                    n.states = n.states,
                    Z = Z,
                    u_mean = 0,
                    u_sd = 0.5,
                    q_alpha = 0.01,
                    q_beta = 0.01,
                    r_alpha = 0.01,
                    r_beta = 0.01,
                    x0_mean = x0_mean,
                    x0_sd = x0_sd)
  
  jags.model <- jags(jags.data, 
                     inits = NULL, 
                     parameters.to.save= jags.params, 
                     model.file=model.loc, 
                     n.chains = MCMC.params$n.chains, 
                     n.burnin = MCMC.params$n.burnin, 
                     n.thin = MCMC.params$n.thin, 
                     n.iter = MCMC.params$n.samples, 
                     DIC = T, 
                     parallel=T, 
                     seed=set.seed)
  
  jm.out <- list(jags.data = jags.data,
                 params = jags.params,
                 out = jags.model)
  
  return(jm.out)
}

run.independentQs_singleU <- function(dat, 
                                      jags.params, 
                                      model.loc, 
                                      MCMC.params, 
                                      seed){
  n.timeseries <- nrow(dat)
  
  # Model-specific parameters
  whichPop <- 1:n.timeseries         # multiple time series -> unique population processes
  n.states <- max(whichPop)
  
  Z <- matrix(0,n.timeseries+1, n.states+1)   # matrix with rows as n.timeseries and cols as n.states (pops)
  Z[n.timeseries+1, ] <- NA                  # add a row of NAs to keep jagsUI from converting single time series matrix into vector
  Z[ , n.states+1] <- NA                     # add a col of NAs to keep jagsUI from converting single state matrix into vector
  for(i in 1:length(whichPop)) Z[i,whichPop[i]] <- 1
  
  jags.data <- list(Y = rbind(dat, NA),
                    n.yrs = ncol(dat),
                    n.timeseries = n.timeseries,
                    n.states = n.states,
                    Z = Z,
                    u_mean = 0,
                    u_sd = 0.5,
                    q_alpha = 0.01,
                    q_beta = 0.01,
                    r_alpha = 0.01,
                    r_beta = 0.01,
                    x0_mean = dat[1,1],
                    x0_sd = 10)
  
  jags.model <- jags(jags.data, 
                     inits = NULL, 
                     parameters.to.save= jags.params, 
                     model.file=model.loc, 
                     n.chains = MCMC.params$n.chains, 
                     n.burnin = MCMC.params$n.burnin, 
                     n.thin = MCMC.params$n.thin, 
                     n.iter = MCMC.params$n.samples, 
                     DIC = T, 
                     parallel=T, 
                     seed=set.seed)
  
  jm.out <- list(jags.data = jags.data,
                 params = jags.params,
                 out = jags.model)
  
}
  
  
run.independentUs_singleQ <- function(dat, 
                                      jags.params, 
                                      model.loc, 
                                      MCMC.params, 
                                      seed){
  
  n.timeseries <- nrow(dat)
  x0_mean <- numeric(length=n.timeseries)
  x0_sd <- numeric(length=n.timeseries)
  for(i in 1:n.timeseries) {
    icol <- min(which(!is.na(dat[i,])))   # use the first non-NA data point in time series i as it's prior mean
    x0_mean[i] <- dat[i, icol]            
    x0_sd[i] <- 10     # we went with wide sd by testing for both JM and W for leathers; 10 works to make it super wide so as to not influence Wermon
  }
  
  
  # Model-specific parameters
  whichPop <- 1:n.timeseries         # multiple time series -> unique population processes
  n.states <- max(whichPop)
  
  Z <- matrix(0,n.timeseries+1,n.states+1)   # matrix with rows as n.timeseries and cols as n.states (pops)
  Z[n.timeseries+1, ] <- NA                  # add a row of NAs to keep jagsUI from converting single time series matrix into vector
  Z[ , n.states+1] <- NA                     # add a col of NAs to keep jagsUI from converting single state matrix into vector
  for(i in 1:length(whichPop)) Z[i,whichPop[i]] <- 1
  
  jags.data <- list(Y = rbind(dat, NA),
                    n.yrs = ncol(dat),
                    n.timeseries = n.timeseries,
                    n.states = n.states,
                    Z = Z,
                    u_mean = 0,
                    u_sd = 0.5,
                    q_alpha = 0.01,
                    q_beta = 0.01,
                    r_alpha = 0.01,
                    r_beta = 0.01,
                    x0_mean = x0_mean,
                    x0_sd = x0_sd)
  
  jags.model <- jags(jags.data, 
                     inits = NULL, 
                     parameters.to.save= jags.params, 
                     model.file=model.loc, 
                     n.chains = MCMC.params$n.chains, 
                     n.burnin = MCMC.params$n.burnin, 
                     n.thin = MCMC.params$n.thin, 
                     n.iter = MCMC.params$n.samples, 
                     DIC = T, 
                     parallel=T, 
                     seed=set.seed)
  
  jm.out <- list(jags.data = jags.data,
                 params = jags.params,
                 out = jags.model)
  
  return(jm.out)
}
