#Dc_Indonesia_nesting_fcns

# these are functions that are used in this project.
# TomosFunctions.R is needed also.

#sysInfo <- Sys.info()
# ifelse(Sys.info()[1] == 'Linux',
#        source('~/Documents/R/tools/TomosFunctions.R'),
#        source('~/R/tools/TomosFunctions.R'))

library(ggplot2)
library(tidyverse)
library(lubridate)
library(loo)
library(dplyr)

save.fig <- F

# 
# # Extracting posterior samples of deviance or any other variable from jags output:
# extract.samples <- function(varname, zm){
#   dev <- unlist(lapply(zm, FUN = function(x) x[, varname]))
#   return(dev)
# }

compute.LOOIC <- function(loglik, data.vector, MCMC.params){
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  loglik.vec <- as.vector(loglik)
  
  # each column corresponds to a data point and rows are MCMC samples
  loglik.mat <- matrix(loglik.vec, nrow = n.per.chain * MCMC.params$n.chains)
  
  # take out the columns that correspond to missing data points
  loglik.mat <- loglik.mat[, !is.na(data.vector)]
  # loglik.mat <- matrix(loglik.vec[!is.na(data.vector)], 
  #                      nrow = MCMC.params$n.chains * n.per.chain)
  
  Reff <- relative_eff(exp(loglik.mat),
                       chain_id = rep(1:MCMC.params$n.chains,
                                      each = n.per.chain),
                       cores = 4)
  
  #
  loo.out <- loo(loglik.mat, 
                 r_eff = Reff, 
                 cores = 4, k_threshold = 0.7)
  
  out.list <- list(Reff = Reff,
                   loo.out = loo.out)
  
  return(out.list)  
}

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

data.extract <- function(location, year.begin, year.end, season.begin = year.begin, season.end = year.end){
  # In March 2019, we received new data for 2018. So, the raw data file
  # has been updated.  
  # On 16 April 2019, the last few data points for 2019 were received
  # so the data files have been updated. 
  if (is.null(season.begin)) season.begin <- year.begin
  if (is.null(season.end)) season.end <- year.end
  
  if (location == "JM"){
    data.0 <- read.csv("data/JM_nests_Sept2020.csv")
    
    data.0 %>% 
      select(Year_begin, Month_begin, JM_Nests) %>%
      mutate(Nests = JM_Nests) -> data.0
    
  } else if (location == "W"){
    data.0 <- read.csv("data/W_nests_Sept2020.csv")
    data.0 %>% 
      select(Year_begin, Month_begin, W_Nests) %>%
      mutate(Nests = W_Nests) -> data.0
  }
  
  # create regularly spaced time series:
  data.frame(Year = rep(year.begin:max(data.0$Year_begin,
                                       na.rm = T),
                        each = 12),
             Month_begin = rep(1:12,
                               max(data.0$Year_begin,
                                   na.rm = T) -
                                 year.begin + 1)) %>%
    mutate(begin_date = as.Date(paste(Year,
                                      Month_begin,
                                      '01', sep = "-"),
                                format = "%Y-%m-%d"),
           Frac.Year = Year + (Month_begin-0.5)/12) %>%
    select(Year, Month_begin, begin_date, Frac.Year) -> data.2

  # also make "nesting season" that starts April and ends March
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
              Nests = Nests,
              Season = ifelse(Month < 4, Year-1, Year),
              Seq.Month = ifelse(Month < 4, Month + 9, Month - 3)) %>%
    reshape::sort_df(.,vars = "Frac.Year") %>%
    filter(Season >= season.begin & Season <= season.end) -> data.1
  
  data.1 %>% filter(Month > 3 & Month < 10) -> data.summer
  data.1 %>% filter(Month > 9 | Month < 4) %>%
    mutate(Seq.Month = Seq.Month - 6) -> data.winter
  
  jags.data <- list(y = log(data.1$Nests),
                    m = data.1$Seq.Month,
                    T = nrow(data.1))
  
  y <- matrix(log(data.1$Nests), ncol = 12, byrow = TRUE)
  
  jags.data2 <- list(y = y,
                     m = matrix(data.1$Seq.Month, 
                                ncol = 12, byrow = TRUE),
                     n.years = nrow(y))
  
  y.summer <- matrix(log(data.summer$Nests),
                     ncol = 6, byrow = TRUE)
  
  y.winter <- matrix(log(data.winter$Nests),
                     ncol = 6, byrow = TRUE)
  
  jags.data2.summer <- list(y = y.summer,
                            m = matrix(data.summer$Seq.Month, 
                                       ncol = 6, byrow = TRUE),
                            n.years = nrow(y.summer))
  
  jags.data2.winter <- list(y = y.winter,
                            m = matrix(data.winter$Seq.Month, 
                                       ncol = 6, byrow = TRUE),
                            n.years = nrow(y.winter))
  
  out <- list(jags.data = jags.data,
              jags.data2 = jags.data2,
              jags.data.summer = jags.data2.summer,
              jags.data.winter = jags.data2.winter,
              data.1 = data.1,
              data.summer = data.summer,
              data.winter = data.winter)
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

# extracts loo statistics, including looic and Pareto-k statistics
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

# same as the previous one but when data are in an array, not a vector
# MCMC output is a 3D array, so that needs to be converted (flattened)
pareto.k.diag.3D <- function(jm, MCMC.params, jags.data){
  
  n.per.chain <- (MCMC.params$n.samples - MCMC.params$n.burnin)/MCMC.params$n.thin
  
  # reduce the dimension by MCMC iterations by year x month
  loglik <- t(apply(jm$sims.list$loglik, 
                    MARGIN = 1, 
                    FUN = function(x) as.vector(t(x))))
  
  # convert the data (y) into a vector also:
  y <- as.vector(t(jags.data$y))
  
  loglik.obs <- loglik[, !is.na(y)]
  
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
              loo.out = loo.out,
              y = y))
}

model.names <- function(){
  # norm.norm.models <- c("model_SSAR1_logY_norm_norm.txt",
  #                       "model_SSAR1_logY_norm_norm_theta.txt",
  #                       "model_SSAR1_logY_norm_norm_var.txt",
  #                       "model_SSAR1_logY_norm_norm_var_theta.txt",
  #                       "model_SSAR1_logY_norm_norm_varM_thetaM.txt",
  #                       "model_SSAR1_logY_norm_norm_varM_theta.txt",
  #                       "model_SSAR1_logY_norm_norm_var_thetaM.txt",
  #                       "model_SSAR1_logY_norm_norm_thetaM.txt",
  #                       "model_SSAR1_logY_norm_norm_varM.txt")
  # 
  # norm.t.models <- c("model_SSAR1_logY_norm_t.txt",
  #                    "model_SSAR1_logY_norm_t_theta.txt",
  #                    "model_SSAR1_logY_norm_t_var.txt",
  #                    "model_SSAR1_logY_norm_t_var_theta.txt",
  #                    "model_SSAR1_logY_norm_t_varM_thetaM.txt",
  #                    "model_SSAR1_logY_norm_t_varM_theta.txt",
  #                    "model_SSAR1_logY_norm_t_var_thetaM.txt",
  #                    "model_SSAR1_logY_norm_t_thetaM.txt",
  #                    "model_SSAR1_logY_norm_t_varM.txt")
  
  Fourier.models <- c("models/model_SSAR1_logY_norm_norm_theta_Four.txt",
                      "models/model_SSAR1_logY_norm_norm_theta_Four_constCV.txt",
                      "models/model_SSAR1_logY_norm_norm_theta_Four_varM.txt",
                      "models/model_SSAR1_logY_norm_t_theta_Four.txt",
                      "models/model_SSAR1_logY_norm_t_theta_Four_constCV.txt",
                      "models/model_SSAR1_logY_norm_t_theta_Four_varM.txt")
  
  #all.models <- c(norm.norm.models, norm.t.models, Fourier.models)
  
  models <- data.frame(number = paste0("M", 1:length(Fourier.models)),
                       names = Fourier.models)
  
  return(models)
  
}

run.singleUQ <- function(dat, 
                         jags.params, 
                         model.loc, 
                         MCMC.params, 
                         seed){
  
  n.timeseries <- nrow(dat)
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


model.Comparison.Fourier <- function(loc, year.begin, year.end, run.date){
  models <- model.names()
        
  M.ys.stats <- M.Xs.stats <- M.jm <- M.loo.out <- M.data.y <- vector(mode = "list", 
                                                          length = nrow(models))
  
  for (k in 1:nrow(models)){
    str.root <- strsplit(strsplit(as.character(models[k, "names"]), "model_")[[1]][2], ".txt")[[1]]
    M <- readRDS(paste0("RData/jagsout_", str.root, "_", loc, "_", 
                        year.begin, "_", year.end, "_",
                        run.date, ".rds"))
    
    M.loo.out[[k]] <- M$loo.out$loo.out
    M.jm[[k]] <- M$jags.out
    M.Xs.stats[[k]] <- M$Xs.stats
    M.ys.stats[[k]] <- M$ys.stats
    M.data.y[[k]] <- na.omit(M$jags.data$y)
    
  }           
  
  out <- list(loo.out = M.loo.out,
              jm = M.jm,
              Xs.stats = M.Xs.stats,
              ys.stats = M.ys.stats,
              data.y = M.data.y,
              models = models)
  
  return(out)
}

PSIS.plots <- function(pareto.k, data.y, m, loc, save.fig = TRUE){
  pareto.df <- data.frame(y = data.y,
                          khat = pareto.k,
                          datapoint = seq(from = 1, to = length(data.y)),
                          k0.7 = cut(pareto.k,
                                     breaks = c(0, 0.7, 1.5),
                                     labels = c("<=0.7", ">0.7")))
  
  p.1 <- ggplot(data = pareto.df) +   
    geom_path(aes(x = datapoint, y = exp(y)), alpha = 0.5) +
    geom_point(aes(x = datapoint, y = exp(y), 
                   size = khat,
                   color = k0.7)) +
    scale_size_continuous(limits = c(0.0, 1.3),
                          range = c(1, 4))+ 
    scale_color_manual(values = c("<=0.7" = "black", 
                                  ">0.7" = "red"))
  
  p.2 <- ggplot(data = pareto.df) +  
    geom_point(aes(x = datapoint, 
                   y = khat,
                   color = -khat)) + 
    geom_hline(yintercept = 0.5,
               color = "red",
               size = 0.5) +
    geom_hline(yintercept = 0.7,
               color = "red",
               size = 0.7) +
    geom_hline(yintercept = 1.0,
               color = "red",
               size = 1.0) +
    theme(legend.position = "none") + 
    labs(x = "Data point", 
         y = "Pareto shape k",
         title = paste0("PSIS diagnostic plot (Model ", m, ")"))
  
  if (save.fig)
    ggsave(p.2, 
           filename = paste0("figures/", loc, "_PSIS_M", m, ".png"),
           device = "png", dpi = 600)
  
  return(list(p1 = p.1, p2 = p.2))
  
}

plot.combined.counts <- function(M.jm, Xs.stats, ys.stats, 
                                 year.begin, year.end, run.date, loc, m,
                                 fill.color =  "darkseagreen",
                                 fill.color.summer = "darksalmon",
                                 fill.color.winter = "gray65",
                                 fill.alpha =  0.65,
                                 line.color =  "darkblue",
                                 line.color.summer = "red3",
                                 line.color.winter ="red3",  #"greenyellow"
                                 data.color = "black",
                                 data.size = 1.5,
                                 obsd.color = "red2"){
  
  Xs.stats %>% mutate(season = ifelse(month < 4, year - 1, year)) %>% 
    mutate(summer = ifelse(month > 3 & month < 10, 1, 0)) -> Xs.stats
  
  ys.stats %>% mutate(season = ifelse(month < 4, year - 1, year)) %>% 
    mutate(summer = ifelse(month > 3 & month < 10, 1, 0)) -> ys.stats
  
  seasons <- as.matrix(unique(Xs.stats$season))
  
  #######################
  X.posterior.seasons <- lapply(apply(seasons, 
                                      MARGIN = 1,
                                      FUN = extract.posterior.jagsUI, 
                                      Xs.stats = Xs.stats, 
                                      samples = M.jm$samples),
                                FUN = function(x){
                                  n.summer <- exp(x$samples[, x$summer == 1]) %>% rowSums() 
                                  n.winter <- exp(x$samples[, x$summer == 0]) %>% rowSums() 
                                  n.season <- exp(x$samples) %>% rowSums()
                                  return(data.frame(summer = n.summer, 
                                                    winter = n.winter,
                                                    all = n.season))
                                } )
  
  
  Xs.season <- as.data.frame(matrix(unlist(lapply(X.posterior.seasons, 
                                                  FUN = function(x){
                                                    qtiles <- apply(x, 
                                                                    MARGIN = 2,
                                                                    FUN = quantile, 
                                                                    probs = c(0.025, 0.5, 0.975))})),
                                    ncol = 9, byrow = T))
  
  colnames(Xs.season) <- c("Summer.low", "Summer.median", "Summer.high",
                           "Winter.low", "Winter.median", "Winter.high",
                           "all.low", "all.median", "all.high")
  Xs.season <- mutate(Xs.season, 
                      season = as.vector(seasons))
  
  ys.stats %>% group_by(season, summer) %>%
    summarize(obs = sum(obsY)) -> tmp.y.season
  
  ys.season <- data.frame(season = seasons,
                          summer = unlist(c(as.vector(filter(tmp.y.season, summer == 1)[, "obs"]))),
                          winter = unlist(filter(tmp.y.season, summer == 0)[, "obs"]))
  
  Xs.season %>% left_join(ys.season, by = "season") -> Xy.season
  
  if (!file.exists(paste0("data/estimatedX_", loc, "_M", 
                          m, "_", run.date, ".csv")))
    write.csv(Xy.season, 
              file = paste0("data/estimatedX_", loc, "_M", 
                            m, "_", run.date, ".csv"),
              quote = F, row.names = F)
  
  p.estimated.counts <- ggplot() + 
    geom_ribbon(data = Xs.season,
                aes(x = season, 
                    ymin = Summer.low, 
                    ymax = Summer.high),
                fill = fill.color.summer,
                alpha = fill.alpha) +
    geom_point(data = Xs.season,
               aes(x = season, 
                   y = Summer.median), 
               color = line.color.summer,
               size = 2) +
    geom_point(data = ys.season,
               aes(x = season, 
                   y = summer),
               color = data.color,
               size = 2) +
    
    geom_ribbon(data = Xs.season,
                aes(x = season, 
                    ymin = Winter.low, 
                    ymax = Winter.high), 
                fill = fill.color.winter,
                alpha = fill.alpha) + 
    geom_point(data = Xs.season,
               aes(x = season, 
                   y = Winter.median), 
               color = line.color.winter,
               shape = 17,
               size = 2)+
    
    geom_point(data = ys.season,
               aes(x = season, 
                   y = winter),
               color = data.color,
               shape = 17,
               size = 2) + 
    scale_x_continuous(breaks = seq(year.begin, year.end, 5),
                       limits = c(year.begin, year.end)) +
    labs(x = '', y = '# nests', 
         title = ifelse(loc == "JM", "Jamursba-Medi", "Wermon"))  +
    theme(axis.text = element_text(size = 12),
          text = element_text(size = 12))
  
  return(p.estimated.counts)
}

plot.imputed <- function(Xs.stats, ys.stats, loc, year.begin, year.end, 
                         fill.color =  "darkseagreen",
                         fill.alpha =  0.65,
                         line.color =  "darkblue",
                         data.color = "black",
                         data.size = 1.5,
                         obsd.color = "red2"){
  
  p.1 <- ggplot() +
    geom_ribbon(data = Xs.stats,
                aes(x = time, 
                    ymin = exp(low_X), 
                    ymax = exp(high_X)),
                fill = fill.color,
                alpha = fill.alpha) +
    #geom_point(data = ys.stats,
    #           aes(x = time, y = mode_y), color = "blue") +
    #geom_line(data = Xs.stats,
    #          aes(x = time, y = mode_X), color = 'blue') +
    # geom_line(data = Xs.stats,
    #           aes(x = time, y = exp(high_X)), 
    #           color = "purple",
    #           linetype = 2) +
    geom_point(data = Xs.stats,
               aes(x = time, y = exp(median_X)), 
               color = line.color,
               size = 1) +
    geom_line(data = Xs.stats,
              aes(x = time, y = exp(median_X)), 
              color = line.color,
              alpha = 0.5) +
    geom_point(data = ys.stats,
               aes(x = time, y = obsY), 
               color = obsd.color,
               alpha = 0.5) + 
    scale_x_continuous(breaks = seq(year.begin, year.end, 5),
                       limits = c(year.begin, year.end)) +
    labs(x = '', y = '# nests', 
         title = ifelse(loc == "JM", "Jamursba-Medi", "Wermon"))  +
    theme(axis.text = element_text(size = 12),
          text = element_text(size = 12))
  
  p.1.log <- ggplot() +
    geom_ribbon(data = Xs.stats,
                aes(x = time, 
                    ymin = (low_X), 
                    ymax = (high_X)),
                fill = fill.color,
                alpha = fill.alpha) +
    #geom_point(data = ys.stats,
    #           aes(x = time, y = mode_y), color = "blue") +
    #geom_line(data = Xs.stats,
    #          aes(x = time, y = mode_X), color = 'blue') +
    # geom_line(data = Xs.stats,
    #           aes(x = time, y = exp(high_X)), 
    #           color = "purple",
    #           linetype = 2) +
    geom_point(data = Xs.stats,
               aes(x = time, y = (median_X)), 
               color = line.color,
               size = 1) +
    geom_line(data = Xs.stats,
              aes(x = time, y = (median_X)), 
              color = line.color,
              alpha = 0.5) +
    geom_point(data = ys.stats,
               aes(x = time, y = log(obsY)), 
               color = obsd.color,
               alpha = 0.5) + 
    scale_x_continuous(breaks = seq(year.begin, year.end, 5),
                       limits = c(year.begin, year.end)) +
    labs(x = '', y = 'log(# nests)', 
         title = ifelse(loc == "JM", "Jamursba-Medi", "Wermon"))  +
    theme(axis.text = element_text(size = 12),
          text = element_text(size = 12))
  
  ps <- list(p.1 = p.1,
             p.1.log = p.1.log)
  
}
