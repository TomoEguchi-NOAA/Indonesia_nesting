

rm(list=ls())

source("Dc_Indonesia_nesting_fcns.R")
library(tidyverse)
library(readxl)
library(jagsUI)
library(coda)
# library(ggplot2)
library(loo)

run.date <- "2023-07-11" #  "2019-04-30" #Sys.Date() #"2019-04-29" #

MCMC.params <- list(n.chains = 5,
                    n.samples = 500000,
                    n.burnin = 350000,
                    n.thin = 50)

loc <- "JM"
if (loc == "JM"){
  year.begin <- 1999
  season.begin <- 1999
  period <- 12   # should be 12 for Jamusrba-Medi
  year.end <- 2023
  season.end <- 2022
  
} else if (loc == "W"){
  year.begin <- 2003  
  season.begin <- 2003
  
  # for model comparison purpose
  #year.begin <- 2006  # for model comparison purpose
  #season.begin <- 2006
  period <- 6   # should be 6 for Wermon 
  
  year.end <- 2023
  season.end <- 2022
}

# year.end <- 2019
# season.end <- 2018

data.jags <- data.extract(location = loc, 
                          year.begin = year.begin, 
                          year.end = year.end,
                          season.begin = season.begin,
                          season.end = season.end)

models.df <- model.names()

# norm.norm.models <- c("models/model_SSAR1_logY_norm_norm_theta_Four.txt",
#                       "models/model_SSAR1_logY_norm_norm_theta_Four_constCV.txt",
#                       "models/model_SSAR1_logY_norm_norm_theta_Four_varM.txt", #)
#                       "models/model_SSAR1_logY_norm_norm_varM_thetaM.txt",
#                       "models/model_SSAR1_W_logY_norm_norm_varM_theta.txt",
#                       "models/model_SSAR1_W_logY_norm_norm_var_thetaM.txt",
#                       "models/model_SSAR1_logY_norm_norm_thetaM.txt",
#                       "models/model_SSAR1_logY_norm_norm_varM.txt")
# # 
# norm.t.models <- c("models/model_SSAR1_logY_norm_t_theta_Four.txt",
#                    "models/model_SSAR1_logY_norm_t_theta_Four_constCV.txt",
#                    "models/model_SSAR1_logY_norm_t_theta_Four_varM.txt", #)
#                    "models/model_SSAR1_W_logY_norm_t_var.txt",
#                    "models/model_SSAR1_W_logY_norm_t_var_theta.txt",
#                    "models/model_SSAR1_logY_norm_t_varM_thetaM.txt",
#                    "models/model_SSAR1_W_logY_norm_t_varM_theta.txt",
#                    "models/model_SSAR1_W_logY_norm_t_var_thetaM.txt",
#                    "models/model_SSAR1_logY_norm_t_thetaM.txt",
#                    "models/model_SSAR1_logY_norm_t_varM.txt")
# # 
# all.models <- c(norm.norm.models, norm.t.models)

jags.params <- c("theta", 'sigma.pro1', "sigma.obs",
                "mu", "y", "X", "theta.beta.sin",
                "theta.beta.cos", "deviance", "loglik")

jags.data <- list(y = data.jags$jags.data$y,
                  T = data.jags$jags.data$T,
                  m = data.jags$jags.data$m,
                  pi = pi,
                  period = period)

k <- 2
for (k in 1:nrow(models.df)){
  print(paste("file", k, "of", nrow(models.df), "models"))
  print(models.df[k, "names"])
  model.name <- models.df[k, "names"]
  
  filename.root <- strsplit(strsplit(models.df[k, "names"], 
                                     'models/model_')[[1]][2], '.txt')[[1]][1]
  
  filename.out <- paste0(filename.root, "_", loc, "_", year.begin, "_", year.end)
  
  if (!file.exists(  paste0('RData/', "jagsout_", 
                            filename.out, "_", run.date, '.rds'))){
    tic <- Sys.time()
    
    if (length(grep("_t_", models.df[k, "names"])) > 0){
      jags.params <- c(jags.params, "df")  # need df for norm-t models
    }
    
    if (length(grep("_constCV", models.df[k, "names"])) > 0){
      jags.params <- c(jags.params, "cv.pro1")
    }
    
    jm <- jags(jags.data,
               inits = NULL,
               parameters.to.save= jags.params,
               model.file = model.name,
               n.chains = MCMC.params$n.chains,
               n.burnin = MCMC.params$n.burnin,
               n.thin = MCMC.params$n.thin,
               n.iter = MCMC.params$n.samples,
               DIC = T, parallel=T)
    
    ys.stats <- data.frame(low_y = jm$q2.5$y,
                           median_y = jm$q50$y,
                           high_y = jm$q97.5$y,
                           obsY = data.jags$data.1$Nests,
                           month = data.jags$data.1$Month,
                           year = data.jags$data.1$Year,
                           Season = data.jags$data.1$Season,
                           time = data.jags$data.1$Frac.Year)
    
    # extract Xs - the state model
    Xs.stats <- data.frame(low_X = jm$q2.5$X,
                           median_X = jm$q50$X,
                           high_X = jm$q97.5$X,
                           obsY = data.jags$data.1$Nests,
                           month = data.jags$data.1$Month,
                           year = data.jags$data.1$Year,
                           Season = data.jags$data.1$Season,
                           time = data.jags$data.1$Frac.Year)
    
    loo.out <- pareto.k.diag(jm, MCMC.params, jags.data)
    
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
    
    p.2 <- ggplot() +
      #geom_point(data = ys.stats,
      #           aes(x = time, y = mode_y), color = "blue") +
      #geom_line(data = Xs.stats,
      #          aes(x = time, y = mode_X), color = 'blue') +
      geom_line(data = Xs.stats,
                aes(x = time, y = (high_X)), color = "red",
                linetype = 2) +
      geom_point(data = Xs.stats,
                 aes(x = time, y = (median_X)), color = "red",
                 alpha = 0.5) +
      geom_line(data = Xs.stats,
                aes(x = time, y = (median_X)), color = "red",
                alpha = 0.5) +
      geom_line(data = Xs.stats,
                aes(x = time, y = (low_X)), color = "red",
                linetype = 2) +
      geom_point(data = ys.stats,
                 aes(x = time, y = log(obsY)), color = "green",
                 alpha = 0.5)
    
    
    #if (length(strsplit(filename.root, paste0("SSAR1_", loc, "_"))[[1]]) == 1){
    
    # } else {
    #   filename.root <- strsplit(filename.root, paste0("SSAR1_", loc, "_"))[[1]][2]
    #   filename.out <- paste0("SSAR1_", filename.root, "_", loc, "_", year.begin, "_", year.end)
    # }
    
    toc <- Sys.time()
    results <- list(data.1 = data.jags$data.1,
                    jags.data = jags.data,
                    jags.out = jm,
                    Xs.stats = Xs.stats,
                    ys.stats = ys.stats,
                    MCMC.params = MCMC.params,
                    loo.out = loo.out,
                    time = toc - tic,
                    system = Sys.getenv())
    
    ggsave(plot = p.1,
           filename = paste0("figures/", "predicted_counts_", filename.out, ".png"),
           dpi = 600)
    
    ggsave(plot = p.2,
           filename = paste0("figures/", "predicted_log_counts_", filename.out, ".png"),
           dpi = 600)
    
    saveRDS(results,
            file = paste0('RData/', "jagsout_", 
                          filename.out, "_", run.date, '.rds'))
    
  }

  

  
}

