
rm(list=ls())

source("Dc_Indonesia_nesting_fcns.R")
library(jagsUI)
library(coda)
# library(ggplot2)
library(loo)


MCMC.n.chains <- 5
MCMC.n.samples <- 100000
MCMC.n.burnin <- 50000
MCMC.n.thin <- 5

MCMC.params <- list(n.chains = MCMC.n.chains,
                    n.samples = MCMC.n.samples,
                    n.burnin = MCMC.n.burnin,
                    n.thin = MCMC.n.thin)

year.begin <- 2001
year.end <- 2023
#year.end <- 2019
season.begin <- 2001
season.end <- 2023
loc <- "JM"
data.jags <- data.extract(location = loc, 
                          year.begin = year.begin, 
                          year.end = year.end,
                          season.begin = season.begin,
                          season.end = season.end)

norm.norm.models <- c("models/model_SSAR1_logY_norm_norm_var_thetaM.txt")

#norm.norm.models <- c("models/model_SSAR1_logY_norm_norm.txt",
#                       "models/model_SSAR1_logY_norm_norm_theta.txt",
#                       "models/model_SSAR1_logY_norm_norm_var.txt",
#                       "models/model_SSAR1_logY_norm_norm_var_theta.txt",
#                       "models/model_SSAR1_logY_norm_norm_varM_thetaM.txt",
#                       "models/model_SSAR1_logY_norm_norm_varM_theta.txt",
#                       "models/model_SSAR1_logY_norm_norm_var_thetaM.txt",
#                       "models/model_SSAR1_logY_norm_norm_thetaM.txt",
#                       "models/model_SSAR1_logY_norm_norm_varM.txt")

norm.t.models <- c()
# norm.t.models <- c("models/model_SSAR1_logY_norm_t.txt",
#                    "models/model_SSAR1_logY_norm_t_theta.txt",
#                    "models/model_SSAR1_logY_norm_t_var.txt",
#                    "models/model_SSAR1_logY_norm_t_var_theta.txt",
#                    "models/model_SSAR1_logY_norm_t_varM_thetaM.txt",
#                    "models/model_SSAR1_logY_norm_t_varM_theta.txt",
#                    "models/model_SSAR1_logY_norm_t_var_thetaM.txt",
#                    "models/model_SSAR1_logY_norm_t_thetaM.txt",
#                    "models/model_SSAR1_logY_norm_t_varM.txt")

all.models <- c(norm.norm.models, norm.t.models)

params.all <- c("theta.1", 'sigma.pro1', "sigma.obs",
                "mu", "y", "X", "deviance", "loglik")

for (k in 1:length(all.models)){
  print(paste("file", k, "of", length(all.models), "models"))
  print(all.models[[k]])
  
  tic <- Sys.time()
  if (k > length(norm.norm.models)){
    m <- k - length(norm.norm.models)
    params.all <- c(params.all, "df")  # need df for norm-t models
  } else {
    m <- k
  }
  
  if (m == 1 | m == 5 | m == 8 | m == 9){ 
    jags.params <- params.all 
  } else if (m == 2 | m == 6){
    jags.params <- c(params.all, "theta.2")
  } else if (m == 3 | m == 7) {
    jags.params <- c(params.all, "sigma.pro2")
  } else if (m == 4){
    jags.params <- c(params.all, "theta.2", "sigma.pro2") 
  } 
  
  jags.out <- run.jagsUI(data.jags$jags.data, 
                         jags.params, 
                         model.file = all.models[[k]], 
                         MCMC.params)
  
  Xs.stats <- jags.out$Xs.stats
  
  Xs.stats$time <- data.jags$data.1$Frac.Year
  Xs.stats$obsY <- data.jags$data.1$Nests
  Xs.stats$month <- data.jags$data.1$Month
  Xs.stats$year <- data.jags$data.1$Year
  
  ys.stats <- jags.out$ys.stats
  ys.stats$time <- data.jags$data.1$Frac.Year
  ys.stats$obsY <- data.jags$data.1$Nests
  ys.stats$month <- data.jags$data.1$Month
  ys.stats$year <- data.jags$data.1$Year
  
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
  
  filename.root <- strsplit(strsplit(all.models[[k]], 
                                     'models/model_')[[1]][2], '.txt')[[1]][1]
  
  if (length(strsplit(filename.root, paste0("SSAR1_", loc, "_"))[[1]]) == 1){
    
    filename.out <- paste0(filename.root, "_", loc, "_", year.begin, "_", year.end)
    
  } else {
    filename.root <- strsplit(filename.root, paste0("SSAR1_", loc, "_"))[[1]][2]
    filename.out <- paste0("SSAR1_", filename.root, "_", loc, "_", year.begin, "_", year.end)
  }
  
  toc <- Sys.time()
  
  results <- list(data.1 = data.jags$data.1,
                  jags.out = jags.out,
                  Xs.stats = Xs.stats,
                  ys.stats = ys.stats,
                  MCMC.params = MCMC.params,
                  time = toc - tic)
  
  ggsave(plot = p.1,
         filename = paste0("figures/", "predicted_counts_", filename.out, ".png"),
         dpi = 600)
  
  saveRDS(results,
          file = paste0('RData/', "jagsout_", 
                        filename.out, "_", Sys.Date(), '.rds'))
  
}

