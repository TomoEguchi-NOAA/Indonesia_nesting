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

save.fig <- F

sum.posterior <- function(yr = 2013, Xs.stats, zm) {
  Xs.stats %>% mutate(var.name = rownames(Xs.stats)) %>%
    filter(year == yr) %>%
    select(var.name) -> Xs.name
  zm.yr <- apply(Xs.name, MARGIN = 1,
                 FUN = extract.samples,
                 zm)

  return(list(samples = zm.yr, var.names = Xs.name))
}

