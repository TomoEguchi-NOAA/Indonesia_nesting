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

