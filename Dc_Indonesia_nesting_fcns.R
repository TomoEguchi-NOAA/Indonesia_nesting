#Dc_Indonesia_nesting_fcns


sysInfo <- Sys.info()
ifelse(sysInfo[1] == 'Linux',
       source('~/Documents/R/TomosFunctions.R'),
       source('~/R/TomosFunctions.R'))

library(ggplot2)
library(tidyverse)
library(lubridate)

mmm2month <- function(x){
  switch(as.character(x),
         "Jan" = 1, "Feb" = 2, "Mar" = 3, "Apr" = 4,
         "May" = 5, "Jun" = 6, "Jul" = 7, "Aug" = 8,
         "Sep" = 9, "Oct" = 10, "Nov" = 11, "Dec" = 12, NA)
}
