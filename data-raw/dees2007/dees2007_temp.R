rm(list = ls())

library(zoo)

setwd("/home/franz/Dropbox/00_Current projects/bgvars/data-raw/dees2007")

load("dees26.RData")

country.data <- c()
for (i in 1:length(data)) {
  i<- 1
  
  zoo(data[[i]], order.by = index(data[[i]]), frequency = 4)
}
