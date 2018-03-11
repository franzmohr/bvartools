rm(list=ls())

setwd("~/Dropbox/00_master/computations") # Set working directory

load("~/Dropbox/00_master/computations/data/pesaran26.rda") # Load GVAR toolbox (2013) data
load("~/Dropbox/00_master/computations/data/gvar_db/gvar_db26.RData") # Load GVAR toolbox (2013) data

c<-"Malaysia";par(mfcol=c(3,2));for (j in dimnames(Data[[which(names(Data)==c)]])[[2]]){
  plot.ts(cbind(Data[[which(names(Data)==c)]][,j],data[[which(names(data)==c)]][,j]),plot.type = "single",col=c("black","blue"),ylab=j)
};par(mfcol=c(1,1))

cbind(colSums(weight)[order(names(colSums(weight)))],colSums(weights)[order(names(colSums(weights)))])
cbind(rowSums(weight)[order(names(rowSums(weight)))],rowSums(weights)[order(names(rowSums(weights)))])
