rm(list=ls())
setwd("~/Dropbox/00_master/computations/data/gvar2013") # install.packages("xlsx")
#setwd("E:/Franz/Dropbox/00_master/computations/data/gvar2013") # install.packages("xlsx")
library(xlsx)

countries <- read.xlsx("country_codes.xls",1,startRow=1,endRow=34)
countries[,1] <- as.character(countries[,1])
countries <- countries[order(countries[,1]),]

vars <- c("RGDP", "CPI", "nomEQ", "FXdol", "Rshort", "Rlong")
global <- c("POILdol", "PRAWMATdol", "PMETALdol")

data <- list()
for (i in 1:33){
  temp <- ts(read.xlsx("country_data.xlsx",sheetIndex=i,startRow=1,endRow=138)[,-1],
            start=1979,frequency = 4)
  pos <- which(apply(temp,2,function(x){all(x==x[1])}))
    if (length(pos)>0) {
      temp <- temp[,-pos]
    }
  temp <- temp[,which(is.element(dimnames(temp)[[2]],vars))]
  data <- c(data,list(temp))
  names(data)[i] <- countries[i,1]
  if (i==33) {
    temp <- read.xlsx("country_data.xlsx",sheetIndex=i,startRow=1,endRow=138)[,-1]
    if (length(global)==1){
      temp <- ts(matrix(temp[,which(is.element(dimnames(temp)[[2]],global))],ncol=1),
                 start=1979,frequency = 4) 
    } else {
      temp <- ts(as.matrix(temp[,which(is.element(dimnames(temp)[[2]],global))]),
                 start=1979,frequency = 4)
    }
    dimnames(temp)[[2]] <- global
    data <- c(data,list(temp))
    names(data)[i+1] <- "exo"
  }
}

#### Construct Europe region
eu <- c("Austria", "Belgium", "Finland", "France", "Germany","Italy", "Netherlands", "Spain")

eu.weights <- read.xlsx("pppgdp.xls",sheetIndex=1,startRow = 1,endRow = 38)[,c(1,34:38)]
eu.weights <- eu.weights[is.element(eu.weights[,1],eu),]
eu.weights <- data.frame("country"=eu,"weight"=rowSums(eu.weights[,-1])/sum(rowSums(eu.weights[,-1])))

pos <- which(is.element(names(data),eu))

var.exist <- c()
for (i in pos) {
  var.exist <- rbind(var.exist,is.element(vars,dimnames(data[[i]])[[2]]))
}
dimnames(var.exist)[[2]] <- vars

eu.data <- c()
for (j in vars){
  w <- eu.weights[var.exist[,j],2]/sum(eu.weights[var.exist[,j],2])
  temp <- c()
  for (i in pos[var.exist[,j]]){
    temp <- cbind(temp,data[[i]][,j])
  }
  temp2 <- 0
  for (i in 1:length(w)){
    temp2 <- temp2 + temp[,i]*w[i]
  }
  eu.data <- cbind(eu.data,temp2)
}
dimnames(eu.data)[[2]] <- vars

dat <- list(eu.data)
names(dat) <- "EuroArea"
for (i in 1:34){
  if (is.element(names(data)[i],eu)) next
  dat <- c(dat,list(data[[i]]))
  names(dat)[length(names(dat))] <- names(data)[i]
}

# Reorder country samples according to country name
data <- list()
names <- names(dat)[-length(dat)]
for (i in order(names)){
  data <- c(data,list(dat[[i]]))
  names(data)[length(data)] <- names[i]
}
data <- c(data,list(dat[[27]]))
names(data)[length(data)] <- "Global"
#save(data,file="gvar_db26.RData")

############################################################################
#################################### Final #################################
############################################################################
dat <- list()
for (i in 1:26){
  temp <- c()
  names <- c()
  # Log GDP
  if (is.element("RGDP",dimnames(data[[i]])[[2]])) {
    tsp.var <- tsp(data[[i]][,"RGDP"])
    temp.data <- (ts(log(data[[i]][,"RGDP"]),start=tsp.var[1],frequency = tsp.var[3]))*100
    temp <- cbind(temp,temp.data);names <- append(names,"y")}
  # Inflation
  if (is.element("CPI",dimnames(data[[i]])[[2]])) {
    tsp.var <- tsp(data[[i]][,"CPI"])
    temp.data <- diff(ts(log(data[[i]][,"CPI"]),start=tsp.var[1],frequency = tsp.var[3]))*100
    temp <- cbind(temp,temp.data);names <- append(names,"Dp")}
  # Short and long term interest rates remain the same
  if (is.element("Rshort",dimnames(data[[i]])[[2]])) {
    tsp.var <- tsp(data[[i]][,"Rshort"])
    temp.data <- .25*log(1+data[[i]][,"Rshort"]/100)*100
    temp <- cbind(temp,temp.data);names <- append(names,"r")}
  
  if (is.element("Rlong",dimnames(data[[i]])[[2]])) {
    tsp.var <- tsp(data[[i]][,"Rlong"])
    temp.data <- .25*log(1+data[[i]][,"Rlong"]/100)*100
    temp <- cbind(temp,temp.data);names <- append(names,"lr")}
  # Real equity prices
  if (is.element("nomEQ",dimnames(data[[i]])[[2]]) && is.element("CPI",dimnames(data[[i]])[[2]])) {
    temp <- cbind(temp,
                 (log(data[[i]][,"nomEQ"]/(data[[i]][,"CPI"])))*100
                 );names <- append(names,"q")}
  # Real exchange rate
  if (is.element("FXdol",dimnames(data[[i]])[[2]]) && is.element("CPI",dimnames(data[[i]])[[2]])) {
    temp <- cbind(temp,
                 (log(data[[i]][,"FXdol"])-log(data[[i]][,"CPI"]))*100
                 #(log(data[[i]][,"FXdol"]/data[[i]][,"CPI"]))#*100
                 );names <- append(names,"ep")}
  temp <- na.omit(temp)
  dimnames(temp)[[2]] <- names
  dat <- c(dat,list(temp))
  names(dat)[i] <- names(data)[i]
  # Append global variables to the list
  if (i==26){
    dat <- c(dat,list(
      ts(matrix(log(data[[27]][-1,])*100,ncol=length(global)), # Minus the first entry due to difference in CPI
                        start=tsp(dat[[1]])[1],frequency=tsp(dat[[1]])[3])
      ))
    dimnames(dat[[27]])[[2]] <- c("poil","prawmat","pmetal")
    names(dat)[i+1] <- "Global"
  }
}
names(dat) <- c("Argentina","Australia","Brazil","Canada","Chile","China","EuroArea","India","Indonesia",
  "Japan","Korea","Malaysia","Mexico","NewZealand","Norway","Peru","Philippines","Safrica",
  "SaudiArabia","Singapore","Sweden","Switzerland","Thailand","Turkey","UK","US","exo")

data <- dat
rm(list=ls()[-which(ls()=="data")])

############################################################################
##################################### Weights ##############################
############################################################################

eu <- c("Austria", "Belgium", "Finland", "France", "Germany","Italy", "Netherlands", "Spain")

weights <- read.xlsx("weights.xls",sheetIndex=4)
row.names(weights) <- weights[,1]
weights <- as.matrix(weights[,-1])

pos <- is.element(dimnames(weights)[[1]],eu)

eu.col <- rowSums(weights[!pos,pos])
eu.row <- colSums(weights[pos,!pos])

weights <- weights[!pos,!pos]

weights <- cbind(eu.col,weights)
weights <- rbind(c(0,eu.row),weights)
weights <- weights/rowSums(weights)

dimnames(weights)[[1]][1] <- "EuroArea"
dimnames(weights)[[2]][1] <- "EuroArea"

weights <- weights[order(dimnames(weights)[[1]]),order(dimnames(weights)[[1]])]
dimnames(weights) <- list(names(data)[-length(names(data))],names(data)[-length(names(data))])

save(list=c("data","weights"),file="gvar_db26.RData")
