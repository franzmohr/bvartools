# Downloaded Feb 19, 2018.
# https://sites.google.com/site/gvarmodelling/gvar-toolbox/download

rm(list = ls())

library(readxl)

# PPP Data
PPP <- as.data.frame(read_xls("data-raw/gvar2016/PPP-GDP WDI (1990-2016).xls", sheet = "WDI"))
PPP <- na.omit(PPP)[, -which(names(PPP) %in% c("Country Code", "Indicator Name", "Indicator Code"))]
nam <- PPP[, "Country Name"]
time <- names(PPP)[-1]
PPP <- data.frame(time, t(PPP[,-1]))
nam[which(nam == "Korea, Rep.")] <- "Korea"
nam[which(nam == "United Kingdom")] <- "UK"
nam[which(nam == "United States")] <- "USA"
names(PPP) <- c("Year", nam)
rownames(PPP) <- NULL
PPP <- zoo::zoo(PPP[, -1], order.by = PPP$Year)

# Country data
variables <- c("y", "Dp", "eq", "ep", "r", "lr")
countries <- as.data.frame(read_xls("data-raw/gvar2016/Country Codes.xls", col_types = "text"))

data <- c()
nam <- c()
for (i in countries[, "Country Name"]) {
  temp <- as.data.frame(read_xls("data-raw/gvar2016/Country Data (1979Q2-2016Q4).xls", sheet = i))
  temp <- zoo::zoo(temp[, which(names(temp)%in%variables)], order.by = zoo::as.yearqtr(temp[, "date"]))
  data <- c(data, list(temp))
  nam <- c(nam, i)
  rm(temp)
}

capcase <- function(x) { # from ?tolower
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

nam <- unlist(lapply(tolower(nam), capcase))
countries <- cbind(countries, "Name" = nam)
countries$Name <- as.character(countries$Name)
countries$Name[which(countries$Name == "United Kingdom")] <- "UK"
countries$Name[which(countries$Name == "Usa")] <- "USA"
names(data) <- countries$Name
country.data <- data

# Global data
global.variables <- c("poil", "pmat", "pmetal")
global.data <- as.data.frame(read_xls("data-raw/gvar2016/Country Data (1979Q2-2016Q4).xls", sheet = 1))
global.data <- zoo::zoo(global.data[, which(names(global.data)%in%global.variables)], order.by = zoo::as.yearqtr(global.data[, "date"]))

# Weight matrix
weights <- array(NA, dim = c(length(country.data), length(country.data), 37))
dimnames(weights) <- list(countries[, "Country Code"], countries[, "Country Code"], as.character(1980:2016))
for (i in countries[, "Country Code"]) {
  temp <- as.data.frame(read_xlsx("data-raw/gvar2016/Trade Flows (1980-2016).xlsx", sheet = i, na = "NaN"))
  dimnames(temp)[[1]] <- temp[, 1]
  temp <- temp[, -(1:2)]
  temp <- temp[ , countries[, "Country Code"]]
  temp[, i] <- 0
  weights[i , ,] <- t(as.matrix(temp))
}
#for (i in 1:dim(weights)[3]) {
#  weights[,, i] <- t(apply(weights[,, i], 1, function(x){s <- sum(x); return(x/s)}))
#}
dimnames(weights)[[1]] <- countries$Name
dimnames(weights)[[2]] <- countries$Name
weight.data <- weights

gvar2016 <- list("country.data" = country.data,
                 "global.data" = global.data,
                 "region.weights" = PPP,
                 "weight.data" = weight.data)
save(gvar2016, file = "data/gvar2016.rda")
