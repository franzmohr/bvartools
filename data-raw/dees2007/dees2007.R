# Downloaded March 1, 2018.
# http://qed.econ.queensu.ca/jae/2007-v22.1/dees-di_mauro-pesaran-smith/

rm(list = ls())

library(readxl)

# Region weights
region.weights <- as.data.frame(read_xls("data-raw/dees2007/pppgdp9901.XLS", sheet = "Sheet1"))
nam <- region.weights[, 1]
region.weights <- t(region.weights[,-1])
nam[which(nam == "Korea, Rep.")] <- "Korea"
nam[which(nam == "United Kingdom")] <- "UK"
nam[which(nam == "United States")] <- "USA"
dimnames(region.weights) <- list(NULL, nam)
region.weights <- as.data.frame(region.weights)
region.weights <- zoo::zoo(PPP, order.by = 1999:2001)

# Country data
global.variables <- "poil"
countries <- as.data.frame(read_xls("data-raw/dees2007/countrycodes.xls", sheet = "Sheet1", col_types = "text"))
capcase <- function(x) { # from ?tolower
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
countries[,1] <- unlist(lapply(tolower(countries[,1]), capcase))
countries[which(countries[,1] == "United Kingdom"), 1] <- "UK"
countries[which(countries[,1] == "Usa"), 1] <- "USA"

data <- c()
nam <- c()
for (i in 1:nrow(countries)) {
  temp <- as.data.frame(read_xls("data-raw/dees2007/countrydata33.xls", sheet = i))
  names(temp) <- tolower(names(temp))
  temp <- temp[, !apply(temp, 2, function(x){length(unique(x)) == 1})]
  if (any(names(temp) %in% global.variables)) {
    temp <- temp[, -which(names(temp)%in%global.variables)]
  }
  if (any(names(temp) %in% "eindex")) {names(temp)[which(names(temp) == "eindex")] <- "e"}
  temp <- zoo::zoo(temp[, -which(names(temp) == "date")], order.by = zoo::as.yearqtr(temp[, "date"]))
  data <- c(data, list(temp))
  nam <- c(nam, countries[i, 1])
  rm(temp)
}

names(data) <- countries[,1]
country.data <- data

# Global data
global.variables <- "POIL"
global.data <- as.data.frame(read_xls("data-raw/dees2007/countrydata33.xls", sheet = 1))
global.data <- zoo::zoo(matrix(global.data[, which(names(global.data)%in%global.variables)]), order.by = zoo::as.yearqtr(global.data[, "date"]))
names(global.data) <- "poil"

# Weight matrix
weights <- array(NA, dim = c(length(country.data), length(country.data), 37))
dimnames(weights) <- list(countries[, "Country Code"], countries[, "Country Code"], as.character(1980:2016))
for (i in countries[, "Country Code"]) {
  temp <- as.data.frame(read_xlsx("data-raw/dees2007/tradematrix8003.xls", sheet = i, na = c("")))
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
