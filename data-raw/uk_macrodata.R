rm(list = ls())

uk <- read.delim("~/Dropbox/R/bvartools/data-raw/data_uk.txt", header = FALSE,
                        col.names = c("date", "u", "r", "Dp"))

uk_macrodata <- ts(uk[, c("Dp", "u", "r")], start = uk[1, 1], frequency = 4)


usethis::use_data(uk_macrodata, overwrite = TRUE)

#plot(uk_macrodata)
