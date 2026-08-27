rm(list = ls())

# Load data
raw <- read.table("data-raw/uk_macrodata.txt")
uk_macrodata <- ts(raw[, -1], start = c(1957, 2), frequency = 4)
dimnames(uk_macrodata) <- list(NULL, c("rs", "rl", "Dp"))
plot(uk_macrodata)

usethis::use_data(uk_macrodata, overwrite = TRUE)