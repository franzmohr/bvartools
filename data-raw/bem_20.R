rm(list = ls())

# Specify URL
#url <- "https://web.ics.purdue.edu/~jltobias/second_edition/Chapter20/code_for_exercise_1/US_macrodata.csv"
url <- "https://web.ics.purdue.edu/~jltobias/second_edition/Chapter20/code_for_exercise_7/US_macrodata1.csv"
# Load data
download.file(url, destfile = "data-raw/US_macrodata1.csv")
data <- read.csv("data-raw/US_macrodata1.csv", col.names = c("date", "Dp", "u", "r"))
# Omit NA values
data <- na.omit(data)
# Transform into time-series object
us_macrodata <- ts(data[, -1], start = c(1959, 2), frequency = 4)

usethis::use_data(us_macrodata, overwrite = TRUE)
