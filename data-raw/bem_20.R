rm(list = ls())

# Specify URL
url <- "https://web.ics.purdue.edu/~jltobias/second_edition/Chapter20/code_for_exercise_1/US_macrodata.csv"
#url <- "https://web.ics.purdue.edu/~jltobias/second_edition/Chapter20/code_for_exercise_7/US_macrodata1.csv"
# Load data
temp <- tempfile()
download.file(url, destfile = temp)
data <- read.csv(temp, col.names = c("date", "Dp", "u", "r"))
unlink(temp)
# Omit NA values
data <- na.omit(data)
# Transform into time-series object
us_macrodata <- ts(data[, -1], start = c(1959, 2), frequency = 4)

usethis::use_data(us_macrodata, overwrite = TRUE)
