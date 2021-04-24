

x <- read.csv(file = "data-raw/bem_dfmdata.csv", header = FALSE)
x[is.na(x)] <- NA

fred <- read.csv(file = "data-raw/fred_qd.csv")[-c(1, 2), -1]

pos <- which(apply(x, 2, function(x) {all(!is.na(x))}))

bem_dfmdata <- ts(x[, pos], start = c(1959, 3), frequency = 4, names = names(fred)[pos])

usethis::use_data(bem_dfmdata, overwrite = TRUE)
