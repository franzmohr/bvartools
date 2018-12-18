rm(list = ls())

download.file("http://www.jmulti.de/download/datasets/e1.dat",
              destfile = "data-raw/e1.dat")
e1 <- ts(read.table("data-raw/e1.dat", skip = 6, header = TRUE),
         start = 1960, frequency = 4)
usethis::use_data(e1, overwrite = TRUE)

download.file("http://www.jmulti.de/download/datasets/e6.dat",
              destfile = "data-raw/e6.dat")
e6 <- ts(read.table("data-raw/e6.dat", skip = 10, header = TRUE),
         start = 1972.25, frequency = 4)
e6 <- e6[, c("R", "Dp")]
usethis::use_data(e6, overwrite = TRUE)
