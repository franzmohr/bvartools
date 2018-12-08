rm(list = ls())

e1 <- ts(read.table("http://www.jmulti.de/download/datasets/e1.dat", skip = 6, header = TRUE), start = 1960, frequency = 4)
usethis::use_data(e1, overwrite = TRUE)
