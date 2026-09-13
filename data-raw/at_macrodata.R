rm(list = ls())

# The Austrian sub-model of a GVAR model on the database of Mohaddes and
# Raissi (2024), which comes with bgvars as data set 'gvar2023'.
data("gvar2023", package = "bgvars")
submodel_data <- gvar2023[["submodel_data"]]
global_data <- gvar2023[["global_data"]]

# Initialize global model
object <- bgvars::create_gvarmodel(submodel_data = submodel_data,
                                   global_data = global_data)

# The foreign (star) variables of Austria are trade weighted averages of the
# series of the other 32 countries. Weights are rolling sums of trade flows
# over the last three years, as in the GVAR literature.
object <- bgvars::add_weight_matrices(object = object,
                                      submodel_data = submodel_data,
                                      period = 3)

submodel <- bgvars::create_varxsubmodel(object, submodel = "AT",
                                        p_endogen = 1, p_exogen = 1,
                                        iterations = 10, burnin = 10)[[1]]

domestic <- submodel[["data"]][["original"]][["endogen"]]
foreign <- submodel[["data"]][["original"]][["exogen"]]

at_macrodata <- cbind(domestic, foreign, global_data)
dimnames(at_macrodata) <- list(NULL, c(dimnames(domestic)[[2]],
                                      dimnames(foreign)[[2]],
                                      dimnames(global_data)[[2]]))
plot(at_macrodata[, 1:8])
plot(at_macrodata[, 9:15])

usethis::use_data(at_macrodata, overwrite = TRUE)
