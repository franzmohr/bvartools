rm(list=ls())

################### Load functions ###################
devtools::load_all(".")
library(bgvars)
library(ggplot2)
################### Load data ###################
data("gvar2016")

region.weights <- gvar2016$region.weights
country.data <- gvar2016$country.data
global.data <- gvar2016$global.data
weight.data <- gvar2016$weight.data
rm(gvar2016)


except <- c("r", "lr")
for (j in names(country.data)) {
  temp <- country.data[[j]]
  for (i in 1:NCOL(temp)) {
    if(dimnames(temp)[[2]][i] %in% except) {
      next
    } else {
      d_temp <- diff(temp[,i])
      #sm <- loess(d_temp ~ matrix(1:length(d_temp)))
      #res <- residuals(sm)
      res <- d_temp
      
      limits <- quantile(res, probs = c(.25, .75))
      iqr <- limits[2] - limits[1]
      correct <- res < limits[1] - 1.5 * iqr
      correct[res > limits[2] + 1.5 * iqr] <- TRUE
      names(correct) <- NULL
      y <- temp[,i]
      t <- length(y)
      pos <- which(correct) + 1
      if (any(pos == t)) {
        pos <- pos[-length(pos)]
      }
      temp[pos,i] <- NA
      temp[, i] <- zoo::na.approx(temp[,i])
      #print(plot.ts(diff(cbind(y, temp[,i])), plot.type = "single", col = c("blue", "red"))) 
    }
  }
  country.data[[j]] <- temp
}

weight.data <- create_weights(weight.data = weight.data, period = 4, country.data = country.data)

# Create regions
regions <- list("EuroArea" = c("Austria", "Belgium", "Finland", "France",
                               "Germany", "Italy", "Netherlands", "Spain"))
regional.data <- create_regions(country.data, regions, 4, region.weights, weight.data)

# Final sample
country.data <- regional.data$country.data
weight.data <- regional.data$weight.data
countries <-  NULL

############### GVAR options ####################
type <- "VEC"
use <- c("y", "Dp", "r", "poil")
iterations <- 30000
burnin <- 15000
thin <- 10

#### Test lag specification with full rank ####
# Basic country configurations
country.specs <- country_specifications(country.data = country.data, global.data = global.data,
                                        p = 1:2, p.star = 1:2, p.global = 1,
                                        case = "IV",
                                        countries = countries,
                                        use = use,
                                        type = type, r = NULL,
                                        structural = FALSE, tvp = c(FALSE, FALSE, FALSE), sv = FALSE)

# Additional specifications for individual countries
country.specs$USA$domestic.variables <- c("poil", "y", "Dp", "r")

# Generate standard priors for all models
prior <- standard_priors()
prior$Omega$constant <- c(0, 0)

# Generate final data for estimation
data <- country_models(country.data, weight.data, global.data, country.specs, prior = prior)

# Estimate
data <- gvar_fit(data, iterations = iterations, burnin = burnin, thin = thin)

# Solve the GVAR model
gvar_aic <- solve_gvar(data, ic = "AIC", t = 149)

# GIRF
temp <- c()
for (i in c("r")) {
  for (j in c("USA", "EuroArea", "Japan")) {
    for (k in c("y", "Dp", "r")) {
      irf_data <- girf(data = gvar_aic, n.ahead = 24, impulse = c("USA", i), response = c(j, k),
                       shock = .0025, ci = c(.9, .68), mean = FALSE, cumulative = FALSE)
      
      irf_temp <- cbind( "Quarter" = 0:(nrow(irf_data$IR) - 1), as.data.frame(irf_data$IR * 100),
                         "impulse_country" = "USA", "impulse_variable" = i,
                         "response_country" = j, "response_variable" = k)
      temp <- rbind(temp, irf_temp) 
    }
  }
}

ggplot(temp, aes(x = Quarter, y = `50%`)) +
  geom_ribbon(aes(ymin = `5%`, ymax = `95%`), fill = "gray70") +
  geom_ribbon(aes(ymin = `16%`, ymax = `84%`), fill = "gray50") +
  geom_hline(yintercept = 0, colour = "black") +
  geom_line(colour = "black") +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(response_variable ~ response_country, scales = "free_y") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(colour = "black"))