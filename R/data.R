#' GVAR Database, 1979Q2-2016Q4
#'
#' Data set containing macroeconomic time series for 33 countries and 3 commodities from 1979 Q2 to 2016 Q4.
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{country.data}{a named list of time series containing 151 quarterly observations of macroeconomic variables for 33 countries.}
#'   \item{global.data}{a time series containing 151 quarterly observations of oil, raw materials and metal prices.}
#'   \item{region.weights}{a time series containing yearly data of PPP-GDP for 33 countries from 1990 to 2016.}
#'   \item{weight.data}{an array of trade weights for 33 countries with yearly observations from 1980 to 2016.}
#' }
#' 
#' @source \url{https://sites.google.com/site/gvarmodelling/data}
#' 
#' @references
#' Mohaddes, K. & Raissi, M. (2018). Compilation, Revision and Updating of the Global VAR (GVAR) Database, 1979Q2-2016Q4. University of Cambridge: Faculty of Economics (mimeo).
"gvar2016"