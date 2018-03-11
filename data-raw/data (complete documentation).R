# ' GVAR Dataset (1979Q1-2003Q4)
# '
# ' The dataset and model specifications used for the GVAR model in Dees et al. (2007).
# '
# ' @format A list containing the following elements:
# ' \describe{
# '   \item{country.data}{a list of 26 "zoo" objects containing 99 observations of macroeconomic variables.}
# '   \item{global.data}{a "zoo" object containing 99 observations of oil price data.}
# '   \item{weight.data}{a matrix of country weights, where each row sums up to 1.}
# '   \item{country.specifications}{a list of country model specifications as used in Dees et al. (2007).}
# ' }
# ' 
# ' @source \url{http://qed.econ.queensu.ca/jae/2007-v22.1/dees-di_mauro-pesaran-smith/}
# ' @references
# ' Dees, S., Mauro, F. d., Pesaran, M. H., & Smith, L. V. (2007). Exploring the international linkages of the euro area: A global VAR analysis. \emph{Journal of Applied Econometrics}, 22(1), 1--38.
#"dees2007"

# ' GVAR Dataset (1979Q1-2013Q1)
# '
# ' A dataset produced from the GVAR dataset (Vintage 2013) of the GVAR Toolbox 2.0.
# '
# ' @format A list containing the following elements:
# ' \describe{
# '   \item{country.data}{a list of 26 "zoo" objects containing 136 observations of macroeconomic variables.}
# '   \item{global.data}{a "zoo" object containing 136 observations oil, raw materials and metal price data.}
# '   \item{weight.data}{a matrix of country weights, where each row sums up to 1.}
# ' }
# ' 
# ' @source \url{https://sites.google.com/site/gvarmodelling/data}
# ' @references
# ' Smith, L. and A. Galesi (2014). GVAR Toolbox 2.0. University of Cambridge: Judge Business School, \url{https://sites.google.com/site/gvarmodelling/gvar-toolbox}.
# "gvar2013"

#' Global VAR (GVAR) Database,  1979Q2-2016Q4
#'
#' GVAR dataset (2016 Vintage) produced from data of the GVAR Toolbox 2.0.
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{country.data}{a list of 33 "zoo" objects containing 151 observations of macroeconomic variables.}
#'   \item{global.data}{a "zoo" object containing 151 observations oil, raw materials and metal price data.}
#'   \item{region.weights}{a "zoo" object containing PPP-GDP time series for 33 countries from 1990 to 2016.}
#'   \item{weight.data}{an array of trade weights for 33 countries and yearly observations from 1980 to 2016.}
#' }
#' 
#' @source \url{https://sites.google.com/site/gvarmodelling/data}
#' @references
#' Mohaddes, K. & Raissi, M. (2018). Compilation, Revision and Updating of the Global VAR (GVAR) Database, 1979Q2-2016Q4. University of Cambridge: Faculty of Economics (mimeo).
"gvar2016"