
#' West German economic time series data
#'
#' The data set contains quarterly, seasonally adjusted time series for West German fixed investment, disposable
#' income, and consumption expenditures in billions of DM from 1960Q1 to 1982Q4. It was produced
#' from file E1 of the data sets associated with Lütkepohl (2006). Raw data are available at
#' \url{http://www.jmulti.de/download/datasets/e1.dat} and were originally obtained from
#' Deutsche Bundesbank.
#' 
#' @usage data("e1")
#' 
#' @format A named time-series object with 92 rows and 3 variables:
#' \describe{
#'   \item{invest}{fixed investment.}
#'   \item{income}{disposable income.}
#'   \item{cons}{consumption expenditures.}
#' }
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
"e1"

#' German interest and inflation rate data
#'
#' The data set contains quarterly, seasonally unadjusted time series for German long-term interest
#' and inflation rates from 1972Q2 to 1998Q4. It was produced from file E6 of the data sets associated
#' with Lütkepohl (2006). Raw data are available at \url{http://www.jmulti.de/download/datasets/e6.dat}
#' and were originally obtained from Deutsche Bundesbank and Deutsches Institut für Wirtschaftsforschung.
#' 
#' @usage data("e6")
#' 
#' @format A named time-series object with 107 rows and 2 variables:
#' \describe{
#'   \item{R}{nominal long-term interest rate (Umlaufsrendite).}
#'   \item{Dp}{\eqn{\Delta} log of GDP deflator.}
#' }
#' 
#' @details The data cover West Germany until 1990Q2 and all of Germany aferwards. The values refer to 
#' the last month of a quarter.
#' 
#' @references
#' 
#' Lütkepohl, H. (2006). \emph{New introduction to multiple time series analysis} (2nd ed.). Berlin: Springer.
#' 
"e6"

#' US macroeconomic data
#'
#' The data set contains quarterly time series for the US CPI inflation rate, unemployment rate, and
#' Fed Funds rate from 1959Q2 to 2007Q4. It was produced from file "US_macrodata.csv" of the data sets associated
#' with Chan, Koop, Poirier and Tobias (2019). Raw data are available at
#' \url{https://web.ics.purdue.edu/~jltobias/second_edition/Chapter20/code_for_exercise_1/US_macrodata.csv}.
#' 
#' @usage data("us_macrodata")
#' 
#' @format A named time-series object with 195 rows and 3 variables:
#' \describe{
#'   \item{Dp}{CPI inflation rate.}
#'   \item{u}{unemployment rate.}
#'   \item{r}{Fed Funds rate.}
#' }
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#' 
"us_macrodata"

#' UK interest and inflation rate data
#'
#' The data set contains quarterly time series for two UK interest rates and inflation
#' from 1957Q4 to 2009Q3. It was produced from the supplementary material to
#' Koop, León-González and Strachan (2011).
#' 
#' @usage data("uk_macrodata")
#' 
#' @format A named time-series object with 210 rows and 3 variables:
#' \describe{
#'   \item{rs}{short term rate: Treasury bills. This should be series 'FITB_PA' of the IMF's International Financial Statistics database.}
#'   \item{rs}{long term rate: Government bonds. This should be series 'FIGB_PA' of the IMF's International Financial Statistics database.}
#'   \item{Dp}{inflation: Quarterly change in the log CPI annualized by multiplying the change by 400.}
#' }
#' 
#' @references
#' 
#' Koop, G., León-González, R., & Strachan R. W. (2011). Bayesian inference in
#' a time varying cointegration model. \emph{Journal of Econometrics, 165}(2), 210--220.
#' \doi{10.1016/j.jeconom.2011.07.007}
#' 
"uk_macrodata"

#' Austrian sub-model of a global VAR
#'
#' The data set contains quarterly time series for Austria, the corresponding foreign
#' (star) variables and global commodity prices from 1979Q2 to 2023Q3. It was produced
#' from data set \code{gvar2023} of package \code{bgvars}, which contains the GVAR database
#' of Mohaddes and Raissi (2024), as the Austrian sub-model of a global VAR model
#' of all 33 countries of the database.
#' 
#' @usage data("at_macrodata")
#' 
#' @format A named time-series object with 178 rows and 15 variables:
#' \describe{
#'   \item{y}{log real GDP.}
#'   \item{Dp}{rate of inflation, the quarterly change in the log CPI.}
#'   \item{eq}{log real equity prices.}
#'   \item{ep}{log exchange rate against the US dollar, deflated by the CPI.}
#'   \item{r}{short-term interest rate, \eqn{0.25 ln(1 + R^{S} / 100)}.}
#'   \item{lr}{long-term interest rate, \eqn{0.25 ln(1 + R^{L} / 100)}.}
#'   \item{y.s, Dp.s, eq.s, ep.s, r.s, lr.s}{foreign counterparts of the domestic variables.}
#'   \item{poil}{log of oil prices.}
#'   \item{pmat}{log of agricultural raw material prices.}
#'   \item{pmetal}{log of metals prices.}
#' }
#' 
#' @details The foreign variables are trade weighted averages of the series of the other
#' countries of the database, for which the respective variable is available.
#' The weights are the shares of the countries in Austria's trade over the current and
#' the two preceding years. Trade data are available from 1980 to 2016, so the weights of the years
#' 1979 to 1982 are those of 1980 to 1982 and the weights from 2016 onwards are those of
#' 2014 to 2016.
#' 
#' In a sub-model of a global VAR the domestic variables are endogenous, and the foreign
#' and global variables are weakly exogenous. They can be passed to argument \code{exogen}
#' of \code{\link{create_bvarmodel}} or \code{\link{create_bvecmodel}}.
#' 
#' @references
#' 
#' Dees, S., di Mauro, F., Pesaran, M. H., & Smith, L. V. (2007). Exploring the international
#' linkages of the euro area: A global VAR analysis. \emph{Journal of Applied Econometrics, 22}(1), 1--38.
#' \doi{10.1002/jae.932}
#' 
#' Mohaddes, K., & Raissi, M. (2024). \emph{Compilation, revision and updating of the global VAR (GVAR)
#' database, 1979Q2--2023Q3} (mimeo). University of Cambridge: Judge Business School.
#' 
"at_macrodata"
