#' FRED-QD data
#'
#' The data set contains quarterly time series for 196 US macroeconomic variables from 1959Q3 to 2015Q3.
#' It was produced from file "freddata_Q.csv" of the data sets associated
#' with Chan, Koop, Poirier and Tobias (2019). Raw data are available at
#' \url{https://web.ics.purdue.edu/~jltobias/second_edition/Chapter18/code_for_exercise_4/freddata_Q.csv}.
#' 
#' @usage data("bem_dfmdata")
#' 
#' @format A named time-series object with 225 rows and 196 variables. A detailed explanation of the
#' variables can be found in McCracken and Ng (2016).
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#' 
#' McCracken, M. W., & Ng, S. (2016). FRED-MD: A monthly database for macroeconomic research.
#' \emph{Journal of Business & Economic Statistics 34}(4), 574-589. \doi{10.1080/07350015.2015.1086655}
#' 
#' 
#' 
"bem_dfmdata"

#' West German economic time series data
#'
#' The data set contains quarterly, seasonally adjusted time series for West German fixed investment, disposable
#' income, and consumption expenditures in billions of DM from 1960Q1 to 1982Q4. It was produced
#' from file E1 of the data sets associated with Lütkepohl (2007). Raw data are available at
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
#' with Lütkepohl (2007). Raw data are available at \url{http://www.jmulti.de/download/datasets/e6.dat}
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