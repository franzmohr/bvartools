#' Stochastic Volatility
#'
#' Produces a draw of log-volatilities.
#'
#' @param y a \eqn{T \times K} matrix containing the time series.
#' @param h a \eqn{T \times K} vector of the current draw of log-volatilities.
#' @param sigma a \eqn{K \times 1} vector of variances of log-volatilities,
#' where the \eqn{i}th element corresponds to the \eqn{i}th column in \code{y}.
#' @param h_init a \eqn{K \times 1} vector of the initial states of log-volatilities,
#' where the \eqn{i}th element corresponds to the \eqn{i}th column in \code{y}.
#' @param constant a \eqn{K \times 1} vector of constants that should be added to \eqn{y^2}
#' before taking the natural logarithm. The \eqn{i}th element corresponds to
#' the \eqn{i}th column in \code{y}. See 'Details'.
#' 
#' @details For each column in \code{y} the function produces a posterior
#' draw of the log-volatility \eqn{h} for the model
#' \deqn{y_{t} = e^{\frac{1}{2}h_t} \epsilon_{t},}
#' where \eqn{\epsilon_t \sim N(0, 1)} and \eqn{h_t} is assumed to evolve according to a random walk
#' \deqn{h_t = h_{t - 1} + u_t,}
#' with \eqn{u_t \sim N(0, \sigma^2)}.
#' 
#' The implementation is based on the algorithm of Kim, Shephard and Chip (1998) and performs the
#' following steps:
#' \enumerate{
#'   \item Perform the transformation \eqn{y_t^* = ln(y_t^2 + constant)}.
#'   \item Obtain a sample from the seven-component normal mixture for
#'   approximating the log-\eqn{\chi_1^2} distribution.
#'   \item Obtain a draw of log-volatilities.
#' }
#' 
#' The implementation follows the code provided on the website to the textbook
#' by Chan, Koop, Poirier, and Tobias (2019).
#' 
#' @return A vector of log-volatility draws.
#' 
#' @references
#' 
#' Chan, J., Koop, G., Poirier, D. J., & Tobias J. L. (2019). \emph{Bayesian econometric methods}
#' (2nd ed.). Cambridge: Cambridge University Press.
#' 
#' Kim, S., Shephard, N., & Chib, S. (1998). Stochastic volatility. Likelihood inference and comparison
#' with ARCH models. \emph{Review of Economic Studies 65}(3), 361--393. \doi{10.1111/1467-937X.00050}
#' 
#' @examples
#' data("us_macrodata")
#' y <- diff(us_macrodata)
#' h_init <- log(diag(var(y)))
#' h <- t(matrix(h_init, 3, nrow(y)))
#' sigma_h <- rep(.05, 3)
#' const <- rep(.0001, 3)
#' stochvol_ksc1998(y, h, sigma_h, h_init, const)
#' 
#' 
#' @export
stochvol_ksc1998 <- function(y, h, sigma, h_init, constant) {
  
  k <- ncol(y)
  tt <- nrow(y)
  
  if (any(is.na(y))) {
    stop("Argument 'y' contains NAs.")
  }
  if (k != ncol(h)) {
    stop("Arguments 'y' and 'h' do not have the same number of columns.")
  }
  if (tt != nrow(h)) {
    stop("Arguments 'y' and 'h' do not have the same number of rows.")
  }
  
  hh <- create_first_difference_matrix(1, tt)
  hh <- t(hh) %*% hh
  
  # Prepare series
  y <- log(y * y + matrix(constant, nrow = tt, ncol = k, byrow = TRUE))
  
  sigs <- Matrix::Diagonal(tt, 1)
  p_i <- matrix(c(0.00730, 0.10556, 0.00002, 0.04395, 0.34001, 0.24566, 0.25750), tt, 7, byrow = TRUE)
  mu <- c(-10.12999, -3.97281, -8.56686, 2.77786, 0.61942, 1.79518, -1.08819) - 1.2704
  mu_large <- matrix(mu, tt, 7, byrow = TRUE)
  sigma2 <- c(5.79596, 2.61369, 5.17950, 0.16735, 0.64009, 0.34023, 1.26261)
  sigma_large <- matrix(sqrt(sigma2), tt, 7, byrow = TRUE)
  
  for (i in 1:k) {
    
    pdf_matrix <- cbind(matrix(y[,i], tt, 7), matrix(h[,i], tt, 7) + mu_large, sigma_large)
    pdf_matrix <- t(apply(pdf_matrix, 1, function(x) {stats::dnorm(x[1:7], x[8:14], x[15:21])}))
    
    q <- p_i * pdf_matrix
    q <- t(apply(q / matrix(rowSums(q), tt, 7), 1, cumsum))
    s <- 8 - rowSums(matrix(stats::runif(tt), tt, 7) < q)
    
    sigh_hh <- hh / sigma[i]
    Matrix::diag(sigs) <- 1 / sigma2[s]
    post_h_v <- sigh_hh + sigs
    post_h_mu <- solve(post_h_v, sigh_hh %*% matrix(1, tt) * h_init[i] + sigs %*% (y[,i] - mu[s]));
    h[, i] <- matrix(post_h_mu + solve(chol(post_h_v), stats::rnorm(tt)))
    
  }

  return(h)
}