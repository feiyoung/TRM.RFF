#' Generate data from an autoregressive threshold regression model
#'
#' @description
#' Generates simulated data from a threshold regression model with a
#' state-dependent regime boundary, conditionally heteroskedastic regressors,
#' and AR(1) disturbances. The default setting corresponds to a linear
#' regime boundary given by \eqn{f_0(Z_t) = Z_t + 0.3}.
#'
#' @param n A positive integer specifying the sample size.
#' @param alpha A numeric vector of length 2 containing the baseline
#' regression coefficients. Default is \code{c(1, 2)}.
#' @param beta A numeric vector of length 2 containing the regime-shift
#' coefficients. Default is \code{c(-1.5, 1.5)}.
#' @param rho_eps AR(1) coefficient of the disturbance process.
#' Default is \code{0.3}.
#' @param sigma_eps Standard deviation of the AR(1) innovations.
#' Default is \code{0.1}.
#' @param boundary_shift Constant shift in the regime boundary
#' \eqn{f_0(Z_t) = Z_t + c}. Default is \code{0.3}.
#' @param burn_in Number of burn-in observations used when generating
#' the AR(1) disturbance. Default is \code{100}.
#' @param seed An integer specifying the random seed for reproducibility.
#' Default is \code{1}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{Y}: simulated response vector of length \code{n}.
#'   \item \code{X}: an \code{n x 2} regression design matrix.
#'   \item \code{Z}: an \code{n x 1} matrix of state variables determining
#'   the regime boundary.
#'   \item \code{U}: observed threshold variable.
#'   \item \code{D0}: true regime indicators.
#'   \item \code{f0}: true boundary values evaluated at \code{Z}.
#'   \item \code{eps}: simulated AR(1) disturbances.
#'   \item \code{alpha}: true baseline regression coefficients.
#'   \item \code{beta}: true regime-shift coefficients.
#'   \item \code{theta}: combined true coefficient vector
#'   \code{c(alpha, beta)}.
#' }
#'
#' @details
#' The data are generated from
#' \deqn{
#' Y_t = X_t^\top \alpha_0 +
#'       (X_t^\top \beta_0)D_{0t} + \varepsilon_t,
#' }
#' where
#' \deqn{
#' D_{0t} = I\{f_0(Z_t) \ge U_t\},
#' \qquad
#' f_0(Z_t) = Z_t + c.
#' }
#' The threshold variable \eqn{U_t} and state variable \eqn{Z_t} are
#' independent standard normal variables. The second component of
#' \eqn{X_t} is generated with a variance depending on \eqn{(U_t,Z_t)},
#' while the disturbance follows the stationary AR(1) process
#' \deqn{
#' \varepsilon_t = \rho_\varepsilon \varepsilon_{t-1} + \nu_t,
#' \qquad
#' \nu_t \sim N(0,\sigma_\varepsilon^2).
#' }
#'
#' @export
#' @importFrom stats rnorm filter
#'
#' @examples
#' dat <- generate_AR(n = 500, seed = 1)
#'
#' str(dat)
#' table(dat$D0)
#'
#' \dontrun{
#' fit <- TRM.RFF(
#'   Y = dat$Y,
#'   X = dat$X,
#'   Z = dat$Z,
#'   U = dat$U,
#'   rff_dim = 15,
#'   h = 0.025,
#'   sigma_h = 1
#' )
#' summary(fit)
#' }
generate_AR <- function(
    n,
    alpha = c(1, 2),
    beta = c(-1.5, 1.5),
    rho_eps = 0.3,
    sigma_eps = 0.1,
    boundary_shift = 0.3,
    burn_in = 100L,
    seed = 1
) {

  if (n <= 0) stop("n must be positive.")
  if (length(alpha) != 2L || length(beta) != 2L)
    stop("alpha and beta must have length 2.")
  if (abs(rho_eps) >= 1)
    stop("rho_eps must satisfy abs(rho_eps) < 1.")

  set.seed(seed)

  ## Threshold and state variables
  ZU <- matrix(rnorm(2 * n), ncol = 2)
  U <- ZU[, 1]
  Z <- matrix(ZU[, 2], ncol = 1)

  ## Regression covariates
  X <- cbind(
    1,
    rnorm(n, sd = 1 / sqrt(1 + rowSums(ZU^2) / 2))
  )

  ## True regime boundary
  f0 <- as.numeric(Z) + boundary_shift
  D0 <- as.integer(f0 >= U)

  ## AR(1) disturbances
  innov <- rnorm(n + burn_in, sd = sigma_eps)
  eps_full <- as.numeric(
    filter(innov, filter = rho_eps, method = "recursive")
  )
  eps <- eps_full[(burn_in + 1L):(n + burn_in)]

  ## Response
  Y <- as.numeric(
    X %*% alpha + (X %*% beta) * D0 + eps
  )

  colnames(X) <- c("Intercept", "X1")
  colnames(Z) <- "Z1"

  list(
    Y = Y,
    X = X,
    Z = Z,
    U = U,
    D0 = D0,
    f0 = f0,
    eps = eps,
    alpha = alpha,
    beta = beta,
    theta = c(alpha, beta)
  )
}
