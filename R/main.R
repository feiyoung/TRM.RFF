# generate man files
# devtools::document()
# R CMD check --as-cran MultiCOAP_1.1.tar.gz
## usethis::use_data(dat_r2_mac)
# pkgdown::build_site()
# pkgdown::build_home()
# pkgdown::build_reference()
# pkgdown::build_article("COAPsimu")
# pkgdown::build_article("ProFASTdlpfc2")

# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('html_document'))

# rmarkdown::render('./vignettes_PDF/COAPsimu.Rmd', output_format=c('pdf_document'), clean = F)



#' Fit a nonlinear threshold regression model using random Fourier features
#'
#' @description
#' Fits a threshold regression model with a nonlinear multivariate regime
#' boundary approximated by random Fourier features (RFF). The procedure first
#' estimates the boundary using a smoothed threshold criterion and then refits
#' the regression coefficients using the estimated hard regime classification.
#' HAC standard errors are computed for the refitted coefficients.
#'
#' @param Y A numeric vector of length \code{n} containing the response variable.
#' @param X An \code{n x d} matrix of regressors.
#' @param Z An \code{n x p} matrix of variables determining the regime boundary.
#' @param U A numeric vector of length \code{n} containing the observed threshold variable.
#' @param rff_dim A positive integer specifying the number of random Fourier features.
#' @param h A positive smoothing bandwidth for the threshold indicator.
#' @param sigma_h A positive bandwidth parameter used to generate the random Fourier features.
#' @param seed An integer specifying the random seed used to generate the random Fourier features.
#' Default is 1.
#' @param lambda A nonnegative regularization parameter for boundary estimation.
#' Default is \code{1e-10}.
#' @param n_starts A positive integer specifying the number of initializations
#' used in the first-stage optimization. Default is 7.
#' @param conf_level Confidence level for the coefficient confidence intervals.
#' Default is 0.95.
#'
#' @return An object of class \code{"TRM.RFF"} containing:
#' \itemize{
#'   \item \code{alpha}: estimated baseline regression coefficients.
#'   \item \code{beta}: estimated regime-shift coefficients.
#'   \item \code{theta}: combined coefficient vector \code{c(alpha, beta)}.
#'   \item \code{se}: HAC standard errors of \code{theta}.
#'   \item \code{confint}: confidence intervals for \code{theta}.
#'   \item \code{D}: estimated hard regime indicators.
#'   \item \code{boundary}: estimated RFF coefficients defining the regime boundary.
#'   \item \code{smooth_fit}: complete first-stage smoothed-threshold fit.
#'   \item \code{refit}: complete second-stage hard-threshold refit.
#'   \item \code{rff}: random Fourier features and associated parameters.
#'   \item \code{tuning}: tuning parameters used in estimation.
#'   \item \code{n}: sample size.
#'   \item \code{call}: matched function call.
#' }
#'
#' @details
#' The fitted model is
#' \deqn{
#' Y_t = X_t^\top \alpha +
#'       X_t^\top \beta I\{f(Z_t) \ge U_t\} + \varepsilon_t,
#' }
#' where the nonlinear boundary \eqn{f} is approximated using random Fourier
#' features. A smooth surrogate of the regime indicator is used only for
#' first-stage boundary estimation. The final regression coefficients are
#' obtained by OLS refitting using the estimated hard regime indicators.
#'
#' @export
#' @importFrom stats qnorm
#'
#' @examples
#' ## Generate data
#' dat <- generate_AR(n = 500, seed = 123)
#'
#' ## Fit the model
#' fit <- TRM.RFF(
#'   Y = dat$Y,
#'   X = dat$X,
#'   Z = dat$Z,
#'   U = dat$U,
#'   rff_dim = 15,
#'   h = 0.025,
#'   sigma_h = 1
#' )
#'
#' ## Estimated coefficients and inference
#' fit$alpha
#' fit$beta
#' fit$se
#' fit$confint
#'
#' ## Estimated regime indicators
#' head(fit$D)

TRM.RFF <- function(Y, X, Z, U, rff_dim, h, sigma_h, seed = 1,
                    lambda = 1e-10, n_starts = 7L, conf_level = 0.95) {

  Y <- as.numeric(Y); U <- as.numeric(U)
  X <- as.matrix(X); Z <- as.matrix(Z)
  n <- length(Y)

  if (nrow(X) != n || nrow(Z) != n || length(U) != n)
    stop("Y, X, Z, and U must have the same number of observations.")
  if (h <= 0 || rff_dim <= 0)
    stop("h and rff_dim must be positive.")

  ## RFF and first-stage boundary estimation
  rff <- get_rffs(Z, rff_dim, sigma_h, seed = seed)
  smooth_fit <- fit_smooth_model_fast(
    Y, X, rff$Z, U, k = 1 / h,
    lambda = lambda, n_starts = n_starts
  )
  D_hat <- as.integer(smooth_fit$D)

  ## Hard-regime refit and HAC inference
  refit <- safe_lm_hac(Y, X, D_hat)
  alpha <- as.numeric(refit$alpha)
  beta  <- as.numeric(refit$beta)
  theta <- c(alpha, beta)
  se <- c(as.numeric(refit$se_alpha), as.numeric(refit$se_beta))

  zcrit <- qnorm(1 - (1 - conf_level) / 2)
  ci <- cbind(lower = theta - zcrit * se,
              upper = theta + zcrit * se)
  xnames <- colnames(X)
  if (is.null(xnames))
    xnames <- paste0("X", seq_len(ncol(X)))

  names(alpha) <- paste0("alpha_", xnames)
  names(beta)  <- paste0("beta_", xnames)
  theta <- c(alpha, beta)
  fit <- list(
    alpha = alpha, beta = beta, theta = theta,
    se = se, confint = ci, D = D_hat,
    boundary = smooth_fit$gamma,
    smooth_fit = smooth_fit, refit = refit, rff = rff,
    vcov = refit$vcov,
    tuning = list(rff_dim = rff_dim, h = h, sigma_h = sigma_h,
                  lambda = lambda, n_starts = n_starts),
    n = n, call = match.call()
  )

  class(fit) <- "TRM.RFF"
  fit
}

#' Extract coefficients from a TRM.RFF model
#'
#' @description
#' Extracts the estimated regression coefficients from a fitted
#' \code{TRM.RFF} model.
#'
#' @param object A fitted object of class \code{"TRM.RFF"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A numeric vector containing the estimated baseline coefficients
#' \code{alpha} followed by the regime-shift coefficients \code{beta}.
#'
#' @export
coef.TRM.RFF <- function(object, ...) {
  object$theta
}

#' Extract the HAC covariance matrix from a TRM.RFF model
#'
#' @description
#' Returns the heteroskedasticity- and autocorrelation-consistent (HAC)
#' covariance matrix of the refitted structural coefficient estimator.
#'
#' @param object A fitted object of class \code{"TRM.RFF"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A covariance matrix for the estimated coefficient vector
#' \code{c(alpha, beta)}.
#'
#' @export
vcov.TRM.RFF <- function(object, ...) {
  object$vcov
}

#' Confidence intervals for a TRM.RFF model
#'
#' @description
#' Computes normal-approximation confidence intervals for the structural
#' coefficients of a fitted \code{TRM.RFF} model using HAC standard errors.
#'
#' @param object A fitted object of class \code{"TRM.RFF"}.
#' @param parm A specification of which coefficients are to be given
#' confidence intervals. The default \code{NULL} returns intervals for all
#' coefficients.
#' @param level Confidence level. Default is \code{0.95}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A matrix containing the lower and upper confidence limits.
#'
#' @export
confint.TRM.RFF <- function(object, parm = NULL, level = 0.95, ...) {

  theta <- coef(object)
  se <- sqrt(diag(vcov(object)))

  if (!is.null(parm)) {
    theta <- theta[parm]
    se <- se[parm]
  }

  zcrit <- qnorm(1 - (1 - level) / 2)

  ci <- cbind(
    lower = theta - zcrit * se,
    upper = theta + zcrit * se
  )

  colnames(ci) <- paste0(
    format(100 * c((1 - level) / 2, 1 - (1 - level) / 2),
           trim = TRUE),
    "%"
  )

  ci
}


#' Summarize a fitted TRM.RFF model
#'
#' @description
#' Produces a summary of a fitted nonlinear threshold regression model,
#' including coefficient estimates, HAC standard errors, test statistics,
#' p-values, confidence intervals, and estimated regime sizes.
#'
#' @param object A fitted object of class \code{"TRM.RFF"}.
#' @param level Confidence level for coefficient intervals.
#' Default is \code{0.95}.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"summary.TRM.RFF"} containing:
#' \itemize{
#'   \item \code{coefficients}: coefficient table with estimates, HAC standard
#'   errors, z statistics, p-values, and confidence intervals.
#'   \item \code{regime}: numbers and proportions of observations assigned
#'   to the two estimated regimes.
#'   \item \code{tuning}: tuning parameters used in model fitting.
#'   \item \code{n}: sample size.
#'   \item \code{call}: original model call.
#' }
#'
#' @export
#' @importFrom stats pnorm printCoefmat
summary.TRM.RFF <- function(object, level = 0.95, ...) {

  theta <- object$theta
  se <- object$se

  zval <- theta / se
  pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)

  coef_table <- cbind(
    Estimate = theta,
    `HAC Std. Error` = se,
    `z value` = zval,
    `Pr(>|z|)` = pval
  )

  n0 <- sum(object$D == 0)
  n1 <- sum(object$D == 1)

  out <- list(
    coefficients = coef_table,
    regime = c(
      regime0 = n0,
      regime1 = n1,
      prop0 = n0 / object$n,
      prop1 = n1 / object$n
    ),
    tuning = object$tuning,
    n = object$n,
    call = object$call
  )

  class(out) <- "summary.TRM.RFF"
  out
}
#' Print a summary of a TRM.RFF model
#'
#' @description
#' Prints a formatted summary of a fitted nonlinear threshold regression
#' model with random Fourier features.
#'
#' @param x An object of class \code{"summary.TRM.RFF"}.
#' @param digits Number of significant digits to display.
#' @param ... Additional arguments.
#'
#' @return The summary object \code{x}, invisibly.
#'
#' @export
print.summary.TRM.RFF <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    ...
) {

  cat("\nNonlinear Threshold Regression with Random Fourier Features\n")
  cat("===========================================================\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients:\n")
  printCoefmat(
    x$coefficients,
    digits = digits,
    P.values = TRUE,
    has.Pvalue = TRUE
  )

  cat("\nEstimated regimes:\n")
  cat(
    "  Regime 0:",
    as.integer(x$regime["regime0"]),
    sprintf("(%.1f%%)", 100 * x$regime["prop0"]),
    "\n"
  )
  cat(
    "  Regime 1:",
    as.integer(x$regime["regime1"]),
    sprintf("(%.1f%%)", 100 * x$regime["prop1"]),
    "\n"
  )

  cat("\nTuning parameters:\n")
  cat("  RFF dimension :", x$tuning$rff_dim, "\n")
  cat("  h             :", x$tuning$h, "\n")
  cat("  sigma_h       :", x$tuning$sigma_h, "\n")
  cat("  lambda        :", x$tuning$lambda, "\n")
  cat("  n_starts      :", x$tuning$n_starts, "\n")

  cat("\nNumber of observations:", x$n, "\n")

  invisible(x)
}
inv_logit <- function(x) plogis(pmax(pmin(x, 35), -35))

get_rffs <- function(X, rff_dim = 5L, sigma = 0.07, seed = 1L) {
  X <- as.matrix(X)
  set.seed(seed)
  W <- matrix(rnorm(rff_dim * ncol(X)), nrow = rff_dim)
  b <- runif(rff_dim, 0, 2 * pi)
  Phi <- sqrt(2 / rff_dim) * cos(sigma * (W %*% t(X)) + b)
  list(Z = cbind(1, t(Phi)), W = W, b = b, sigma=sigma)
}



#' Predict from a fitted TRM.RFF model
#'
#' @description
#' Computes predictions from a fitted nonlinear threshold regression model
#' with random Fourier features. Predictions can be returned for the response,
#' the estimated regime classification, or the estimated nonlinear regime
#' boundary.
#'
#' @param object A fitted object of class \code{"TRM.RFF"}.
#' @param newdata An optional list containing new observations. It should
#' contain \code{X}, \code{Z}, and \code{U}. If \code{newdata = NULL},
#' predictions are computed for the observations used to fit the model.
#' @param type Character string specifying the type of prediction. One of
#' \code{"response"}, \code{"regime"}, or \code{"boundary"}.
#' Default is \code{"response"}.
#' @param ... Additional arguments (currently unused).
#'
#' @return Depending on \code{type}:
#' \itemize{
#'   \item \code{"response"}: a numeric vector of predicted response values.
#'   \item \code{"regime"}: an integer vector containing the estimated hard
#'   regime indicators.
#'   \item \code{"boundary"}: a numeric vector containing the estimated
#'   boundary values \eqn{\hat f(Z)}.
#' }
#'
#' @details
#' For a new observation with regressors \eqn{X}, state variables \eqn{Z},
#' and threshold variable \eqn{U}, the estimated regime is
#' \deqn{
#' \hat D = I\{\hat f(Z) \ge U\},
#' }
#' where \eqn{\hat f} is the nonlinear boundary estimated using random
#' Fourier features. The response prediction is
#' \deqn{
#' \hat Y =
#' X^\top \hat\alpha +
#' (X^\top \hat\beta)\hat D.
#' }
#'
#' The same random Fourier frequencies and phases used during model fitting
#' are reused when evaluating the boundary at new observations.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dat <- generate_AR(n = 500, seed = 123)
#'
#' fit <- TRM.RFF(
#'   Y = dat$Y, X = dat$X,
#'   Z = dat$Z, U = dat$U,
#'   rff_dim = 15,
#'   h = 0.025,
#'   sigma_h = 1
#' )
#'
#' ## Fitted responses
#' predict(fit)
#'
#' ## Estimated regimes
#' predict(fit, type = "regime")
#'
#' ## Estimated boundary values
#' predict(fit, type = "boundary")
#'
#' ## Prediction for new observations
#' newdat <- generate_AR(n = 100, seed = 456)
#'
#' predict(
#'   fit,
#'   newdata = list(
#'     X = newdat$X,
#'     Z = newdat$Z,
#'     U = newdat$U
#'   )
#' )
#' }
predict.TRM.RFF <- function(
    object,
    newdata = NULL,
    type = c("response", "regime", "boundary"),
    ...
) {

  type <- match.arg(type)

  ## Training data
  if (is.null(newdata)) {

    X <- object$smooth_fit$X
    Z <- object$smooth_fit$Z
    U <- object$smooth_fit$U

  } else {

    if (!is.list(newdata) ||
        !all(c("X", "Z", "U") %in% names(newdata)))
      stop("newdata must be a list containing X, Z, and U.")

    X <- as.matrix(newdata$X)
    Z <- as.matrix(newdata$Z)
    U <- as.numeric(newdata$U)

    if (nrow(X) != nrow(Z) || nrow(X) != length(U))
      stop("X, Z, and U must have the same number of observations.")
  }

  ## Evaluate RFF at Z using the original random features
  Psi <- transform_rffs(Z, object$rff)

  ## Estimated nonlinear boundary
  f_hat <- as.numeric(Psi %*% object$boundary)

  if (type == "boundary")
    return(f_hat)

  ## Hard regime classification
  D_hat <- as.integer(f_hat >= U)

  if (type == "regime")
    return(D_hat)

  ## Response prediction
  as.numeric(
    X %*% object$alpha +
      (X %*% object$beta) * D_hat
  )
}

transform_rffs <- function(X, rff) {

  X <- as.matrix(X)

  if (ncol(X) != ncol(rff$W))
    stop("X has an incompatible number of columns.")

  rff_dim <- nrow(rff$W)

  Phi <- sqrt(2 / rff_dim) *
    cos(
      rff$sigma * (rff$W %*% t(X)) + rff$b
    )

  cbind(1, t(Phi))
}

safe_lm <- function(Y, X, D, ridge = 1e-8) {
  W <- cbind(X, X * as.numeric(D))
  n <- nrow(W)
  cp <- crossprod(W)
  rhs <- crossprod(W, Y)

  # 求解回归系数
  inv_cp <- tryCatch(solve(cp), error = function(e) solve(cp + ridge * diag(ncol(W))))
  th <- as.numeric(inv_cp %*% rhs)

  # 计算残差与三明治协方差矩阵 V_hat = inv(W'W) * (W' diag(e^2) W) * inv(W'W)
  fitted <- as.numeric(W %*% th)
  resid <- Y - fitted
  meat <- crossprod(W * resid)
  vcov_mat <- inv_cp %*% meat %*% inv_cp
  se <- sqrt(pmax(diag(vcov_mat), 0))

  d <- ncol(X)
  list(alpha = th[seq_len(d)],
       beta  = th[d + seq_len(d)],
       se_alpha = se[seq_len(d)],
       se_beta  = se[d + seq_len(d)],
       fitted = fitted)
}

# ==============================================================================
# 1. 弱相依下的 HAC 重拟合与稳健方差估计函数
# ==============================================================================
safe_lm_hac <- function(Y, X, D, lag_order = NULL) {
  require(sandwich)
  T_obs <- length(Y)
  p_x <- ncol(X)

  # 构建分块设计矩阵 W = [X, X * D]
  W <- cbind(X, X * D)

  # 避免极少数区制极度不平衡或完全共线性导致的奇异值
  fit <- tryCatch({
    lm(Y ~ W - 1)
  }, error = function(e) NULL)

  if (is.null(fit) || any(is.na(coef(fit)))) {
    return(list(
      alpha     = rep(NA_real_, p_x),
      beta      = rep(NA_real_, p_x),
      se_alpha  = rep(NA_real_, p_x),
      se_beta   = rep(NA_real_, p_x),
      residuals = rep(NA_real_, T_obs)
    ))
  }

  # 确定 Newey-West HAC 滞后阶数（默认按理论公式 mT = floor(4 * (T/100)^(2/9))）
  if (is.null(lag_order)) {
    lag_order <- max(1L, floor(4 * (T_obs / 100)^(2 / 9)))
  }

  # 计算 Newey-West HAC 协方差矩阵
  V_hac <- tryCatch({
    sandwich::NeweyWest(fit, lag = lag_order, prewhite = FALSE, adjust = TRUE)
  }, error = function(e) {
    # 极端情况下若 HAC 计算失败，回退到普通稳健白噪声方差
    sandwich::vcovHC(fit, type = "HC1")
  })

  se_all <- sqrt(pmax(diag(V_hac), 0))

  list(
    alpha     = as.numeric(coef(fit)[1:p_x]),
    beta      = as.numeric(coef(fit)[(p_x + 1):(2 * p_x)]),
    se_alpha  = as.numeric(se_all[1:p_x]),
    se_beta   = as.numeric(se_all[(p_x + 1):(2 * p_x)]),
    residuals = as.numeric(residuals(fit)),
    vcov = V_hac
  )
}


#' @importFrom  stats optim
fit_smooth_model_fast <- function(Y, X, boundary_design, U,
                                  k = 40, lambda = 1e-10, lower = -10, upper = 10,
                                  maxit = 400L, n_starts = 5L) {
  X <- as.matrix(X); boundary_design <- as.matrix(boundary_design)
  d <- ncol(X); q <- ncol(boundary_design)

  # 联合目标函数与解析梯度
  obj_and_grad <- function(par) {
    alpha <- par[seq_len(d)]
    beta  <- par[d + seq_len(d)]
    gamma <- par[2 * d + seq_len(q)]

    eta <- k * (as.numeric(boundary_design %*% gamma) - U)
    eta_clamped <- pmax(pmin(eta, 35), -35)
    prob <- plogis(eta_clamped)

    X_alpha <- as.numeric(X %*% alpha)
    X_beta  <- as.numeric(X %*% beta)
    pred <- X_alpha + X_beta * prob
    resid <- Y - pred

    # 目标函数值
    val <- sum(resid^2) + lambda * sum(gamma^2)

    # 解析梯度
    # dval/dalpha = -2 * X' resid
    grad_alpha <- -2 * crossprod(X, resid)
    # dval/dbeta = -2 * (X * prob)' resid
    grad_beta  <- -2 * crossprod(X * prob, resid)
    # dval/dgamma = -2 * boundary_design' (resid * X_beta * prob * (1 - prob) * k) + 2 * lambda * gamma
    dprob_deta <- prob * (1 - prob)
    chain_gamma <- resid * X_beta * dprob_deta * k
    grad_gamma <- -2 * crossprod(boundary_design, chain_gamma) + 2 * lambda * gamma

    attr(val, "gradient") <- c(as.numeric(grad_alpha), as.numeric(grad_beta), as.numeric(grad_gamma))
    val
  }

  objective <- function(par) {
    res <- obj_and_grad(par)
    res
  }

  gradient <- function(par) {
    attr(obj_and_grad(par), "gradient")
  }

  # 生成初始点 (可适当将 n_starts 降为 3~5 个即可稳定捕获好解)
  gamma_starts <- vector("list", n_starts)
  gamma_starts[[1]] <- rep(0, q)
  if (n_starts > 1L) {
    for (s in 2:n_starts) {
      gamma_starts[[s]] <- 0.75 * sin(seq_len(q) * (s - 1))
    }
  }

  starts <- lapply(gamma_starts, function(gamma_init) {
    make_data_driven_init(Y = Y, X = X, boundary_design = boundary_design,
                          U = U, gamma_init = gamma_init, k = k)
  })

  opts <- lapply(starts, function(start) {
    tryCatch(
      optim(start, fn = objective, gr = gradient, method = "L-BFGS-B",
            lower = rep(lower, length(start)), upper = rep(upper, length(start)),
            control = list(maxit = maxit, factr = 1e9, pgtol = 1e-5)),
      error = function(e) NULL
    )
  })

  ok <- vapply(opts, function(z) !is.null(z) && is.finite(z$value), logical(1))
  if (!any(ok)) stop("All L-BFGS-B starts failed.")

  opt <- opts[[which.min(vapply(opts, function(z) if (is.null(z)) Inf else z$value, numeric(1)))]]
  gamma <- opt$par[2 * d + seq_len(q)]
  Dhat <- as.numeric(boundary_design %*% gamma >= U)

  list(alpha = opt$par[seq_len(d)], beta = opt$par[d + seq_len(d)],
       gamma = gamma, D = Dhat, refit = safe_lm(Y, X, Dhat),
       convergence = opt$convergence, value = opt$value, counts = opt$counts,
       message = if (is.null(opt$message)) "" else opt$message)
}

make_data_driven_init <- function(Y, X, boundary_design, U,
                                  gamma_init, k = 40,
                                  ridge = 1e-7) {
  gamma_init <- as.numeric(gamma_init)

  prob_init <- plogis(
    pmax(
      pmin(
        k * (as.numeric(boundary_design %*% gamma_init) - U),
        35
      ),
      -35
    )
  )

  W_init <- cbind(X, X * prob_init)

  theta_init <- tryCatch(
    qr.solve(crossprod(W_init), crossprod(W_init, Y)),
    error = function(e) {
      solve(
        crossprod(W_init) + ridge * diag(ncol(W_init)),
        crossprod(W_init, Y)
      )
    }
  )

  c(as.numeric(theta_init), gamma_init)
}


