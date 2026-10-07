
#' Sup-Wald test for regime effects in nonlinear threshold regression
#'
#' @description
#' Performs a Sup-Wald test for the existence of regime effects in a
#' nonlinear threshold regression model with a regime boundary approximated
#' by random Fourier features. Critical values and p-values are obtained
#' using a dependent wild bootstrap.
#'
#' @param Y A numeric response vector of length \code{n}.
#' @param X An \code{n x d} matrix of regression covariates.
#' @param Z An \code{n x p} matrix of state variables entering the nonlinear
#' regime boundary.
#' @param U A numeric threshold variable of length \code{n}.
#' @param rff_dim Number of random Fourier features. Default is \code{20}.
#' @param sigma_h Scale parameter used in the random Fourier feature mapping.
#' Default is \code{0.5}.
#' @param n_grid Number of candidate boundary functions used to approximate
#' the supremum. Default is \code{25}.
#' @param B Number of dependent wild bootstrap replications.
#' Default is \code{399}.
#' @param hac_lag Nonnegative integer specifying the truncation lag used in
#' the HAC covariance estimator. Default is \code{2}.
#' @param dwb_bandwidth Positive bandwidth parameter for the dependent wild
#' bootstrap multipliers. Default is \code{hac_lag + 1}.
#' @param trim Trimming proportion used to exclude candidate regime
#' partitions with too few observations in either regime. Default is
#' \code{0.05}.
#' @param rff_seed Random seed used to generate the random Fourier features.
#' Default is \code{1}.
#' @param grid_seed Random seed used to generate the candidate boundary grid.
#' Default is \code{999}.
#' @param bootstrap_seed Optional random seed for the bootstrap. If
#' \code{NULL}, the current random-number state is used.
#' @param parallel Logical indicating whether the bootstrap replications
#' should be computed in parallel. Default is \code{FALSE}.
#' @param n_cores Number of parallel workers. If \code{NULL}, the number of
#' available cores minus one is used. This argument is ignored when
#' \code{parallel = FALSE}.
#' @param keep_boot Logical indicating whether bootstrap statistics should
#' be included in the returned object. Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"supWald.RFF"} containing:
#' \itemize{
#'   \item \code{statistic}: the observed Sup-Wald statistic.
#'   \item \code{p.value}: the dependent-wild-bootstrap p-value.
#'   \item \code{parameter}: the number of tested regime-shift coefficients.
#'   \item \code{B}: the number of bootstrap replications.
#'   \item \code{valid_grid}: the number of valid candidate regime partitions.
#'   \item \code{max_index}: the candidate-grid index attaining the observed
#'   supremum.
#'   \item \code{wald}: the candidate-specific Wald statistics.
#'   \item \code{boot_stats}: bootstrap Sup-Wald statistics when
#'   \code{keep_boot = TRUE}.
#'   \item \code{parallel}: whether parallel computation was used.
#'   \item \code{n_cores}: number of workers used for the bootstrap.
#'   \item \code{call}: the matched function call.
#' }
#'
#' @details
#' The null hypothesis is
#' \deqn{
#' H_0: \beta_0 = 0,
#' }
#' corresponding to the absence of regime effects. Under the null, the
#' nonlinear regime boundary is not identified. The test therefore computes
#' Wald statistics over a collection of candidate regime boundaries and uses
#' their maximum as the test statistic.
#'
#' Candidate boundaries are constructed using random Fourier features.
#' For each candidate boundary, observations are classified into two regimes,
#' and a HAC covariance estimator is used to construct the corresponding
#' Wald statistic for the regime-shift coefficients.
#'
#' Because the boundary is unidentified under the null, the null distribution
#' of the Sup-Wald statistic is approximated using a dependent wild bootstrap.
#' Bootstrap samples are generated from the null-restricted regression
#' \deqn{
#' Y_t = X_t^\top \alpha_0 + \varepsilon_t.
#' }
#'
#' Parallel computation is implemented using a multisession backend and is
#' therefore supported on both Windows and Unix-like operating systems.
#'
#' @export
#'
#' @examples
#' dat <- generate_AR(n = 500, seed = 123)
#'
#' ## Sequential computation
#' test <- supWald.RFF(
#'   Y = dat$Y,
#'   X = dat$X,
#'   Z = dat$Z,
#'   U = dat$U,
#'   rff_dim = 15,
#'   sigma_h = 1,
#'   n_grid = 25,
#'   B = 99,
#'   parallel = FALSE
#' )
#'
#' test
#'
#' \dontrun{
#' ## Parallel computation using 4 workers
#' test_parallel <- supWald.RFF(
#'   Y = dat$Y,
#'   X = dat$X,
#'   Z = dat$Z,
#'   U = dat$U,
#'   rff_dim = 15,
#'   sigma_h = 1,
#'   n_grid = 25,
#'   B = 399,
#'   parallel = TRUE,
#'   n_cores = 4,
#'   bootstrap_seed = 123
#' )
#'
#' test_parallel
#' }
supWald.RFF <- function(
    Y,
    X,
    Z,
    U,
    rff_dim = 20L,
    sigma_h = 0.5,
    n_grid = 25L,
    B = 399L,
    hac_lag = 2L,
    dwb_bandwidth = hac_lag + 1,
    trim = 0.05,
    rff_seed = 1L,
    grid_seed = 999L,
    bootstrap_seed = NULL,
    parallel = FALSE,
    n_cores = NULL,
    keep_boot = FALSE,
    ...
) {

  ## ------------------------------------------------------------
  ## Input checks
  ## ------------------------------------------------------------

  Y <- as.numeric(Y)
  U <- as.numeric(U)
  X <- as.matrix(X)
  Z <- as.matrix(Z)

  n <- length(Y)
  p_x <- ncol(X)

  if (nrow(X) != n ||
      nrow(Z) != n ||
      length(U) != n) {
    stop("Y, X, Z, and U must have the same number of observations.")
  }

  if (rff_dim <= 0L)
    stop("rff_dim must be positive.")

  if (n_grid <= 0L)
    stop("n_grid must be positive.")

  if (B <= 0L)
    stop("B must be positive.")

  if (hac_lag < 0L)
    stop("hac_lag must be nonnegative.")

  if (trim < 0 || trim >= 0.5)
    stop("trim must lie in [0, 0.5).")

  if (dwb_bandwidth <= 0)
    stop("dwb_bandwidth must be positive.")

  if (!is.logical(parallel) || length(parallel) != 1L)
    stop("parallel must be TRUE or FALSE.")

  ## ------------------------------------------------------------
  ## Number of workers
  ## ------------------------------------------------------------

  if (parallel) {

    if (!requireNamespace("future", quietly = TRUE) ||
        !requireNamespace("future.apply", quietly = TRUE)) {
      stop(
        "Packages 'future' and 'future.apply' are required ",
        "for parallel computation."
      )
    }

    if (is.null(n_cores)) {

      available <- future::availableCores()

      n_cores <- max(1L, as.integer(available) - 1L)

    } else {

      if (length(n_cores) != 1L ||
          !is.finite(n_cores) ||
          n_cores < 1L) {
        stop("n_cores must be a positive integer.")
      }

      n_cores <- as.integer(n_cores)
    }

  } else {

    n_cores <- 1L
  }

  ## ------------------------------------------------------------
  ## Random Fourier features
  ## ------------------------------------------------------------

  rff <- get_rffs(
    Z,
    rff_dim = rff_dim,
    sigma = sigma_h,
    seed = rff_seed
  )

  Psi <- rff$Z
  L_features <- ncol(Psi)

  ## ------------------------------------------------------------
  ## Candidate boundary grid
  ## ------------------------------------------------------------

  set.seed(grid_seed)

  candidate_gamma <- matrix(
    stats::rnorm(L_features * n_grid),
    nrow = L_features,
    ncol = n_grid
  )

  gamma_norm <- sqrt(colSums(candidate_gamma^2))

  candidate_gamma <- sweep(
    candidate_gamma,
    2L,
    gamma_norm,
    "/"
  )

  candidate_D_list <- lapply(
    seq_len(n_grid),
    function(k) {

      f_k <- as.numeric(
        Psi %*% candidate_gamma[, k]
      )

      as.integer(f_k >= U)
    }
  )

  ## ------------------------------------------------------------
  ## Precompute candidate regressions
  ## ------------------------------------------------------------

  proj_info <- prepare_candidate_projections(
    X,
    candidate_D_list,
    trim = trim
  )

  valid_grid <- sum(
    vapply(
      proj_info,
      function(x) isTRUE(x$valid),
      logical(1)
    )
  )

  if (valid_grid == 0L)
    stop("No valid candidate regime partitions were found.")

  ## ------------------------------------------------------------
  ## Observed Sup-Wald statistic
  ## ------------------------------------------------------------

  observed <- calc_sup_wald_fast(
    Y = Y,
    proj_info = proj_info,
    n = n,
    hac_lag = hac_lag
  )

  ## ------------------------------------------------------------
  ## Null-restricted model
  ## ------------------------------------------------------------

  fit_null <- stats::lm.fit(X, Y)

  res_null <- fit_null$residuals
  pred_null <- fit_null$fitted.values

  ## ------------------------------------------------------------
  ## Bootstrap function
  ## ------------------------------------------------------------

  bootstrap_one <- function(b) {

    eta <- generate_dwb_multipliers(
      n,
      bandwidth = dwb_bandwidth
    )

    Y_star <- pred_null + res_null * eta

    calc_sup_wald_fast(
      Y = Y_star,
      proj_info = proj_info,
      n = n,
      hac_lag = hac_lag
    )$statistic
  }

  ## ------------------------------------------------------------
  ## Dependent wild bootstrap
  ## ------------------------------------------------------------

  if (!is.null(bootstrap_seed))
    set.seed(bootstrap_seed)

  if (parallel && n_cores > 1L) {

    ## Save current future strategy
    old_plan <- future::plan()

    ## Restore it when leaving this function
    on.exit(
      future::plan(old_plan),
      add = TRUE
    )

    ## Multisession works on Windows, Linux, and macOS
    future::plan(
      future::multisession,
      workers = n_cores
    )

    boot_stats <- unlist(
      future.apply::future_lapply(
        seq_len(B),
        bootstrap_one,
        future.seed = TRUE
      ),
      use.names = FALSE
    )

  } else {

    boot_stats <- unlist(
      lapply(
        seq_len(B),
        bootstrap_one
      ),
      use.names = FALSE
    )
  }

  ## ------------------------------------------------------------
  ## Finite-bootstrap p-value
  ## ------------------------------------------------------------

  p_value <- (
    1 + sum(boot_stats >= observed$statistic)
  ) / (B + 1)

  ## ------------------------------------------------------------
  ## Output
  ## ------------------------------------------------------------

  out <- list(
    statistic = observed$statistic,
    p.value = p_value,
    parameter = p_x,
    B = B,
    valid_grid = valid_grid,
    max_index = observed$index,
    wald = observed$wald,
    candidate_gamma = candidate_gamma,
    rff = rff,
    boot_stats = if (keep_boot) boot_stats else NULL,
    parallel = parallel && n_cores > 1L,
    n_cores = n_cores,
    method = paste(
      "Sup-Wald test for regime effects",
      "with dependent wild bootstrap"
    ),
    call = match.call()
  )

  class(out) <- "supWald.RFF"

  out
}



#' Print a Sup-Wald test
#'
#' @param x An object of class \code{"supWald.RFF"}.
#' @param digits Number of significant digits to print.
#' @param ... Additional arguments (currently unused).
#'
#' @return The object \code{x}, invisibly.
#'
#' @export
print.supWald.RFF <- function(
    x,
    digits = max(3L, getOption("digits") - 3L),
    ...
) {

  cat("\n")
  cat(x$method, "\n")
  cat("====================================================\n\n")

  cat("Call:\n")
  print(x$call)

  cat("\nNull hypothesis:\n")
  cat("  No regime effects (beta = 0)\n\n")

  cat(
    "Sup-Wald statistic:",
    format(x$statistic, digits = digits),
    "\n"
  )

  cat(
    "Bootstrap p-value:",
    format.pval(x$p.value, digits = digits),
    "\n"
  )

  cat(
    "Bootstrap replications:",
    x$B,
    "\n"
  )

  cat(
    "Valid candidate boundaries:",
    x$valid_grid,
    "\n"
  )

  invisible(x)
}



# ==============================================================================
# 0. Dependent Wild Bootstrap (DWB) multiplier generating function (Shao, 2010).
# ==============================================================================
generate_dwb_multipliers <- function(T_obs, bandwidth = 3.0) {
  # Generate dependent normal/transformed multipliers based on a Bartlett/Parzen-type autocovariance matrix.
  # Ensure E*[eta_t * eta_s] = k((t - s) / bandwidth)
  t_idx <- 1:T_obs
  dist_mat <- abs(outer(t_idx, t_idx, "-"))

  # Bartlett (Newey-West) kernel-based dependence structure.
  kernel_weights <- pmax(0, 1 - dist_mat / bandwidth)

  # Perform Cholesky decomposition on the kernel matrix (adding a tiny ridge term to prevent numerical singularity).
  R <- chol(kernel_weights + diag(1e-8, T_obs))

  # Generate multipliers with the corresponding autocovariance.
  as.numeric(crossprod(R, rnorm(T_obs)))
}

# ==============================================================================
# 1. Pre-extract the projection matrix information of the candidate grid (to greatly accelerate the bootstrap).
# ==============================================================================
#' Prepare candidate regime projections
#'
#' @keywords internal
prepare_candidate_projections <- function(
    X,
    candidate_D_list,
    trim = 0.05
) {

  n <- nrow(X)
  p <- ncol(X)

  proj_info <- vector("list", length(candidate_D_list))

  for (k in seq_along(candidate_D_list)) {

    D_k <- candidate_D_list[[k]]
    p_regime <- mean(D_k)

    if (p_regime < trim || p_regime > 1 - trim) {
      proj_info[[k]] <- list(valid = FALSE)
      next
    }

    W_k <- cbind(X, X * D_k)

    A_mat <- crossprod(W_k) / n

    A_inv <- tryCatch(
      solve(A_mat),
      error = function(e) NULL
    )

    if (is.null(A_inv)) {
      proj_info[[k]] <- list(valid = FALSE)
      next
    }

    P_W <- (A_inv %*% t(W_k)) / n

    proj_info[[k]] <- list(
      valid = TRUE,
      W_k = W_k,
      A_inv = A_inv,
      P_W = P_W,
      beta_idx = (p + 1L):(2L * p)
    )
  }

  proj_info
}

# ==============================================================================
# 2. Quickly compute the Sup-Wald statistic for a given Y.
# ==============================================================================
#' Compute the Sup-Wald statistic
#'
#' @keywords internal
calc_sup_wald_fast <- function(
    Y,
    proj_info,
    n,
    hac_lag = 2L
) {

  wald_vals <- rep(NA_real_, length(proj_info))

  for (k in seq_along(proj_info)) {

    info <- proj_info[[k]]

    if (!info$valid)
      next

    ## OLS coefficients
    coef_hat <- as.numeric(info$P_W %*% Y)
    beta_hat <- coef_hat[info$beta_idx]

    ## Residuals and scores
    res <- Y - as.numeric(info$W_k %*% coef_hat)
    scores <- info$W_k * res

    ## HAC long-run covariance
    Omega_hat <- crossprod(scores) / n

    if (hac_lag > 0L) {

      for (l in seq_len(hac_lag)) {

        w_l <- 1 - l / (hac_lag + 1)

        Gamma_l <- crossprod(
          scores[(l + 1L):n, , drop = FALSE],
          scores[seq_len(n - l), , drop = FALSE]
        ) / n

        Omega_hat <-
          Omega_hat +
          w_l * (Gamma_l + t(Gamma_l))
      }
    }

    ## Covariance of regression coefficients
    V_hac <-
      info$A_inv %*%
      Omega_hat %*%
      info$A_inv / n

    V_beta <- V_hac[
      info$beta_idx,
      info$beta_idx,
      drop = FALSE
    ]

    V_beta_inv <- tryCatch(
      solve(V_beta),
      error = function(e) NULL
    )

    if (is.null(V_beta_inv))
      next

    wald_vals[k] <- as.numeric(
      crossprod(beta_hat, V_beta_inv %*% beta_hat)
    )
  }

  if (all(is.na(wald_vals)))
    stop("No valid candidate regime partitions were found.")

  list(
    statistic = max(wald_vals, na.rm = TRUE),
    wald = wald_vals,
    index = which.max(replace(wald_vals, is.na(wald_vals), -Inf))
  )
}
