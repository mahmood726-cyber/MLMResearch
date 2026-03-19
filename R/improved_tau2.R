#' Improved Heterogeneity (tau²) Estimation
#'
#' Better methods for estimating between-study heterogeneity in multilevel
#' meta-analysis, especially for small samples and sparse data.
#'
#' @name improved_tau2
NULL

#' Improved Tau² Estimator
#'
#' Computes tau² using a combination of estimators with small-sample
#' corrections.
#'
#' @param yi Vector of effect sizes
#' @param vi Vector of sampling variances
#' @param method Primary estimation method ("REML", "DL", "SJ", "PM")
#' @return List with tau² estimate and method details
#' @export
improved_tau2_estimator <- function(yi, vi, method = "REML") {

  n <- length(yi)

  # Fit standard rma model for baseline
  fit <- tryCatch({
    metafor::rma.uni(yi = yi, vi = vi, method = method)
  }, error = function(e) NULL)

  if (is.null(fit)) {
    return(list(
      tau2 = NA,
      se_tau2 = NA,
      method = method,
      success = FALSE,
      message = "Model fitting failed"
    ))
  }

  tau2_base <- fit$tau2

  # Apply small-sample correction if n < 20
  if (n < 20 && tau2_base > 0) {
    # Use Sidik-Jonkman estimator as alternative
    fit_sj <- tryCatch({
      metafor::rma.uni(yi = yi, vi = vi, method = "SJ")
    }, error = function(e) NULL)

    if (!is.null(fit_sj)) {
      tau2_sj <- fit_sj$tau2

      # Average the estimators (simple ensemble)
      tau2_improved <- (tau2_base + tau2_sj) / 2

      # Use Paule-Mandel as third check
      fit_pm <- tryCatch({
        metafor::rma.uni(yi = yi, vi = vi, method = "PM")
      }, error = function(e) NULL)

      if (!is.null(fit_pm)) {
        tau2_pm <- fit_pm$tau2
        # Median of three estimators (more robust)
        tau2_improved <- median(c(tau2_base, tau2_sj, tau2_pm))
      }

      return(list(
        tau2 = tau2_improved,
        se_tau2 = fit$se.tau2,
        method = paste0(method, "-improved"),
        tau2_base = tau2_base,
        tau2_sj = tau2_sj,
        tau2_pm = if (!is.null(fit_pm)) fit_pm$tau2 else NA,
        success = TRUE,
        n_studies = n,
        message = "Small-sample correction applied"
      ))
    }
  }

  # Standard result for adequate samples
  list(
    tau2 = tau2_base,
    se_tau2 = fit$se.tau2,
    method = method,
    success = TRUE,
    n_studies = n,
    message = "Standard estimation"
  )
}

#' Profile Likelihood Confidence Interval for Tau²
#'
#' Computes profile likelihood CI for tau², which is more accurate than
#' Wald intervals for heterogeneity parameters.
#'
#' @param yi Vector of effect sizes
#' @param vi Vector of sampling variances
#' @param level Confidence level (default: 0.95)
#' @return List with tau² estimate and CI
#' @export
profile_tau2_ci <- function(yi, vi, level = 0.95) {

  fit <- tryCatch({
    metafor::rma.uni(yi = yi, vi = vi, method = "REML")
  }, error = function(e) NULL)

  if (is.null(fit)) {
    return(list(
      tau2 = NA,
      ci_lb = NA,
      ci_ub = NA,
      success = FALSE
    ))
  }

  # Get profile likelihood CI
  tryCatch({
    prof <- metafor::profile(fit)
    ci <- metafor::confint(prof, level = level)

    list(
      tau2 = fit$tau2,
      ci_lb = ci$tau2[1],
      ci_ub = ci$tau2[2],
      method = "profile",
      success = TRUE
    )
  }, error = function(e) {
    # Fallback to Q-profile method
    list(
      tau2 = fit$tau2,
      ci_lb = fit$ci.lb.tau2,
      ci_ub = fit$ci.ub.tau2,
      method = "q-profile",
      success = TRUE
    )
  })
}

#' QE-based Tau² Estimator
#'
#' Uses the Q-statistic to estimate tau², which can be more stable
#' in small samples.
#'
#' @param yi Vector of effect sizes
#' @param vi Vector of sampling variances
#' @return List with tau² estimate
#' @export
qe_tau2_estimator <- function(yi, vi) {

  n <- length(yi)

  # Compute weighted mean
  wi <- 1 / vi
  yi_bar <- sum(wi * yi) / sum(wi)

  # Q statistic
  Q <- sum(wi * (yi - yi_bar)^2)

  # Degrees of freedom
  df <- n - 1

  # DerSimonian-Laird estimator
  if (Q > df) {
    tau2 <- (Q - df) / (sum(wi) - sum(wi^2) / sum(wi))
  } else {
    tau2 <- 0
  }

  # Improved DL with truncation at 0
  tau2 <- max(0, tau2)

  list(
    tau2 = tau2,
    Q = Q,
    df = df,
    method = "QE-based",
    success = TRUE
  )
}

#' Multilevel Tau² Estimation
#'
#' Estimates tau² at multiple levels in a multilevel model.
#'
#' @param df Data frame with TE, seTE, and clustering variables
#' @param cluster_var Primary clustering variable
#' @return List with multilevel tau² estimates
#' @export
multilevel_tau2 <- function(df, cluster_var = "review_id") {

  fit <- tryCatch({
    V <- diag(df$seTE^2, nrow = nrow(df))
    metafor::rma.mv(
      yi = df$TE,
      V = V,
      random = list(~ 1 | df[[cluster_var]], ~ 1 | df$study_id),
      method = "REML"
    )
  }, error = function(e) NULL)

  if (is.null(fit)) {
    return(list(
      tau2 = NA,
      sigma2 = NA,
      success = FALSE,
      message = "Multilevel model fitting failed"
    ))
  }

  # Extract variance components
  tau2 <- fit$tau2
  sigma2 <- fit$sigma2

  # Calculate I² for multilevel model
  total_var <- tau2 + mean(sigma2) + mean(df$seTE^2)
  i2_total <- (tau2 + mean(sigma2)) / total_var * 100

  # I² at each level
  i2_level1 <- tau2 / total_var * 100
  i2_level2 <- mean(sigma2) / total_var * 100

  list(
    tau2 = tau2,
    sigma2 = sigma2,
    i2_total = i2_total,
    i2_level1 = i2_level1,
    i2_level2 = i2_level2,
    success = TRUE,
    fit = fit
  )
}

#' Bootstrap Tau² Confidence Interval
#'
#' Uses parametric bootstrap to estimate CI for tau².
#'
#' @param yi Vector of effect sizes
#' @param vi Vector of sampling variances
#' @param n_bootstrap Number of bootstrap samples (default: 1000)
#' @param seed Random seed
#' @param level Confidence level (default: 0.95)
#' @return List with tau² estimate and bootstrap CI
#' @export
bootstrap_tau2_ci <- function(yi, vi, n_bootstrap = 1000, seed = 123, level = 0.95) {

  set.seed(seed)

  # Fit base model
  fit <- tryCatch({
    metafor::rma.uni(yi = yi, vi = vi, method = "REML")
  }, error = function(e) NULL)

  if (is.null(fit)) {
    return(list(
      tau2 = NA,
      ci_lb = NA,
      ci_ub = NA,
      success = FALSE,
      message = "Base model fitting failed"
    ))
  }

  tau2_base <- fit$tau2
  n <- length(yi)

  # Bootstrap iterations
  tau2_bootstrap <- numeric(n_bootstrap)

  for (i in seq_len(n_bootstrap)) {
    # Generate bootstrap sample
    yi_boot <- yi + rnorm(n, 0, sqrt(vi))

    # Fit model to bootstrap sample
    fit_boot <- tryCatch({
      metafor::rma.uni(yi = yi_boot, vi = vi, method = "REML")
    }, error = function(e) NULL)

    if (!is.null(fit_boot)) {
      tau2_bootstrap[i] <- fit_boot$tau2
    } else {
      tau2_bootstrap[i] <- NA
    }
  }

  # Remove failed iterations
  tau2_valid <- tau2_bootstrap[!is.na(tau2_bootstrap)]

  if (length(tau2_valid) < n_bootstrap / 2) {
    return(list(
      tau2 = tau2_base,
      ci_lb = NA,
      ci_ub = NA,
      success = FALSE,
      message = "Too many bootstrap failures"
    ))
  }

  # Calculate bootstrap CI
  alpha <- 1 - level
  ci_lb <- quantile(tau2_valid, alpha / 2)
  ci_ub <- quantile(tau2_valid, 1 - alpha / 2)

  # Bias-corrected estimate
  tau2_bc <- 2 * tau2_base - mean(tau2_valid)

  list(
    tau2 = tau2_base,
    tau2_bias_corrected = tau2_bc,
    ci_lb = ci_lb,
    ci_ub = ci_ub,
    ci_method = "percentile",
    n_bootstrap = n_bootstrap,
    n_valid = length(tau2_valid),
    success = TRUE
  )
}

#' Compare Tau² Estimators
#'
#' Compares different tau² estimation methods on the same data.
#'
#' @param df Data frame with TE and seTE
#' @param methods Methods to compare (default: all available)
#' @return Data frame with comparison
#' @export
compare_tau2_methods <- function(df,
                                  methods = c("REML", "DL", "SJ", "PM", "HE")) {

  yi <- df$TE
  vi <- df$seTE^2

  results <- list()

  for (m in methods) {
    fit <- tryCatch({
      metafor::rma.uni(yi = yi, vi = vi, method = m)
    }, error = function(e) NULL)

    if (!is.null(fit)) {
      results[[m]] <- list(
        method = m,
        tau2 = fit$tau2,
        se = fit$se,
        i2 = fit$I2
      )
    } else {
      results[[m]] <- list(
        method = m,
        tau2 = NA,
        se = NA,
        i2 = NA,
        error = TRUE
      )
    }
  }

  # Convert to data frame
  comparison <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      method = x$method,
      tau2 = round(x$tau2, 4),
      i2 = round(x$i2, 2),
      error = !is.null(x$error),
      stringsAsFactors = FALSE
    )
  }))

  # Add range and coefficient of variation
  valid_tau2 <- comparison$tau2[!is.na(comparison$tau2) & comparison$tau2 > 0]
  tau2_cv <- if (length(valid_tau2) > 1) sd(valid_tau2) / mean(valid_tau2) else NA

  list(
    comparison = comparison,
    tau2_range = range(valid_tau2),
    tau2_cv = tau2_cv,
    n_studies = length(yi)
  )
}

#' Get Recommended Tau² Method
#'
#' Returns recommendation for tau² estimation based on sample size.
#'
#' @param n_studies Number of studies
#' @return Recommended method and reasoning
#' @export
get_tau2_recommendation <- function(n_studies) {

  if (n_studies >= 40) {
    list(
      method = "REML",
      confidence = "high",
      reasoning = "Large sample - REML is optimal"
    )
  } else if (n_studies >= 20) {
    list(
      method = "REML or DL",
      confidence = "medium",
      reasoning = "Moderate sample - REML preferred, DL acceptable"
    )
  } else if (n_studies >= 10) {
    list(
      method = "SJ or improved ensemble",
      confidence = "medium",
      reasoning = "Small sample - consider Sidik-Jonkman or ensemble"
    )
  } else {
    list(
      method = "Use caution - consider Bayesian",
      confidence = "low",
      reasoning = "Very small sample - tau² estimates are unreliable"
    )
  }
}
