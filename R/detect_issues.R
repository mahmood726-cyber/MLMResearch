#' Issue Detection Functions for Multilevel Meta-Analysis
#'
#' Functions to detect and categorize issues that arise when fitting
#' multilevel meta-analysis models to Cochrane review data.
#'
#' @name detect_issues
NULL

#' Check for Convergence Problems
#'
#' @param fit A fitted metafor model (rma.mv) or test result
#' @return List with convergence diagnostics
#' @export
check_convergence_problems <- function(fit) {

  issues <- list()

  # If fit is a test result list, extract the actual fit
  if (is.list(fit) && !is.null(fit$fit)) {
    fit <- fit$fit
  }

  # Check if fit exists
  if (is.null(fit) || !inherits(fit, "rma.mv")) {
    return(list(
      has_issues = TRUE,
      issue_type = "no_fit",
      severity = "critical",
      message = "Model fitting failed - no fit object returned"
    ))
  }

  # Check for NaN or infinite estimates
  if (any(!is.finite(fit$beta))) {
    issues$nan_estimate <- TRUE
  }

  # Check for negative variance components
  if (!is.null(fit$tau2) && fit$tau2 < 0) {
    issues$negative_tau2 <- TRUE
  }

  if (!is.null(fit$sigma2) && any(fit$sigma2 < 0)) {
    issues$negative_sigma2 <- TRUE
  }

  # Check for near-zero variance
  if (!is.null(fit$tau2) && fit$tau2 < 0.0001) {
    issues$near_zero_tau2 <- TRUE
  }

  # Check for extreme SE
  if (any(fit$se > 10)) {
    issues$large_se <- TRUE
  }

  # Check for warnings
  if (length(fit$warnings) > 0) {
    issues$warnings <- fit$warnings
  }

  has_issues <- length(issues) > 0

  list(
    has_issues = has_issues,
    issues = issues,
    issue_type = if (has_issues) names(issues)[1] else "none",
    severity = ifelse(
      has_issues && any(names(issues) %in% c("nan_estimate", "negative_tau2")),
      "high", "low"
    ),
    message = if (has_issues) paste("Issues detected:", paste(names(issues), collapse = ", ")) else "No convergence issues"
  )
}

#' Check for Small Sample Issues
#'
#' @param df Data frame with the analysis data
#' @param fit Fitted model or test result
#' @return List with small sample diagnostics
#' @export
check_small_sample_issues <- function(df, fit = NULL) {

  n_effects <- nrow(df)
  n_studies <- length(unique(df$study_id))
  n_clusters <- length(unique(df$review_id))

  issues <- list()

  # Check if too few studies
  if (n_studies < 5) {
    issues$few_studies <- TRUE
  }

  # Check if too few clusters for RVE
  if (n_clusters < 10) {
    issues$few_clusters_rve <- TRUE
  }
  if (n_clusters < 20) {
    issues$few_clusters_rve_recommended <- TRUE
  }

  # Check for sparse data within studies
  study_counts <- table(df$study_id)
  if (any(study_counts == 1)) {
    issues$single_effect_studies <- TRUE
  }

  # Check if fit was provided and has issues
  if (!is.null(fit)) {
    conv_check <- check_convergence_problems(fit)
    if (conv_check$has_issues) {
      issues$convergence_issues <- TRUE
    }
  }

  has_issues <- length(issues) > 0

  list(
    has_issues = has_issues,
    issues = issues,
    n_effects = n_effects,
    n_studies = n_studies,
    n_clusters = n_clusters,
    issue_type = if (has_issues) "small_sample" else "none",
    severity = ifelse(
      n_studies < 5 || n_clusters < 10,
      "high", "medium"
    ),
    message = if (has_issues) {
      paste("Small sample issues:", paste(names(issues), collapse = ", "))
    } else {
      "Sample size adequate"
    }
  )
}

#' Check for Heterogeneity Problems
#'
#' @param fit Fitted model or test result
#' @param df Data frame (optional, for additional checks)
#' @return List with heterogeneity diagnostics
#' @export
check_heterogeneity_problems <- function(fit, df = NULL) {

  # If fit is a test result list, extract the actual fit
  if (is.list(fit) && !is.null(fit$fit)) {
    fit <- fit$fit
  }

  issues <- list()

  # Check if fit exists
  if (is.null(fit) || !inherits(fit, "rma.mv")) {
    return(list(
      has_issues = TRUE,
      issue_type = "no_fit",
      severity = "critical",
      message = "Cannot check heterogeneity - no fit object"
    ))
  }

  # Extract heterogeneity statistics
  tau2 <- fit$tau2
  i2 <- fit$I2

  # Check for zero heterogeneity
  if (!is.null(tau2) && tau2 < 0.0001) {
    issues$zero_heterogeneity <- TRUE
  }

  # Check for very high heterogeneity
  if (!is.null(i2) && i2 > 95) {
    issues$extreme_heterogeneity <- TRUE
  }

  # Check for questionable I2 estimates
  if (!is.null(i2)) {
    if (i2 < 0) {
      issues$negative_i2 <- TRUE
    }
    if (i2 > 100) {
      issues$i2_exceeds_100 <- TRUE
    }
  }

  # Check if df provided for Q-test
  if (!is.null(df) && !is.null(fit$QE)) {
    # Q-test p-value
    q_pval <- fit$QEp
    if (!is.null(q_pval) && q_pval > 0.05) {
      issues$not_significant_q <- TRUE
    }
  }

  has_issues <- length(issues) > 0

  list(
    has_issues = has_issues,
    issues = issues,
    tau2 = tau2,
    i2 = i2,
    issue_type = if (has_issues) "heterogeneity" else "none",
    severity = ifelse(
      !is.null(issues$extreme_heterogeneity),
      "high", "medium"
    ),
    message = if (has_issues) {
      paste("Heterogeneity issues:", paste(names(issues), collapse = ", "))
    } else {
      "Heterogeneity within normal range"
    }
  )
}

#' Check for Variance Collapse
#'
#' Detect when variance components collapse to zero or become problematic.
#'
#' @param fit Fitted model or test result
#' @return List with variance diagnostics
#' @export
check_variance_collapse <- function(fit) {

  # If fit is a test result list, extract the actual fit
  if (is.list(fit) && !is.null(fit$fit)) {
    fit <- fit$fit
  }

  issues <- list()

  # Check if fit exists
  if (is.null(fit) || !inherits(fit, "rma.mv")) {
    return(list(
      has_issues = TRUE,
      issue_type = "no_fit",
      severity = "critical",
      message = "Cannot check variance - no fit object"
    ))
  }

  # Check variance components
  if (!is.null(fit$tau2)) {
    if (fit$tau2 <= 0) {
      issues$tau2_nonpositive <- TRUE
    }
    if (fit$tau2 < 1e-10) {
      issues$tau2_near_zero <- TRUE
    }
  }

  if (!is.null(fit$sigma2)) {
    if (any(fit$sigma2 <= 0)) {
      issues$sigma2_nonpositive <- TRUE
    }
    if (any(fit$sigma2 < 1e-10)) {
      issues$sigma2_near_zero <- TRUE
    }
  }

  # Check if variance ratio is extreme
  if (!is.null(fit$tau2) && !is.null(fit$sigma2) && length(fit$sigma2) > 0) {
    ratio <- fit$tau2 / fit$sigma2[1]
    if (ratio > 1000) {
      issues$extreme_variance_ratio <- TRUE
    }
    if (ratio < 0.001) {
      issues$extreme_variance_ratio_low <- TRUE
    }
  }

  has_issues <- length(issues) > 0

  list(
    has_issues = has_issues,
    issues = issues,
    tau2 = fit$tau2,
    sigma2 = fit$sigma2,
    issue_type = if (has_issues) "variance_collapse" else "none",
    severity = ifelse(
      !is.null(issues$tau2_nonpositive),
      "high", "medium"
    ),
    message = if (has_issues) {
      paste("Variance issues:", paste(names(issues), collapse = ", "))
    } else {
      "Variance components adequate"
    }
  )
}

#' Detect All Issues in a Test Result
#'
#' Run all issue detection checks on a single test result.
#'
#' @param test_result Result from test_reml_basic, test_rve_cr2, etc.
#' @param df Original data frame (optional, for additional checks)
#' @return Comprehensive issue report
#' @export
detect_all_issues <- function(test_result, df = NULL) {

  issues <- list()

  # Convergence check
  conv <- check_convergence_problems(test_result)
  if (conv$has_issues) {
    issues$convergence <- conv
  }

  # Variance collapse check
  var_collapse <- check_variance_collapse(test_result)
  if (var_collapse$has_issues) {
    issues$variance <- var_collapse
  }

  # Heterogeneity check
  het <- check_heterogeneity_problems(test_result, df)
  if (het$has_issues) {
    issues$heterogeneity <- het
  }

  # Small sample check (requires df)
  if (!is.null(df)) {
    small_sample <- check_small_sample_issues(df, test_result)
    if (small_sample$has_issues) {
      issues$small_sample <- small_sample
    }
  }

  # Warning check
  if (length(test_result$warnings) > 0) {
    issues$warnings <- test_result$warnings
  }

  # Success check
  if (!test_result$success) {
    issues$fit_failed <- list(
      message = test_result$error %||% "Unknown error",
      severity = "critical"
    )
  }

  has_any_issues <- length(issues) > 0
  max_severity <- if (has_any_issues) {
    severities <- sapply(issues, function(x) x$severity %||% "medium")
    if (any(severities == "critical")) "critical" else if (any(severities == "high")) "high" else "medium"
  } else "none"

  list(
    has_issues = has_any_issues,
    severity = max_severity,
    issue_count = length(issues),
    issues = issues,
    method = test_result$method,
    summary = if (has_any_issues) {
      paste(sprintf("%s (%s)", names(issues), sapply(issues, function(x) x$message)), collapse = "; ")
    } else "No issues detected"
  )
}

#' Categorize Issue by Type
#'
#' @param issue_report Report from detect_all_issues()
#' @return Character vector of issue categories
#' @export
categorize_issue <- function(issue_report) {

  if (!issue_report$has_issues) {
    return("none")
  }

  categories <- character()

  if (!is.null(issue_report$issues$convergence)) {
    categories <- c(categories, "convergence")
  }
  if (!is.null(issue_report$issues$variance)) {
    categories <- c(categories, "variance")
  }
  if (!is.null(issue_report$issues$heterogeneity)) {
    categories <- c(categories, "heterogeneity")
  }
  if (!is.null(issue_report$issues$small_sample)) {
    categories <- c(categories, "small_sample")
  }
  if (!is.null(issue_report$issues$fit_failed)) {
    categories <- c(categories, "fit_failed")
  }

  paste(categories, collapse = ",")
}
