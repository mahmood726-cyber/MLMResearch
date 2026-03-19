#' Test Framework for Multilevel Meta-Analysis Methods
#'
#' A comprehensive testing framework to evaluate different meta-analysis methods
#' on Cochrane review datasets. Tracks convergence, warnings, and performance.
#'
#' @name test_framework
NULL

#' Test Basic REML Multilevel Model (metafor)
#'
#' @param df Data frame with columns: TE, seTE, review_id, analysis_id, study_id
#' @param cluster_var Variable to use for clustering (default: "review_id")
#' @param silent Suppress messages (default: FALSE)
#' @return List with fit object, convergence status, warnings, timing
#' @export
#' @examples
#' \dontrun{
#' result <- test_reml_basic(df)
#' print(result$summary)
#' }
test_reml_basic <- function(df, cluster_var = "review_id", silent = FALSE) {

  start_time <- Sys.time()
  warnings_list <- list()
  convergence <- "unknown"

  # Check required columns
  required_cols <- c("TE", "seTE", cluster_var)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    return(list(
      method = "REML",
      success = FALSE,
      error = paste("Missing columns:", paste(missing_cols, collapse = ", ")),
      fit = NULL,
      convergence = NA,
      warnings = list(),
      time_seconds = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    ))
  }

  # Try to fit the model
  fit <- tryCatch({

    # Build variance-covariance matrix
    V <- diag(df$seTE^2, nrow = nrow(df))

    # Fit multilevel model
    result <- withCallingHandlers({
      metafor::rma.mv(
        yi = df$TE,
        V = V,
        random = list(~ 1 | df[[cluster_var]], ~ 1 | df$study_id),
        method = "REML",
        silent = silent
      )
    },
    warning = function(w) {
      warnings_list[[length(warnings_list) + 1]] <<- w$message
      invokeRestart("muffleWarning")
    })

    # Check convergence
    convergence <- ifelse(!is.null(result$tau2) && is.finite(result$tau2), "converged", "failed")

    result
  },
  error = function(e) {
    convergence <<- "error"
    list(error = e$message)
  })

  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Check if fit contains error
  if (is.list(fit) && !is.null(fit$error)) {
    return(list(
      method = "REML",
      success = FALSE,
      error = fit$error,
      fit = NULL,
      convergence = "error",
      warnings = warnings_list,
      time_seconds = elapsed,
      n_studies = length(unique(df$study_id)),
      n_effects = nrow(df)
    ))
  }

  # Extract results
  summary_stats <- NULL
  if (!is.null(fit) && inherits(fit, "rma.mv")) {
    summary_stats <- list(
      estimate = fit$beta,
      se = fit$se,
      ci_lb = fit$ci.lb,
      ci_ub = fit$ci.ub,
      tau2 = fit$tau2,
      sigma2 = fit$sigma2,
      i2 = fit$I2,
      pval = fit$pval
    )
  }

  list(
    method = "REML",
    success = !is.null(fit) && inherits(fit, "rma.mv"),
    fit = fit,
    convergence = convergence,
    warnings = warnings_list,
    time_seconds = elapsed,
    summary = summary_stats,
    n_studies = length(unique(df$study_id)),
    n_effects = nrow(df),
    n_clusters = length(unique(df[[cluster_var]]))
  )
}

#' Test RVE-CR2 Method (clubSandwich)
#'
#' @param df Data frame with meta-analysis data
#' @param cluster_var Variable to use for clustering
#' @param silent Suppress messages
#' @return List with fit object, convergence status, warnings, timing
#' @export
test_rve_cr2 <- function(df, cluster_var = "review_id", silent = FALSE) {

  start_time <- Sys.time()
  warnings_list <- list()

  # First fit the base model
  base_fit <- test_reml_basic(df, cluster_var = cluster_var, silent = TRUE)

  if (!base_fit$success) {
    return(list(
      method = "RVE-CR2",
      success = FALSE,
      error = base_fit$error %||% "Base model failed",
      fit = NULL,
      convergence = "error",
      warnings = warnings_list,
      time_seconds = base_fit$time_seconds,
      n_clusters = base_fit$n_clusters
    ))
  }

  # Apply clubSandwich RVE
  rve_result <- tryCatch({

    V_club <- clubSandwich::vcovCR(
      base_fit$fit,
      cluster = df[[cluster_var]],
      type = "CR2"
    )

    # Get robust coefficient tests
    coef_test <- clubSandwich::coef_test(
      base_fit$fit,
      vcov = V_club,
      cluster = df[[cluster_var]]
    )

    list(
      vcov = V_club,
      coef_test = coef_test,
      robust_se = coef_test$SE,
      robust_pval = coef_test$pval,
      robust_ci_lb = coef_test$CI.L,
      robust_ci_ub = coef_test$CI.U
    )
  },
  error = function(e) {
    list(error = e$message)
  })

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  if (!is.null(rve_result$error)) {
    return(list(
      method = "RVE-CR2",
      success = FALSE,
      error = rve_result$error,
      fit = base_fit$fit,
      convergence = "error",
      warnings = warnings_list,
      time_seconds = elapsed,
      n_clusters = base_fit$n_clusters
    ))
  }

  list(
    method = "RVE-CR2",
    success = TRUE,
    fit = base_fit$fit,
    rve = rve_result,
    convergence = "converged",
    warnings = warnings_list,
    time_seconds = elapsed,
    n_clusters = base_fit$n_clusters,
    n_studies = base_fit$n_studies,
    n_effects = base_fit$n_effects
  )
}

#' Test RVE-CR1 Method (clubSandwich)
#'
#' @param df Data frame with meta-analysis data
#' @param cluster_var Variable to use for clustering
#' @param silent Suppress messages
#' @return List with fit object, convergence status, warnings, timing
#' @export
test_rve_cr1 <- function(df, cluster_var = "review_id", silent = FALSE) {

  start_time <- Sys.time()
  warnings_list <- list()

  # First fit the base model
  base_fit <- test_reml_basic(df, cluster_var = cluster_var, silent = TRUE)

  if (!base_fit$success) {
    return(list(
      method = "RVE-CR1",
      success = FALSE,
      error = base_fit$error %||% "Base model failed",
      fit = NULL,
      convergence = "error",
      warnings = warnings_list,
      time_seconds = base_fit$time_seconds
    ))
  }

  # Apply clubSandwich RVE with CR1
  rve_result <- tryCatch({

    V_club <- clubSandwich::vcovCR(
      base_fit$fit,
      cluster = df[[cluster_var]],
      type = "CR1"
    )

    coef_test <- clubSandwich::coef_test(
      base_fit$fit,
      vcov = V_club,
      cluster = df[[cluster_var]]
    )

    list(
      vcov = V_club,
      coef_test = coef_test,
      robust_se = coef_test$SE,
      robust_pval = coef_test$pval
    )
  },
  error = function(e) {
    list(error = e$message)
  })

  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  if (!is.null(rve_result$error)) {
    return(list(
      method = "RVE-CR1",
      success = FALSE,
      error = rve_result$error,
      fit = base_fit$fit,
      convergence = "error",
      warnings = warnings_list,
      time_seconds = elapsed
    ))
  }

  list(
    method = "RVE-CR1",
    success = TRUE,
    fit = base_fit$fit,
    rve = rve_result,
    convergence = "converged",
    warnings = warnings_list,
    time_seconds = elapsed,
    n_clusters = base_fit$n_clusters
  )
}

#' Test All Methods on a Single Dataset
#'
#' @param df Data frame with meta-analysis data
#' @param cluster_var Variable to use for clustering
#' @param methods Character vector of methods to test ("all" for all available)
#' @param silent Suppress messages
#' @return Data frame with results for all methods
#' @export
test_all_methods <- function(df, cluster_var = "review_id",
                             methods = c("REML", "RVE-CR2", "RVE-CR1"),
                             silent = FALSE) {

  if (!silent) {
    cat(sprintf("Testing %d methods on %d effects from %d studies...\n",
                length(methods), nrow(df), length(unique(df$study_id))))
  }

  results <- list()

  for (method in methods) {
    if (!silent) cat(sprintf("  Testing %s...", method))

    result <- switch(method,
      "REML" = test_reml_basic(df, cluster_var, silent = TRUE),
      "RVE-CR2" = test_rve_cr2(df, cluster_var, silent = TRUE),
      "RVE-CR1" = test_rve_cr1(df, cluster_var, silent = TRUE),
      list(success = FALSE, error = paste("Unknown method:", method))
    )

    results[[method]] <- result
    if (!silent) cat(sprintf(" %s\n", ifelse(result$success, "OK", "FAILED")))
  }

  # Convert to summary data frame
  summary_df <- do.call(rbind, lapply(names(results), function(m) {
    r <- results[[m]]
    data.frame(
      method = m,
      success = r$success,
      convergence = r$convergence,
      time_seconds = r$time_seconds,
      n_warnings = length(r$warnings),
      n_clusters = r$n_clusters %||% NA,
      n_studies = r$n_studies %||% NA,
      n_effects = r$n_effects %||% NA,
      stringsAsFactors = FALSE
    )
  }))

  list(
    results = results,
    summary = summary_df
  )
}

#' Summarize Results from Multiple Method Tests
#'
#' @param test_results Result from test_all_methods()
#' @return Summary statistics
#' @export
summarize_results <- function(test_results) {

  summary_df <- test_results$summary

  list(
    n_methods_tested = nrow(summary_df),
    n_successful = sum(summary_df$success),
    success_rate = mean(summary_df$success),
    n_converged = sum(summary_df$convergence == "converged"),
    convergence_rate = mean(summary_df$convergence == "converged"),
    total_warnings = sum(summary_df$n_warnings),
    avg_time = mean(summary_df$time_seconds, na.rm = TRUE),
    fastest_method = summary_df$method[which.min(summary_df$time_seconds)],
    slowest_method = summary_df$method[which.max(summary_df$time_seconds)],
    details = summary_df
  )
}

# Null coalescing operator helper
`%||%` <- function(x, y) if (is.null(x)) y else x
