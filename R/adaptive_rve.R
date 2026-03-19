#' Adaptive Robust Variance Estimation
#'
#' Improved RVE methods that adapt based on cluster count and data characteristics.
#' Addresses small-sample issues with standard CR2 method.
#'
#' @name adaptive_rve
NULL

#' Adaptive RVE Method
#'
#' Automatically selects the best RVE approach based on cluster count
#' and data characteristics.
#'
#' @param df Data frame with meta-analysis data
#' @param cluster_var Variable to use for clustering (default: "review_id")
#' @param min_clusters_cr2 Minimum clusters for CR2 (default: 15)
#' @param min_clusters_cr1 Minimum clusters for CR1 (default: 8)
#' @param fallback_method Method to use if RVE is not suitable (default: "REML")
#' @return List with fit and chosen method
#' @export
adaptive_rve <- function(df,
                         cluster_var = "review_id",
                         min_clusters_cr2 = 15,
                         min_clusters_cr1 = 8,
                         fallback_method = "REML") {

  n_clusters <- length(unique(df[[cluster_var]]))
  n_effects <- nrow(df)
  n_studies <- length(unique(df$study_id))

  # Determine appropriate method
  chosen_method <- if (n_clusters >= min_clusters_cr2) {
    "RVE-CR2"
  } else if (n_clusters >= min_clusters_cr1) {
    "RVE-CR1"
  } else {
    fallback_method
  }

  result <- switch(chosen_method,
    "RVE-CR2" = test_rve_cr2(df, cluster_var, silent = TRUE),
    "RVE-CR1" = test_rve_cr1(df, cluster_var, silent = TRUE),
    "REML" = test_reml_basic(df, cluster_var, silent = TRUE),
    list(success = FALSE, error = paste("Unknown method:", chosen_method))
  )

  # Add metadata about adaptive choice
  result$adaptive_choice <- list(
    n_clusters = n_clusters,
    n_effects = n_effects,
    n_studies = n_studies,
    chosen_method = chosen_method,
    reason = if (n_clusters < min_clusters_cr1) {
      "Fewer clusters than minimum for RVE"
    } else if (n_clusters < min_clusters_cr2) {
      "Insufficient clusters for CR2, using CR1"
    } else {
      "Sufficient clusters for optimal CR2"
    }
  )

  result$method <- paste0("Adaptive-", chosen_method)

  result
}

#' Small-Sample Adjusted RVE
#'
#' Applies additional small-sample corrections to RVE when cluster count
#' is limited.
#'
#' @param df Data frame with meta-analysis data
#' @param cluster_var Variable to use for clustering
#' @param type RVE type ("CR2", "CR1", or "auto")
#' @return List with fit and adjustments
#' @export
small_sample_rve <- function(df, cluster_var = "review_id", type = "auto") {

  n_clusters <- length(unique(df[[cluster_var]]))

  # Auto-select type based on clusters
  if (type == "auto") {
    type <- if (n_clusters >= 15) "CR2" else "CR1"
  }

  # Get base RVE result
  base_result <- if (type == "CR2") {
    test_rve_cr2(df, cluster_var, silent = TRUE)
  } else {
    test_rve_cr1(df, cluster_var, silent = TRUE)
  }

  if (!base_result$success) {
    return(base_result)
  }

  # Apply small-sample adjustments
  # For very small clusters, use a more conservative degrees of freedom
  if (n_clusters < 15) {
    # Use clubSandwich's constrain_df option for small samples
    tryCatch({
      V_adj <- clubSandwich::vcovCR(
        base_result$fit,
        cluster = df[[cluster_var]],
        type = type,
        constrain_df = FALSE  # Don't constrain degrees of freedom
      )

      coef_test_adj <- clubSandwich::coef_test(
        base_result$fit,
        vcov = V_adj,
        cluster = df[[cluster_var]]
      )

      base_result$rve$adjusted_vcov <- V_adj
      base_result$rve$adjusted_coef_test <- coef_test_adj
      base_result$small_sample_adjustment <- "Unconstrained degrees of freedom"
    }, error = function(e) {
      base_result$small_sample_adjustment <- paste("Adjustment failed:", e$message)
    })
  }

  base_result$method <- paste0("SS-RVE-", type)

  base_result
}

#' Bootstrap RVE for Small Samples
#'
#' Uses residual bootstrap to improve inference with very small samples.
#'
#' @param df Data frame with meta-analysis data
#' @param cluster_var Variable to use for clustering
#' @param n_bootstrap Number of bootstrap iterations (default: 500)
#' @param seed Random seed
#' @return List with bootstrap results
#' @export
bootstrap_rve <- function(df,
                          cluster_var = "review_id",
                          n_bootstrap = 500,
                          seed = 123) {

  set.seed(seed)

  # Fit base model
  base_fit <- test_reml_basic(df, cluster_var, silent = TRUE)

  if (!base_fit$success) {
    return(list(
      success = FALSE,
      error = base_fit$error,
      method = "Bootstrap-RVE"
    ))
  }

  # Get unique clusters
  clusters <- unique(df[[cluster_var]])
  n_clusters <- length(clusters)

  if (n_clusters < 5) {
    return(list(
      success = FALSE,
      error = "Too few clusters for bootstrap",
      method = "Bootstrap-RVE",
      n_clusters = n_clusters
    ))
  }

  # Bootstrap iterations
  bootstrap_estimates <- numeric(n_bootstrap)
  bootstrap_se <- numeric(n_bootstrap)

  for (i in seq_len(n_bootstrap)) {
    # Resample clusters with replacement
    boot_clusters <- sample(clusters, replace = TRUE)

    # Create bootstrap sample
    boot_df <- do.call(rbind, lapply(boot_clusters, function(c) {
      df[df[[cluster_var]] == c, , drop = FALSE]
    }))

    # Fit model to bootstrap sample
    boot_fit <- tryCatch({
      V <- diag(boot_df$seTE^2, nrow = nrow(boot_df))
      metafor::rma.mv(
        yi = boot_df$TE,
        V = V,
        random = list(~ 1 | boot_df[[cluster_var]], ~ 1 | boot_df$study_id),
        method = "REML",
        silent = TRUE
      )
    }, error = function(e) NULL)

    if (!is.null(boot_fit) && inherits(boot_fit, "rma.mv")) {
      bootstrap_estimates[i] <- boot_fit$beta
      bootstrap_se[i] <- boot_fit$se
    } else {
      bootstrap_estimates[i] <- NA
      bootstrap_se[i] <- NA
    }
  }

  # Calculate bootstrap statistics
  valid_estimates <- bootstrap_estimates[!is.na(bootstrap_estimates)]
  valid_se <- bootstrap_se[!is.na(bootstrap_se)]

  list(
    success = length(valid_estimates) > n_bootstrap / 2,
    method = "Bootstrap-RVE",
    n_bootstrap = n_bootstrap,
    n_valid = length(valid_estimates),
    bootstrap_mean = mean(valid_estimates),
    bootstrap_se = sd(valid_estimates),
    bootstrap_ci = stats::quantile(valid_estimates, c(0.025, 0.975)),
    base_estimate = base_fit$fit$beta,
    base_se = base_fit$fit$se,
    base_fit = base_fit
  )
}

#' Compare RVE Methods
#'
#' Run all RVE variants and compare performance.
#'
#' @param df Data frame with meta-analysis data
#' @param cluster_var Variable to use for clustering
#' @return Data frame with comparison results
#' @export
compare_rve_methods <- function(df, cluster_var = "review_id") {

  n_clusters <- length(unique(df[[cluster_var]]))

  # Run all RVE methods
  methods <- list(
    CR2 = test_rve_cr2(df, cluster_var, silent = TRUE),
    CR1 = test_rve_cr1(df, cluster_var, silent = TRUE),
    Adaptive = adaptive_rve(df, cluster_var, silent = TRUE),
    SmallSample = small_sample_rve(df, cluster_var, type = "auto")
  )

  # Optionally run bootstrap if clusters > 5
  if (n_clusters >= 5 && n_clusters <= 20) {
    methods$Bootstrap <- bootstrap_rve(df, cluster_var, n_bootstrap = 200)
  }

  # Compile comparison
  comparison <- data.frame(
    method = names(methods),
    success = sapply(methods, function(x) x$success),
    time_seconds = sapply(methods, function(x) x$time_seconds),
    n_clusters = n_clusters,
    stringsAsFactors = FALSE
  )

  # Add estimates for successful methods
  estimates <- sapply(methods, function(x) {
    if (x$success && !is.null(x$fit)) {
      if (!is.null(x$rve)) {
        x$rve$robust_se[1]
      } else {
        x$fit$se
      }
    } else {
      NA
    }
  })
  comparison$robust_se <- estimates

  list(
    comparison = comparison,
    details = methods,
    recommended = if (n_clusters >= 15) "CR2" else if (n_clusters >= 8) "CR1" else "Adaptive"
  )
}

#' Get Adaptive RVE Recommendation
#'
#' Returns a recommendation for which RVE method to use based on data.
#'
#' @param df Data frame with meta-analysis data
#' @param cluster_var Variable to use for clustering
#' @return List with recommendation and reasoning
#' @export
get_rve_recommendation <- function(df, cluster_var = "review_id") {

  n_clusters <- length(unique(df[[cluster_var]]))
  n_effects <- nrow(df)
  n_studies <- length(unique(df$study_id))

  # Determine recommendation
  recommendation <- if (n_clusters >= 20) {
    list(
      method = "CR2",
      confidence = "high",
      reasoning = sprintf("Large number of clusters (%d) - CR2 is optimal", n_clusters)
    )
  } else if (n_clusters >= 10) {
    list(
      method = "CR2",
      confidence = "medium",
      reasoning = sprintf("Adequate clusters (%d) - CR2 recommended with caution", n_clusters)
    )
  } else if (n_clusters >= 5) {
    list(
      method = "CR1 or Bootstrap",
      confidence = "medium",
      reasoning = sprintf("Few clusters (%d) - consider CR1 or bootstrap methods", n_clusters)
    )
  } else {
    list(
      method = "REML with caution",
      confidence = "low",
      reasoning = sprintf("Very few clusters (%d) - RVE not recommended", n_clusters)
    )
  }

  list(
    n_clusters = n_clusters,
    n_effects = n_effects,
    n_studies = n_studies,
    recommendation = recommendation
  )
}
