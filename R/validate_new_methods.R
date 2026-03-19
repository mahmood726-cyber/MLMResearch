#' Method Validation Framework
#'
#' Functions to validate new and improved meta-analysis methods against
#' standard approaches using cross-validation and bootstrapping.
#'
#' @name validate_new_methods
NULL

#' Validate Method Comparison
#'
#' Compare a new method against standard methods on multiple datasets.
#'
#' @param datasets List of data frames, each with TE, seTE, and clustering
#' @param new_method Function to apply (signature: function(df) -> result)
#' @param comparison_methods Named list of comparison functions
#' @param metrics Metrics to compute ("bias", "mse", "coverage", "time")
#' @return Validation results
#' @export
validate_method_comparison <- function(datasets,
                                       new_method,
                                       comparison_methods = list(
                                         REML = function(df) test_reml_basic(df),
                                         CR2 = function(df) test_rve_cr2(df)
                                       ),
                                       metrics = c("bias", "mse", "coverage", "time")) {

  n_datasets <- length(datasets)
  all_methods <- c("New", names(comparison_methods))

  # Initialize results storage
  results <- list()

  for (i in seq_along(datasets)) {
    df <- datasets[[i]]
    results[[i]] <- list(
      dataset_id = names(datasets)[i],
      n_effects = nrow(df),
      n_studies = length(unique(df$study_id))
    )

    # Run new method
    new_result <- tryCatch({
      new_method(df)
    }, error = function(e) {
      list(success = FALSE, error = e$message)
    })

    # Run comparison methods
    for (m in names(comparison_methods)) {
      comp_result <- tryCatch({
        comparison_methods[[m]](df)
      }, error = function(e) {
        list(success = FALSE, error = e$message)
      })
      results[[i]][[m]] <- comp_result
    }

    results[[i]]$New <- new_result
  }

  # Summarize results
  summary <- summarize_validation(results, all_methods, metrics)

  list(
    results = results,
    summary = summary,
    n_datasets = n_datasets,
    metrics_computed = metrics
  )
}

#' Summarize Validation Results
#'
#' @param results Results from validate_method_comparison()
#' @param methods Method names
#' @param metrics Metrics to compute
#' @return Summary data frame
#' @export
summarize_validation <- function(results, methods, metrics) {

  n_methods <- length(methods)
  n_datasets <- length(results)

  # Success rate by method
  success_by_method <- sapply(methods, function(m) {
    successes <- sapply(results, function(r) {
      !is.null(r[[m]]) && r[[m]]$success
    })
    mean(successes, na.rm = TRUE)
  })

  # Average time by method
  time_by_method <- sapply(methods, function(m) {
    times <- sapply(results, function(r) {
      if (!is.null(r[[m]]) && !is.null(r[[m]]$time_seconds)) {
        r[[m]]$time_seconds
      } else {
        NA
      }
    })
    mean(times, na.rm = TRUE)
  })

  # Convergence rate by method
  convergence_by_method <- sapply(methods, function(m) {
    converged <- sapply(results, function(r) {
      !is.null(r[[m]]) && !is.null(r[[m]]$convergence) &&
        r[[m]]$convergence == "converged"
    })
    mean(converged, na.rm = TRUE)
  })

  data.frame(
    method = methods,
    success_rate = round(success_by_method, 3),
    convergence_rate = round(convergence_by_method, 3),
    avg_time = round(time_by_method, 2),
    stringsAsFactors = FALSE
  )
}

#' Bootstrap Validation
#'
#' Validate methods using bootstrap resampling to estimate variability.
#'
#' @param df Original data frame
#' @param methods Named list of methods to validate
#' @param n_bootstrap Number of bootstrap samples (default: 500)
#' @param seed Random seed
#' @return Bootstrap validation results
#' @export
bootstrap_validation <- function(df,
                                  methods = list(
                                    REML = function(d) test_reml_basic(d),
                                    CR2 = function(d) test_rve_cr2(d),
                                    Adaptive = function(d) adaptive_rve(d)
                                  ),
                                  n_bootstrap = 500,
                                  seed = 123) {

  set.seed(seed)

  n_methods <- length(methods)
  bootstrap_results <- list()

  # Get unique studies for resampling
  studies <- unique(df$study_id)
  n_studies <- length(studies)

  for (i in seq_len(n_bootstrap)) {
    # Bootstrap sample: resample studies with replacement
    boot_studies <- sample(studies, replace = TRUE)
    boot_df <- do.call(rbind, lapply(boot_studies, function(s) {
      df[df$study_id == s, , drop = FALSE]
    }))

    # Run each method on bootstrap sample
    for (m_name in names(methods)) {
      result <- tryCatch({
        methods[[m_name]](boot_df)
      }, error = function(e) {
        list(success = FALSE, error = e$message)
      })

      # Store estimate if successful
      if (result$success && !is.null(result$fit)) {
        if (!is.null(result$rve)) {
          estimate <- result$fit$beta
          se <- result$rve$robust_se[1]
        } else {
          estimate <- result$fit$beta
          se <- result$fit$se
        }
      } else {
        estimate <- NA
        se <- NA
      }

      if (is.null(bootstrap_results[[m_name]])) {
        bootstrap_results[[m_name]] <- list(
          estimates = numeric(n_bootstrap),
          ses = numeric(n_bootstrap)
        )
      }

      bootstrap_results[[m_name]]$estimates[i] <- estimate
      bootstrap_results[[m_name]]$ses[i] <- se
    }
  }

  # Compute bootstrap statistics
  summary <- lapply(bootstrap_results, function(res) {
    valid_est <- res$estimates[!is.na(res$estimates)]
    valid_se <- res$ses[!is.na(res$ses)]

    list(
      n_valid = length(valid_est),
      mean_estimate = mean(valid_est),
      sd_estimate = sd(valid_est),
      mean_se = mean(valid_se),
      ci_estimate = quantile(valid_est, c(0.025, 0.975)),
      se_coverage = mean(valid_est >= quantile(valid_est, 0.025) &
                          valid_est <= quantile(valid_est, 0.975))
    )
  })

  list(
    summary = summary,
    n_bootstrap = n_bootstrap,
    n_studies = n_studies,
    bootstrap_results = bootstrap_results
  )
}

#' Cross-Validation for Method Selection
#'
#' Use k-fold cross-validation to compare method predictive performance.
#'
#' @param df Data frame with all data
#' @param methods Named list of methods to compare
#' @param k_folds Number of folds (default: 5)
#' @param seed Random seed
#' @return Cross-validation results
#' @export
cv_method_comparison <- function(df,
                                 methods = list(
                                   REML = function(d) test_reml_basic(d),
                                   CR2 = function(d) test_rve_cr2(d)
                                 ),
                                 k_folds = 5,
                                 seed = 123) {

  set.seed(seed)

  # Get unique studies for folding
  studies <- unique(df$study_id)
  n_studies <- length(studies)

  # Create folds
  fold_ids <- sample(rep(1:k_folds, length.out = n_studies))
  names(fold_ids) <- studies

  # Storage for CV results
  cv_results <- list()

  for (fold in 1:k_folds) {
    # Split data
    test_studies <- names(fold_ids)[fold_ids == fold]
    train_studies <- names(fold_ids)[fold_ids != fold]

    train_df <- df[df$study_id %in% train_studies, ]
    test_df <- df[df$study_id %in% test_studies, ]

    # Fit each method on training data
    for (m_name in names(methods)) {
      fit <- tryCatch({
        methods[[m_name]](train_df)
      }, error = function(e) {
        list(success = FALSE, error = e$message)
      })

      if (is.null(cv_results[[m_name]])) {
        cv_results[[m_name]] <- list(
          prediction_errors = list(),
          fold_estimates = numeric(k_folds)
        )
      }

      if (fit$success && !is.null(fit$fit)) {
        # Get estimate
        estimate <- fit$fit$beta

        # Compute prediction error on test set
        test_effects <- test_df$TE
        prediction_errors <- (test_effects - estimate)^2
        mse <- mean(prediction_errors)

        cv_results[[m_name]]$prediction_errors[[fold]] <- prediction_errors
        cv_results[[m_name]]$fold_estimates[fold] <- estimate
      } else {
        cv_results[[m_name]]$prediction_errors[[fold]] <- NA
        cv_results[[m_name]]$fold_estimates[fold] <- NA
      }
    }
  }

  # Compute CV statistics
  cv_summary <- lapply(cv_results, function(res) {
    all_errors <- unlist(res$prediction_errors)
    valid_errors <- all_errors[!is.na(all_errors)]

    list(
      mean_mse = mean(valid_errors),
      sd_mse = sd(valid_errors),
      n_valid_folds = sum(!is.na(res$fold_estimates)),
      fold_estimates = res$fold_estimates
    )
  })

  # Rank methods by MSE
  mses <- sapply(cv_summary, function(x) x$mean_mse)
  ranking <- names(sort(mses))

  list(
    cv_summary = cv_summary,
    ranking = ranking,
    best_method = ranking[1],
    k_folds = k_folds
  )
}

#' Simulation Study for Method Validation
#'
#' Run simulation studies to validate method properties under known conditions.
#'
#' @param scenarios List of simulation scenarios
#' @param methods Methods to test
#' @param n_sim Number of simulations per scenario (default: 1000)
#' @param seed Random seed
#' @return Simulation results
#' @export
simulation_study <- function(scenarios,
                             methods = list(
                               REML = function(d) test_reml_basic(d),
                               CR2 = function(d) test_rve_cr2(d)
                             ),
                             n_sim = 1000,
                             seed = 123) {

  set.seed(seed)

  results <- list()

  for (s_name in names(scenarios)) {
    scenario <- scenarios[[s_name]]
    cat(sprintf("Running scenario: %s (%d simulations)\n", s_name, n_sim))

    scenario_results <- list()

    for (m_name in names(methods)) {
      estimates <- numeric(n_sim)
      ses <- numeric(n_sim)
      coverage <- logical(n_sim)

      for (i in seq_len(n_sim)) {
        # Generate data under scenario
        df <- generate_simulated_data(scenario)

        # Fit method
        fit <- tryCatch({
          methods[[m_name]](df)
        }, error = function(e) {
          list(success = FALSE)
        })

        if (fit$success && !is.null(fit$fit)) {
          estimates[i] <- fit$fit$beta
          ses[i] <- fit$fit$se

          # Check coverage
          if (!is.null(fit$fit$ci.lb) && !is.null(fit$fit$ci.ub)) {
            coverage[i] <- (scenario$true_effect >= fit$fit$ci.lb) &&
                           (scenario$true_effect <= fit$fit$ci.ub)
          } else {
            coverage[i] <- NA
          }
        } else {
          estimates[i] <- NA
          ses[i] <- NA
          coverage[i] <- NA
        }
      }

      # Compute statistics
      valid <- !is.na(estimates)

      scenario_results[[m_name]] <- list(
        bias = mean(estimates[valid] - scenario$true_effect, na.rm = TRUE),
        empirical_se = sd(estimates[valid], na.rm = TRUE),
        avg_se = mean(ses[valid], na.rm = TRUE),
        coverage_rate = mean(coverage[valid], na.rm = TRUE),
        convergence_rate = mean(valid),
        n_sim = n_sim
      )
    }

    results[[s_name]] <- scenario_results
  }

  results
}

#' Generate Simulated Data for Validation
#'
#' @param scenario List with simulation parameters
#' @return Simulated data frame
#' @export
generate_simulated_data <- function(scenario) {

  n_studies <- scenario$n_studies %||% 20
  n_clusters <- scenario$n_clusters %||% 5
  true_effect <- scenario$true_effect %||% 0
  tau2 <- scenario$tau2 %||% 0.1
  sigma2 <- scenario$sigma2 %||% 0.05

  # Assign studies to clusters
  cluster_assignments <- sample(1:n_clusters, n_studies, replace = TRUE)

  # Generate effects
  study_effects <- true_effect +
    rnorm(n_clusters, 0, sqrt(tau2))[cluster_assignments] +
    rnorm(n_studies, 0, sqrt(sigma2))

  # Generate sampling variances
  se <- scenario$se %||% runif(n_studies, 0.1, 0.3)
  vi <- se^2

  # Observed effects
  yi <- study_effects + rnorm(n_studies, 0, se)

  data.frame(
    TE = yi,
    seTE = se,
    study_id = paste0("S", 1:n_studies),
    review_id = paste0("C", cluster_assignments)
  )
}

#' Generate Validation Report
#'
#' @param validation_results Results from validate_method_comparison()
#' @param output_path Optional path to save report
#' @return Report summary
#' @export
generate_validation_report <- function(validation_results, output_path = NULL) {

  summary <- validation_results$summary

  cat("\n=== Method Validation Report ===\n\n")
  cat(sprintf("Datasets tested: %d\n\n", validation_results$n_datasets))

  cat("Performance Summary:\n")
  print(summary, row.names = FALSE)

  # Find best method by success rate
  best_by_success <- summary$method[which.max(summary$success_rate)]
  best_by_time <- summary$method[which.min(summary$avg_time)]

  cat(sprintf("\nBest method (success rate): %s (%.1f%%)\n",
              best_by_success,
              max(summary$success_rate) * 100))
  cat(sprintf("Fastest method: %s (%.2f sec)\n",
              best_by_time,
              min(summary$avg_time)))

  # Save if requested
  if (!is.null(output_path)) {
    utils::write.csv(summary, file.path(output_path, "validation_summary.csv"), row.names = FALSE)
    cat(sprintf("\nReport saved to: %s\n", output_path))
  }

  invisible(summary)
}

#' Null coalescing operator (internal)
`%||%` <- function(x, y) if (is.null(x)) y else x
