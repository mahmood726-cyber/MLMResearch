#' Run Cohort Tests on Dichotomous Outcomes
#'
#' Functions to run systematic tests on the dichotomous outcomes (logOR) cohort
#' from the MLM501 Cochrane reviews dataset.
#'
#' @name run_cohort_tests
NULL

#' Get Dichotomous Cohort from MLM501 Data
#'
#' Loads and filters the MLM501 effects table for dichotomous outcomes (logOR).
#'
#' @param mlm_data_path Path to mlm_effects.csv (NULL to use MLM501 package)
#' @param min_n_effects Minimum number of effects per analysis (default: 3)
#' @param min_n_studies Minimum number of studies per analysis (default: 2)
#' @return Filtered data frame
#' @export
get_dichotomous_cohort <- function(mlm_data_path = NULL,
                                    min_n_effects = 3,
                                    min_n_studies = 2) {

  # Try to load from MLM501 package first
  if (is.null(mlm_data_path)) {
    if (requireNamespace("MLM501", quietly = TRUE)) {
      df <- MLM501::read_mlm_effects()
    } else {
      # Try to find in current directory
      possible_paths <- c(
        "inst/extdata/mlm_effects.csv",
        "../MLM501/inst/extdata/mlm_effects.csv",
        "../../MLM501/inst/extdata/mlm_effects.csv"
      )
      for (p in possible_paths) {
        if (file.exists(p)) {
          df <- utils::read.csv(p, stringsAsFactors = FALSE)
          mlm_data_path <- p
          break
        }
      }
      if (is.null(mlm_data_path)) {
        stop("Could not find mlm_effects.csv. Provide path or install MLM501 package.")
      }
    }
  } else {
    df <- utils::read.csv(mlm_data_path, stringsAsFactors = FALSE)
  }

  # Filter for dichotomous outcomes with logOR measure
  coh <- df[df$outcome_type == "DICH" & df$measure == "logOR", ]

  # Drop rows with NA or infinite values
  coh <- coh[is.finite(coh$TE) & is.finite(coh$seTE) & coh$seTE > 0, ]

  # Filter by size criteria
  if (!is.null(min_n_effects) || !is.null(min_n_studies)) {
    analysis_stats <- aggregate(
      cbind(n_effects = rep(1, nrow(coh)),
            n_studies = as.numeric(!duplicated(coh$study_id))),
      by = list(analysis_id = coh$analysis_id),
      FUN = sum
    )

    valid_analyses <- analysis_stats$analysis_id
    if (!is.null(min_n_effects)) {
      valid_analyses <- analysis_stats$analysis_id[analysis_stats$n_effects >= min_n_effects]
    }
    if (!is.null(min_n_studies)) {
      valid_analyses <- analysis_stats$analysis_id[analysis_stats$n_studies >= min_n_studies]
    }

    coh <- coh[coh$analysis_id %in% valid_analyses, ]
  }

  cat(sprintf("Loaded dichotomous cohort: %d effects from %d analyses (%d studies)\n",
              nrow(coh), length(unique(coh$analysis_id)), length(unique(coh$study_id))))

  coh
}

#' Run Tests on a Single Analysis
#'
#' @param df Data frame for a single analysis_id
#' @param methods Methods to test
#' @return Result list
#' @export
run_single_analysis_tests <- function(df, methods = c("REML", "RVE-CR2", "RVE-CR1")) {

  analysis_id <- unique(df$analysis_id)

  result <- list(
    analysis_id = analysis_id,
    review_id = unique(df$review_id),
    n_effects = nrow(df),
    n_studies = length(unique(df$study_id)),
    n_clusters = length(unique(df$review_id)),
    test_results = test_all_methods(df, methods = methods, silent = TRUE)
  )

  result$summary <- result$test_results$summary
  result$summary$analysis_id <- analysis_id
  result$summary$review_id <- unique(df$review_id)

  result
}

#' Run Dichotomous Tests on All Analyses
#'
#' @param mlm_data_path Path to mlm_effects.csv
#' @param max_analyses Maximum number of analyses to test (NULL for all)
#' @param subsample_fraction Fraction of analyses to randomly sample (NULL for all)
#' @param methods Methods to test
#' @param progress Show progress bar
#' @return List of all test results
#' @export
run_dichotomous_tests <- function(mlm_data_path = NULL,
                                   max_analyses = NULL,
                                   subsample_fraction = NULL,
                                   methods = c("REML", "RVE-CR2", "RVE-CR1"),
                                   progress = TRUE) {

  # Load cohort
  coh <- get_dichotomous_cohort(mlm_data_path)

  analysis_ids <- unique(coh$analysis_id)

  # Subsample if requested
  if (!is.null(subsample_fraction) && subsample_fraction < 1) {
    n_sample <- ceiling(length(analysis_ids) * subsample_fraction)
    analysis_ids <- sample(analysis_ids, n_sample)
    cat(sprintf("Subsampled to %d analyses (%.0f%%)\n",
                n_sample, subsample_fraction * 100))
  }

  # Limit max analyses
  if (!is.null(max_analyses) && length(analysis_ids) > max_analyses) {
    analysis_ids <- analysis_ids[1:max_analyses]
    cat(sprintf("Limited to %d analyses\n", max_analyses))
  }

  cat(sprintf("Testing %d analyses with %d methods...\n",
              length(analysis_ids), length(methods)))

  # Run tests
  all_results <- list()
  all_summaries <- list()

  for (i in seq_along(analysis_ids)) {
    aid <- analysis_ids[i]

    if (progress) {
      cat(sprintf("\r[%d/%d] Testing analysis %s...", i, length(analysis_ids), aid))
    }

    df_i <- coh[coh$analysis_id == aid, ]

    result <- run_single_analysis_tests(df_i, methods = methods)

    all_results[[aid]] <- result
    all_summaries[[i]] <- result$summary
  }

  if (progress) cat("\n")

  # Combine all summaries
  combined_summary <- do.call(rbind, all_summaries)

  list(
    cohort_data = coh,
    n_analyses_tested = length(analysis_ids),
    analysis_ids = analysis_ids,
    results = all_results,
    summary = combined_summary,
    methods_tested = methods
  )
}

#' Generate Test Report
#'
#' @param test_results Result from run_dichotomous_tests()
#' @param output_path Path to save report (NULL to print)
#' @return Summary statistics
#' @export
generate_test_report <- function(test_results, output_path = NULL) {

  summary_df <- test_results$summary

  # Overall statistics
  total_tests <- nrow(summary_df)
  successful_tests <- sum(summary_df$success)
  converged_tests <- sum(summary_df$convergence == "converged")

  # By method
  by_method <- lapply(unique(summary_df$method), function(m) {
    df_m <- summary_df[summary_df$method == m, ]
    list(
      method = m,
      n_tests = nrow(df_m),
      n_success = sum(df_m$success),
      success_rate = mean(df_m$success),
      n_converged = sum(df_m$convergence == "converged"),
      convergence_rate = mean(df_m$convergence == "converged"),
      avg_time = mean(df_m$time_seconds, na.rm = TRUE),
      n_warnings = sum(df_m$n_warnings)
    )
  })

  by_method_df <- do.call(rbind, lapply(by_method, function(x) {
    data.frame(
      method = x$method,
      n_tests = x$n_tests,
      n_success = x$n_success,
      success_rate = round(x$success_rate, 3),
      n_converged = x$n_converged,
      convergence_rate = round(x$convergence_rate, 3),
      avg_time = round(x$avg_time, 2),
      n_warnings = x$n_warnings,
      stringsAsFactors = FALSE
    )
  }))

  report <- list(
    overall = list(
      total_analyses = test_results$n_analyses_tested,
      total_tests = total_tests,
      successful_tests = successful_tests,
      overall_success_rate = round(successful_tests / total_tests, 3),
      converged_tests = converged_tests,
      overall_convergence_rate = round(converged_tests / total_tests, 3)
    ),
    by_method = by_method_df,
    problematic_analyses = summary_df[!summary_df$success, ]
  )

  # Print report
  cat("\n=== Multilevel Meta-Analysis Test Report ===\n\n")
  cat("Overall Results:\n")
  cat(sprintf("  Analyses tested: %d\n", report$overall$total_analyses))
  cat(sprintf("  Total tests: %d\n", report$overall$total_tests))
  cat(sprintf("  Success rate: %.1f%%\n", report$overall$overall_success_rate * 100))
  cat(sprintf("  Convergence rate: %.1f%%\n\n", report$overall$overall_convergence_rate * 100))

  cat("Results by Method:\n")
  print(by_method_df, row.names = FALSE)

  if (nrow(report$problematic_analyses) > 0) {
    cat(sprintf("\nProblematic analyses: %d\n", nrow(report$problematic_analyses)))
  }

  # Save if path provided
  if (!is.null(output_path)) {
    utils::write.csv(report$by_method, file.path(output_path, "method_performance.csv"), row.names = FALSE)
    utils::write.csv(report$problematic_analyses, file.path(output_path, "problematic_analyses.csv"), row.names = FALSE)
    cat(sprintf("\nReport saved to: %s\n", output_path))
  }

  invisible(report)
}

#' Create Cohort Subsample for Quick Testing
#'
#' @param full_df Full cohort data frame
#' @param n_analyses Number of analyses to sample
#' @param seed Random seed
#' @return Subsampled data frame
#' @export
create_downsampled_cohort <- function(full_df, n_analyses = 50, seed = 123) {
  set.seed(seed)

  analysis_ids <- unique(full_df$analysis_id)
  if (length(analysis_ids) > n_analyses) {
    selected_ids <- sample(analysis_ids, n_analyses)
    cat(sprintf("Downsampled from %d to %d analyses\n",
                length(analysis_ids), n_analyses))
  } else {
    selected_ids <- analysis_ids
    cat(sprintf("Using all %d analyses\n", length(analysis_ids)))
  }

  full_df[full_df$analysis_id %in% selected_ids, ]
}
