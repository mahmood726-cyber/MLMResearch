#' Pattern Analysis Functions for Meta-Analysis Testing
#'
#' Functions to analyze patterns in test results and identify conditions
#' that lead to method failures or issues.
#'
#' @name analyze_patterns
NULL

#' Analyze Method Performance Patterns
#'
#' @param test_results Result from run_dichotomous_tests()
#' @return List with performance analysis
#' @export
analyze_method_patterns <- function(test_results) {

  summary_df <- test_results$summary

  # Overall success rates by method
  method_stats <- summary_df %>%
    stats::aggregate(list(
      n_tests = method,
      n_success = success,
      n_converged = (convergence == "converged"),
      avg_time = time_seconds
    ), by = list(method = method), FUN = function(x) {
      if (is.logical(x)) sum(x) else mean(x, na.rm = TRUE)
    })

  method_stats$success_rate <- method_stats$n_success / method_stats$n_tests
  method_stats$convergence_rate <- method_stats$n_converged / method_stats$n_tests

  # Find worst and best performing methods
  method_stats <- method_stats[order(-method_stats$success_rate), ]
  best_method <- method_stats$method[1]
  worst_method <- method_stats$method[nrow(method_stats)]

  # Method agreement analysis
  successful_methods <- summary_df[summary_df$success, ]
  if (nrow(successful_methods) > 0) {
    agreement <- stats::aggregate(
      n_success ~ analysis_id,
      data = successful_methods,
      FUN = length
    )
    names(agreement)[2] <- "n_successful_methods"
    agreement_rate <- mean(agreement$n_successful_methods) / length(unique(summary_df$method))
  } else {
    agreement <- data.frame(analysis_id = character(), n_successful_methods = integer())
    agreement_rate <- NA
  }

  list(
    method_stats = method_stats,
    best_method = best_method,
    worst_method = worst_method,
    method_agreement_rate = agreement_rate,
    agreement_details = agreement
  )
}

#' Analyze Failure Patterns
#'
#' Identify conditions that lead to method failures.
#'
#' @param test_results Result from run_dichomotous_tests()
#' @param cohort_data Original cohort data frame
#' @return List with failure pattern analysis
#' @export
analyze_failure_patterns <- function(test_results, cohort_data) {

  summary_df <- test_results$summary
  failures <- summary_df[!summary_df$success, ]

  if (nrow(failures) == 0) {
    return(list(
      has_failures = FALSE,
      message = "No failures detected"
    ))
  }

  # Get analysis details for failures
  analysis_stats <- aggregate(
    cbind(n_effects = rep(1, nrow(cohort_data)),
          n_studies = as.numeric(!duplicated(cohort_data$study_id))),
    by = list(analysis_id = cohort_data$analysis_id),
    FUN = sum
  )

  # Merge with failures
  failure_details <- merge(
    failures,
    analysis_stats,
    by = "analysis_id",
    all.x = TRUE
  )

  # Analyze failure by data characteristics
  fail_by_n_effects <- stats::cut(failure_details$n_effects,
                                    breaks = c(0, 5, 10, 20, Inf),
                                    labels = c("<5", "5-10", "10-20", ">20"))
  fail_by_n_studies <- stats::cut(failure_details$n_studies,
                                   breaks = c(0, 3, 5, 10, Inf),
                                   labels = c("<3", "3-5", "5-10", ">10"))
  fail_by_n_clusters <- stats::cut(failure_details$n_clusters,
                                    breaks = c(0, 5, 10, 20, Inf),
                                    labels = c("<5", "5-10", "10-20", ">20"))

  # Failure rates by characteristic
  failure_by_n_effects <- table(fail_by_n_effects)
  failure_by_n_studies <- table(fail_by_n_studies)
  failure_by_n_clusters <- table(fail_by_n_clusters)

  # Method-specific failures
  failures_by_method <- table(failures$method)

  list(
    has_failures = TRUE,
    n_failures = nrow(failures),
    failure_rate = nrow(failures) / nrow(summary_df),
    failure_details = failure_details,
    failures_by_method = failures_by_method,
    failure_by_n_effects = failure_by_n_effects,
    failure_by_n_studies = failure_by_n_studies,
    failure_by_n_clusters = failure_by_n_clusters,
    summary = list(
      most_failed_method = names(which.max(failures_by_method)),
      common_failure_size = names(which.max(failure_by_n_effects)),
      common_failure_clusters = names(which.max(failure_by_n_clusters))
    )
  )
}

#' Analyze Computational Performance
#'
#' @param test_results Result from run_dichotomous_tests()
#' @return List with performance analysis
#' @export
analyze_performance <- function(test_results) {

  summary_df <- test_results$summary

  # Performance by method
  perf_by_method <- stats::aggregate(
    time_seconds ~ method,
    data = summary_df,
    FUN = function(x) c(
      mean = mean(x, na.rm = TRUE),
      median = median(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE)
    )
  )

  # Overall performance
  total_time <- sum(summary_df$time_seconds, na.rm = TRUE)
  avg_time_per_analysis <- mean(summary_df$time_seconds, na.rm = TRUE)

  # Warning analysis
  total_warnings <- sum(summary_df$n_warnings, na.rm = TRUE)

  list(
    total_time_seconds = total_time,
    avg_time_per_analysis = avg_time_per_analysis,
    total_warnings = total_warnings,
    performance_by_method = perf_by_method,
    fastest_method = summary_df$method[which.min(summary_df$time_seconds)],
    slowest_method = summary_df$method[which.max(summary_df$time_seconds)]
  )
}

#' Generate Comparison Plot Data
#'
#' Prepare data for plotting method comparisons.
#'
#' @param test_results Result from run_dichotomous_tests()
#' @return Data frame ready for plotting
#' @export
get_plot_data <- function(test_results) {

  summary_df <- test_results$summary

  # Success rate by method
  success_data <- aggregate(
    success ~ method,
    data = summary_df,
    FUN = mean
  )
  names(success_data)[2] <- "success_rate"

  # Convergence rate by method
  conv_data <- aggregate(
    (convergence == "converged") ~ method,
    data = summary_df,
    FUN = mean
  )
  names(conv_data)[2] <- "convergence_rate"

  # Avg time by method
  time_data <- aggregate(
    time_seconds ~ method,
    data = summary_df,
    FUN = mean
  )

  # Merge
  plot_data <- merge(success_data, conv_data, by = "method")
  plot_data <- merge(plot_data, time_data, by = "method")

  plot_data
}

#' Identify Problematic Analyses
#'
#' Find analyses that consistently fail across methods.
#'
#' @param test_results Result from run_dichotomous_tests()
#' @param threshold Minimum number of methods that must fail for an analysis to be flagged
#' @return Data frame of problematic analyses
#' @export
identify_problematic_analyses <- function(test_results, threshold = 2) {

  summary_df <- test_results$summary

  # Count failures per analysis
  failures_per_analysis <- aggregate(
    !success ~ analysis_id + review_id,
    data = summary_df,
    FUN = sum
  )
  names(failures_per_analysis)[3] <- "n_failures"

  # Count total methods per analysis
  total_methods <- aggregate(
    method ~ analysis_id,
    data = summary_df,
    FUN = length
  )
  names(total_methods)[2] <- "n_methods"

  # Merge
  problematic <- merge(failures_per_analysis, total_methods, by = "analysis_id")

  # Flag those exceeding threshold
  problematic$is_problematic <- problematic$n_failures >= threshold

  # Get details
  problematic_details <- problematic[problematic$is_problematic, ]
  problematic_details <- problematic_details[order(-problematic_details$n_failures), ]

  list(
    n_problematic = sum(problematic$is_problematic),
    problematic_rate = mean(problematic$is_problematic),
    details = problematic_details,
    flagged_analyses = problematic_details$analysis_id
  )
}

#' Correlate Issues with Data Characteristics
#'
#' @param test_results Result from run_dichotomous_tests()
#' @param cohort_data Original cohort data frame
#' @return List with correlations
#' @export
analyze_issue_correlations <- function(test_results, cohort_data) {

  summary_df <- test_results$summary

  # Get analysis characteristics
  analysis_stats <- aggregate(
    cbind(
      n_effects = rep(1, nrow(cohort_data)),
      n_studies = as.numeric(!duplicated(cohort_data$study_id))
    ),
    by = list(analysis_id = cohort_data$analysis_id),
    FUN = sum
  )

  # Get cluster info
  cluster_info <- unique(cohort_data[, c("analysis_id", "review_id")])
  analysis_stats <- merge(analysis_stats, cluster_info, by = "analysis_id")

  # Count clusters per analysis
  analysis_stats$n_clusters <- sapply(analysis_stats$review_id, function(x) {
    sum(cohort_data$review_id == x)
  })
  analysis_stats$n_unique_clusters <- length(unique(analysis_stats$review_id))

  # Merge with test results
  merged <- merge(summary_df, analysis_stats, by = "analysis_id")

  # Calculate correlations
  success_numeric <- as.numeric(merged$success)
  converged_numeric <- as.numeric(merged$convergence == "converged")

  cor_success_effects <- cor(success_numeric, merged$n_effects, use = "complete.obs")
  cor_success_studies <- cor(success_numeric, merged$n_studies, use = "complete.obs")
  cor_success_clusters <- cor(success_numeric, merged$n_unique_clusters, use = "complete.obs")

  cor_time_effects <- cor(merged$time_seconds, merged$n_effects, use = "complete.obs")

  list(
    success_vs_n_effects = cor_success_effects,
    success_vs_n_studies = cor_success_studies,
    success_vs_n_clusters = cor_success_clusters,
    time_vs_n_effects = cor_time_effects,
    interpretation = list(
      success_effects = ifelse(cor_success_effects < -0.1, "More effects -> lower success", "No clear pattern"),
      success_studies = ifelse(cor_success_studies < -0.1, "More studies -> lower success", "No clear pattern"),
      success_clusters = ifelse(cor_success_clusters < -0.1, "More clusters -> lower success", "No clear pattern")
    )
  )
}

#' Generate Comprehensive Analysis Report
#'
#' @param test_results Result from run_dichotomous_tests()
#' @param cohort_data Original cohort data frame
#' @param output_path Optional path to save report
#' @return List with all analyses
#' @export
analyze_all_patterns <- function(test_results, cohort_data = NULL, output_path = NULL) {

  if (is.null(cohort_data)) {
    cohort_data <- test_results$cohort_data
  }

  # Run all analyses
  method_patterns <- analyze_method_patterns(test_results)
  failure_patterns <- analyze_failure_patterns(test_results, cohort_data)
  performance <- analyze_performance(test_results)
  problematic <- identify_problematic_analyses(test_results)
  correlations <- analyze_issue_correlations(test_results, cohort_data)
  plot_data <- get_plot_data(test_results)

  # Compile report
  report <- list(
    method_patterns = method_patterns,
    failure_patterns = failure_patterns,
    performance = performance,
    problematic_analyses = problematic,
    correlations = correlations,
    plot_data = plot_data
  )

  # Print summary
  cat("\n=== Comprehensive Pattern Analysis Report ===\n\n")

  cat("Method Performance:\n")
  cat(sprintf("  Best method: %s (%.1f%% success)\n",
              method_patterns$best_method,
              method_patterns$method_stats$success_rate[1] * 100))
  cat(sprintf("  Worst method: %s (%.1f%% success)\n",
              method_patterns$worst_method,
              method_patterns$method_stats$success_rate[nrow(method_patterns$method_stats)] * 100))
  cat(sprintf("  Method agreement: %.1f%%\n\n", method_patterns$method_agreement_rate * 100))

  if (failure_patterns$has_failures) {
    cat("Failure Patterns:\n")
    cat(sprintf("  Total failures: %d (%.1f%%)\n", failure_patterns$n_failures, failure_patterns$failure_rate * 100))
    cat(sprintf("  Most failed method: %s\n", failure_patterns$summary$most_failed_method))
    cat(sprintf("  Common failure size: %s effects\n\n", failure_patterns$summary$common_failure_size))
  }

  cat("Performance:\n")
  cat(sprintf("  Total time: %.1f seconds\n", performance$total_time_seconds))
  cat(sprintf("  Avg per analysis: %.2f seconds\n", performance$avg_time_per_analysis))
  cat(sprintf("  Total warnings: %d\n\n", performance$total_warnings))

  cat("Problematic Analyses:\n")
  cat(sprintf("  %d analyses flagged (%.1f%%)\n",
              problematic$n_problematic, problematic$problematic_rate * 100))

  # Save if requested
  if (!is.null(output_path)) {
    utils::write.csv(method_patterns$method_stats,
                     file.path(output_path, "method_performance.csv"), row.names = FALSE)
    utils::write.csv(failure_patterns$failure_details,
                     file.path(output_path, "failure_details.csv"), row.names = FALSE)
    utils::write.csv(performance$performance_by_method,
                     file.path(output_path, "performance_by_method.csv"), row.names = FALSE)
    utils::write.csv(plot_data,
                     file.path(output_path, "plot_data.csv"), row.names = FALSE)
    cat(sprintf("\nReport saved to: %s\n", output_path))
  }

  invisible(report)
}
