#' Issue Tracking System for Meta-Analysis Testing
#'
#' Functions to log, track, and report issues found during meta-analysis
#' method testing on Cochrane review datasets.
#'
#' @name issue_tracker
NULL

# Global issue database (in-memory storage)
.issue_database <- new.env(parent = emptyenv())
.issue_database$issues <- data.frame(
  timestamp = character(),
  review_id = character(),
  analysis_id = character(),
  method = character(),
  issue_type = character(),
  severity = character(),
  n_effects = integer(),
  n_studies = integer(),
  n_clusters = integer(),
  message = character(),
  details = character(),
  stringsAsFactors = FALSE
)

#' Log an Issue to the Database
#'
#' @param review_id Review identifier
#' @param analysis_id Analysis identifier
#' @param method Method name (e.g., "REML", "RVE-CR2")
#' @param issue_type Type of issue (e.g., "convergence", "variance", "small_sample")
#' @param severity Severity level ("critical", "high", "medium", "low")
#' @param message Brief description of the issue
#' @param details Additional details (can be a list or string)
#' @param n_effects Number of effects in the analysis
#' @param n_studies Number of studies
#' @param n_clusters Number of clusters
#' @return Invisible TRUE
#' @export
log_issue <- function(review_id, analysis_id, method, issue_type, severity,
                      message, details = "", n_effects = NA, n_studies = NA,
                      n_clusters = NA) {

  # Convert details to character if it's a list
  if (is.list(details)) {
    details <- paste(names(details), unlist(details), sep = "=", collapse = "; ")
  }

  new_issue <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    review_id = as.character(review_id),
    analysis_id = as.character(analysis_id),
    method = as.character(method),
    issue_type = as.character(issue_type),
    severity = as.character(severity),
    n_effects = as.integer(n_effects),
    n_studies = as.integer(n_studies),
    n_clusters = as.integer(n_clusters),
    message = as.character(message),
    details = as.character(details),
    stringsAsFactors = FALSE
  )

  .issue_database$issues <- rbind(.issue_database$issues, new_issue)

  invisible(TRUE)
}

#' Log Issues from Test Result
#'
#' Automatically extract and log issues from a test result.
#'
#' @param test_result Result from test_reml_basic, test_rve_cr2, etc.
#' @param review_id Review identifier
#' @param analysis_id Analysis identifier
#' @param df Original data frame (for additional info)
#' @return Number of issues logged
#' @export
log_test_result <- function(test_result, review_id, analysis_id, df = NULL) {

  # Run issue detection
  issue_report <- detect_all_issues(test_result, df)

  n_logged <- 0

  if (issue_report$has_issues) {
    # Log each issue category
    for (issue_name in names(issue_report$issues)) {
      issue_data <- issue_report$issues[[issue_name]]

      log_issue(
        review_id = review_id,
        analysis_id = analysis_id,
        method = test_result$method,
        issue_type = issue_name,
        severity = issue_data$severity,
        message = issue_data$message,
        details = paste(names(issue_data$issues), collapse = ","),
        n_effects = test_result$n_effects,
        n_studies = test_result$n_studies,
        n_clusters = test_result$n_clusters
      )
      n_logged <- n_logged + 1
    }
  }

  n_logged
}

#' Get All Logged Issues
#'
#' @return Data frame of all logged issues
#' @export
get_all_issues <- function() {
  .issue_database$issues
}

#' Get Issues by Method
#'
#' @param method Method name (e.g., "REML", "RVE-CR2")
#' @return Data frame of issues for the specified method
#' @export
get_issues_by_method <- function(method) {
  issues <- .issue_database$issues
  issues[issues$method == method, ]
}

#' Get Issues by Severity
#'
#' @param severity Severity level ("critical", "high", "medium", "low")
#' @return Data frame of issues with the specified severity
#' @export
get_issues_by_severity <- function(severity) {
  issues <- .issue_database$issues
  issues[issues$severity == severity, ]
}

#' Get Issues by Type
#'
#' @param issue_type Issue type (e.g., "convergence", "variance")
#' @return Data frame of issues of the specified type
#' @export
get_issues_by_type <- function(issue_type) {
  issues <- .issue_database$issues
  # Handle comma-separated types
  issues[sapply(issues$issue_type, function(x) grepl(issue_type, x)), ]
}

#' Get Critical Issues
#'
#' @return Data frame of critical and high severity issues
#' @export
get_critical_issues <- function() {
  issues <- .issue_database$issues
  issues[issues$severity %in% c("critical", "high"), ]
}

#' Clear the Issue Database
#'
#' @return Invisible NULL
#' @export
clear_issues <- function() {
  .issue_database$issues <- data.frame(
    timestamp = character(),
    review_id = character(),
    analysis_id = character(),
    method = character(),
    issue_type = character(),
    severity = character(),
    n_effects = integer(),
    n_studies = integer(),
    n_clusters = integer(),
    message = character(),
    details = character(),
    stringsAsFactors = FALSE
  )
  invisible(NULL)
}

#' Generate Issue Report
#'
#' Create a comprehensive report of all logged issues.
#'
#' @param output_path Optional path to save the report
#' @return List with report statistics and tables
#' @export
generate_issue_report <- function(output_path = NULL) {

  issues <- .issue_database$issues

  if (nrow(issues) == 0) {
    cat("No issues logged.\n")
    return(
      list(
        total_issues = 0,
        by_severity = data.frame(),
        by_method = data.frame(),
        by_type = data.frame(),
        details = data.frame()
      )
    )
  }

  # Summary by severity
  by_severity <- as.data.frame(table(issues$severity))
  names(by_severity) <- c("severity", "count")
  by_severity <- by_severity[order(-by_severity$count), ]

  # Summary by method
  by_method <- as.data.frame(table(issues$method))
  names(by_method) <- c("method", "count")
  by_method <- by_method[order(-by_method$count), ]

  # Summary by issue type
  by_type <- as.data.frame(table(issues$issue_type))
  names(by_type) <- c("issue_type", "count")
  by_type <- by_type[order(-by_type$count), ]

  # Statistics
  n_critical <- sum(issues$severity == "critical")
  n_high <- sum(issues$severity == "high")
  n_affected_analyses <- length(unique(issues$analysis_id))

  report <- list(
    total_issues = nrow(issues),
    n_critical = n_critical,
    n_high = n_high,
    n_affected_analyses = n_affected_analyses,
    by_severity = by_severity,
    by_method = by_method,
    by_type = by_type,
    details = issues
  )

  # Print report
  cat("\n=== Issue Report ===\n\n")
  cat(sprintf("Total issues logged: %d\n", report$total_issues))
  cat(sprintf("Critical issues: %d\n", report$n_critical))
  cat(sprintf("High severity: %d\n", report$n_high))
  cat(sprintf("Affected analyses: %d\n\n", report$n_affected_analyses))

  cat("By Severity:\n")
  print(by_severity, row.names = FALSE)
  cat("\n")

  cat("By Method:\n")
  print(by_method, row.names = FALSE)
  cat("\n")

  cat("By Issue Type:\n")
  print(by_type, row.names = FALSE)
  cat("\n")

  # Save if path provided
  if (!is.null(output_path)) {
    utils::write.csv(issues, file.path(output_path, "all_issues.csv"), row.names = FALSE)
    utils::write.csv(by_severity, file.path(output_path, "issues_by_severity.csv"), row.names = FALSE)
    utils::write.csv(by_method, file.path(output_path, "issues_by_method.csv"), row.names = FALSE)
    utils::write.csv(by_type, file.path(output_path, "issues_by_type.csv"), row.names = FALSE)
    cat(sprintf("Report saved to: %s\n", output_path))
  }

  invisible(report)
}

#' Save Issue Database to File
#'
#' @param file_path Path to save the .rda file
#' @return Invisible TRUE
#' @export
save_issue_database <- function(file_path = "issues_database.rda") {
  save(.issue_database$issues, file = file_path)
  invisible(TRUE)
}

#' Load Issue Database from File
#'
#' @param file_path Path to the .rda file
#' @return Invisible TRUE
#' @export
load_issue_database <- function(file_path = "issues_database.rda") {
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(invisible(FALSE))
  }
  load(file_path)
  # The loaded object is named `.issue_database$issues`
  invisible(TRUE)
}

#' Batch Log Issues from Test Results
#'
#' Process multiple test results and log all issues.
#'
#' @param test_results List of test results from run_dichotomous_tests()
#' @return Summary of logged issues
#' @export
batch_log_issues <- function(test_results) {

  n_total <- 0
  n_logged <- 0

  for (aid in names(test_results$results)) {
    result <- test_results$results[[aid]]

    for (method in names(result$test_results$results)) {
      method_result <- result$test_results$results[[method]]
      n_total <- n_total + 1

      n <- log_test_result(
        method_result,
        result$review_id,
        result$analysis_id,
        NULL  # df not available in stored results
      )
      n_logged <- n_logged + n
    }
  }

  list(
    total_tests = n_total,
    issues_logged = n_logged,
    issue_rate = n_logged / n_total
  )
}
