#' Multilevel Meta-Analysis Methods Registry
#'
#' A comprehensive registry of modern multilevel meta-analysis methods
#' with their R package implementations, parameters, and known issues.
#'
#' @return A data.frame containing method information
#' @export
#' @examples
#' registry <- get_methods_registry()
#' print(registry)
get_methods_registry <- function() {

  methods <- data.frame(
    category = c(
      # Frequentist Methods
      rep("Frequentist", 4),
      # Robust Variance Estimation
      rep("RVE", 3),
      # Network Meta-Analysis
      rep("Network", 2),
      # Bayesian
      rep("Bayesian", 2)
    ),
    method_name = c(
      # Frequentist
      "REML (Restricted Maximum Likelihood)",
      "ML (Maximum Likelihood)",
      "DL (DerSimonian-Laird)",
      "HE (Hedges estimator)",
      # RVE
      "CR2 (Bias-reduced linearization)",
      "CR1 (HC3 sandwich)",
      "CR0 (White sandwich)",
      # Network
      "Netmeta (frequentist network)",
      "Netmeta (fixed effect network)",
      # Bayesian
      "Hierarchical (brms/Stan)",
      "Hierarchical (rstanarm)"
    ),
    r_function = c(
      "metafor::rma.mv(method='REML')",
      "metafor::rma.mv(method='ML')",
      "metafor::rma.mv(method='DL')",
      "metafor::rma.mv(method='HE')",
      "clubSandwich::vcovCR(type='CR2')",
      "clubSandwich::vcovCR(type='CR1')",
      "clubSandwich::vcovCR(type='CR0')",
      "netmeta::netmeta(comb.fixed=FALSE)",
      "netmeta::netmeta(comb.fixed=TRUE)",
      "brms::brm()",
      "rstanarm::stan_glmer()"
    ),
    package = c(
      "metafor", "metafor", "metafor", "metafor",
      "clubSandwich", "clubSandwich", "clubSandwich",
      "netmeta", "netmeta",
      "brms", "rstanarm"
    ),
    primary_parameters = c(
      "method='REML', random structure",
      "method='ML', random structure",
      "method='DL', random structure",
      "method='HE', random structure",
      "type='CR2', cluster",
      "type='CR1', cluster",
      "type='CR0', cluster",
      "comb.fixed, smallcorrection",
      "comb.fixed=TRUE",
      "prior, iter, chains",
      "prior, iter, chains"
    ),
    known_issues = c(
      "Small sample bias, negative variance possible",
      "Underestimates heterogeneity",
      "Can produce negative variance, less accurate with few studies",
      "Less common, limited validation",
      "Requires 10+ clusters, conservative with few clusters",
      "Less accurate than CR2 for small samples",
      "No small-sample correction, liberal tests",
      "Assumes transitivity, inconsistency detection difficult",
      "Stronger transitivity assumption",
      "High prior sensitivity with small samples, slow",
      "Less flexible than brms, prior sensitivity"
    ),
    recommended_n_clusters = c(
      "5+ for reasonable performance",
      "10+ preferred",
      "10+ preferred",
      "10+ preferred",
      "20+ ideal, 10+ minimum",
      "20+ ideal",
      "30+ ideal",
      "Depends on network structure",
      "Depends on network structure",
      "Any (but priors matter with few clusters)",
      "Any (but priors matter with few clusters)"
    ),
    recommended_outcomes = c(
      "All (logOR, SMD, MD, GIV)",
      "All (logOR, SMD, MD, GIV)",
      "All (logOR, SMD, MD, GIV)",
      "All (logOR, SMD, MD, GIV)",
      "All (especially with dependent effects)",
      "All (especially with dependent effects)",
      "All (especially with dependent effects)",
      "Multiple treatments, logOR preferred",
      "Multiple treatments, logOR preferred",
      "All (logOR, SMD, MD, GIV)",
      "All (logOR, SMD, MD, GIV)"
    ),
    computational_speed = c(
      "Fast",
      "Fast",
      "Fast",
      "Fast",
      "Medium",
      "Medium",
      "Medium",
      "Fast-Medium",
      "Fast",
      "Slow",
      "Medium-Slow"
    ),
    stringsAsFactors = FALSE
  )

  # Add method IDs for easy reference
  methods$method_id <- paste0("m", sprintf("%02d", seq_len(nrow(methods))))
  methods$priority <- ifelse(methods$package %in% c("metafor", "clubSandwich"), "High", "Medium")

  # Reorder columns
  methods <- methods[, c("method_id", "category", "method_name", "r_function",
                         "package", "primary_parameters", "known_issues",
                         "recommended_n_clusters", "recommended_outcomes",
                         "computational_speed", "priority")]

  return(methods)
}

#' Get method recommendations based on data characteristics
#'
#' @param n_effects Number of effect sizes in the analysis
#' @param n_clusters Number of independent clusters
#' @param outcome_type Type of outcome ("DICH", "CONT", "GENSUM")
#' @param has_dependencies Whether there are dependent effect sizes
#' @return Character vector of recommended methods
#' @export
get_method_recommendations <- function(n_effects, n_clusters, outcome_type = "DICH",
                                       has_dependencies = TRUE) {

  registry <- get_methods_registry()

  # Filter by outcome compatibility
  compatible <- registry

  # Small sample rules
  if (n_clusters < 10) {
    # Avoid CR2 with very few clusters
    compatible <- compatible[compatible$method_name != "CR2 (Bias-reduced linearization)" |
                             !grepl("CR2", compatible$known_issues), ]
    message("Warning: Fewer than 10 clusters. RVE methods may be unreliable.")
  }

  # Computational efficiency for large datasets
  if (n_effects > 500) {
    compatible <- compatible[compatible$computational_speed != "Slow", ]
  }

  # Dependency handling
  if (has_dependencies) {
    # Prioritize RVE methods
    compatible$priority[grepl("RVE", compatible$category)] <- "High"
  }

  # Return recommended methods
  recommended <- compatible[compatible$priority == "High" |
                            compatible$method_name %in% c("REML (Restricted Maximum Likelihood)",
                                                         "CR2 (Bias-reduced linearization)"), ]

  return(recommended$method_name)
}

#' Get method details by ID
#'
#' @param method_id Method identifier (e.g., "m01", "m05")
#' @return List with method details
#' @export
get_method_details <- function(method_id) {
  registry <- get_methods_registry()
  method <- registry[registry$method_id == method_id, ]

  if (nrow(method) == 0) {
    stop("Method ID not found: ", method_id)
  }

  list(
    id = method$method_id,
    name = method$method_name,
    category = method$category,
    function_call = method$r_function,
    package = method$package,
    parameters = method$primary_parameters,
    issues = method$known_issues,
    clusters_required = method$recommended_n_clusters,
    outcomes = method$recommended_outcomes,
    speed = method$computational_speed
  )
}

#' Print method comparison table
#'
#' @param subset_methods Character vector of method names to compare (NULL for all)
#' @export
print_method_comparison <- function(subset_methods = NULL) {
  registry <- get_methods_registry()

  if (!is.null(subset_methods)) {
    registry <- registry[registry$method_name %in% subset_methods, ]
  }

  cat("\n=== Multilevel Meta-Analysis Methods Comparison ===\n\n")
  cat(sprintf("%-4s | %-20s | %-15s | %-10s | %-8s\n",
              "ID", "Method", "Package", "Speed", "Priority"))
  cat(sprintf("%s\n", paste(rep("-", 75), collapse = "")))

  for (i in seq_len(nrow(registry))) {
    cat(sprintf("%-4s | %-20s | %-15s | %-10s | %-8s\n",
                registry$method_id[i],
                substr(registry$method_name[i], 1, 20),
                registry$package[i],
                registry$computational_speed[i],
                registry$priority[i]))
  }
  cat("\n")
}

# Internal: Helper to check if package is available
check_package_available <- function(package_name) {
  requireNamespace(package_name, quietly = TRUE)
}

# Internal: Helper to get install command for packages
get_install_commands <- function() {
  registry <- get_methods_registry()
  unique_packages <- unique(registry$package)

  cat("Install required packages:\n\n")
  cat("# CRAN packages\n")
  cat("install.packages(c(", paste0('"', unique_packages, '"'), "))\n\n")
}
