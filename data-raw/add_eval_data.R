# MLMResearch Evaluation Dataset (Simulated Results for 50 Meta-Analyses)
# This mimics the output of run_dichotomous_tests() for demonstration
set.seed(42)

mlm_eval_results <- data.frame(
  analysis_id = paste0("Analysis_", 1:50),
  k_studies = sample(5:60, 50, replace = TRUE),
  n_effects = sample(10:200, 50, replace = TRUE),
  
  # REML Results
  reml_estimate = rnorm(50, 0.2, 0.1),
  reml_se = runif(50, 0.05, 0.15),
  reml_tau2 = pmax(0, rnorm(50, 0.05, 0.08)), # Some will collapse to 0
  reml_omega2 = pmax(0, rnorm(50, 0.02, 0.04)), # Many will collapse to 0
  reml_converged = sample(c(TRUE, FALSE), 50, replace = TRUE, prob = c(0.9, 0.1)),
  
  # RVE-CR1 Results
  rve_cr1_se = runif(50, 0.06, 0.18),
  rve_cr1_df = runif(50, 2, 45),
  
  # RVE-CR2 (Small-sample adjusted)
  rve_cr2_se = runif(50, 0.08, 0.22),
  rve_cr2_df = runif(50, 3, 50)
)

# Flag issues based on the data
mlm_eval_results$variance_collapse <- mlm_eval_results$reml_tau2 == 0 | mlm_eval_results$reml_omega2 == 0
mlm_eval_results$small_sample <- mlm_eval_results$k_studies < 40

# Save to data directory
if(!dir.exists("data")) dir.create("data")
save(mlm_eval_results, file = "data/mlm_eval_results.rda")
cat("✓ mlm_eval_results.rda created successfully.
")
