#' Mediation Analysis Functions for Multi-omics Mendelian Randomization
#'
#' This module contains functions for performing mediation analysis to test
#' causal pathways through intermediate variables (e.g., Gene → Brain Imaging → MDD)
#'
#' @author Xuejiao Hou
#' @export

# Load required libraries
suppressPackageStartupMessages({
  library(boot)
  library(dplyr)
  library(data.table)
})

#' Calculate mediation effects using Sobel test and bootstrap confidence intervals
#'
#' @description Performs mediation analysis to decompose total effects into
#'              direct and indirect (mediated) effects with statistical testing
#'
#' @param exposure_mediator_results MR results for exposure → mediator pathway
#' @param mediator_outcome_results MR results for mediator → outcome pathway  
#' @param exposure_outcome_results MR results for exposure → outcome pathway (total effect)
#' @param bootstrap_n Number of bootstrap samples (default: 10000)
#' @param confidence_level Confidence level for intervals (default: 0.95)
#' @param set_seed Random seed for reproducibility (default: 123)
#'
#' @return List containing:
#'   - mediation_effects: Data frame with effect estimates
#'   - sobel_test: Sobel test results 
#'   - bootstrap_results: Bootstrap confidence intervals
#'   - proportion_mediated: Proportion of total effect that is mediated
#'
#' @examples
#' \dontrun{
#' # Calculate mediation effects
#' med_results <- calculate_mediation_effect(
#'   exposure_mediator_results = gene_protein_results,
#'   mediator_outcome_results = protein_mdd_results,
#'   exposure_outcome_results = gene_mdd_results
#' )
#' 
#' # View results
#' print(med_results$mediation_effects)
#' }
#'
#' @export
calculate_mediation_effect <- function(exposure_mediator_results,
                                      mediator_outcome_results,
                                      exposure_outcome_results,
                                      bootstrap_n = 10000,
                                      confidence_level = 0.95,
                                      set_seed = 123) {
  
  # Set seed for reproducibility
  set.seed(set_seed)
  
  # Extract effect estimates (use IVW results)
  beta_xm <- extract_ivw_estimate(exposure_mediator_results, "beta")
  se_xm <- extract_ivw_estimate(exposure_mediator_results, "se")
  
  beta_my <- extract_ivw_estimate(mediator_outcome_results, "beta")
  se_my <- extract_ivw_estimate(mediator_outcome_results, "se")
  
  beta_total <- extract_ivw_estimate(exposure_outcome_results, "beta")
  se_total <- extract_ivw_estimate(exposure_outcome_results, "se")
  
  cat("Calculating mediation effects...\n")
  cat("Exposure → Mediator: β =", round(beta_xm, 4), "± SE =", round(se_xm, 4), "\n")
  cat("Mediator → Outcome: β =", round(beta_my, 4), "± SE =", round(se_my, 4), "\n")
  cat("Total Effect: β =", round(beta_total, 4), "± SE =", round(se_total, 4), "\n")
  
  # Calculate mediation effects
  effects <- calculate_effects(beta_xm, se_xm, beta_my, se_my, beta_total)
  
  # Perform Sobel test
  sobel_results <- perform_sobel_test(beta_xm, se_xm, beta_my, se_my)
  
  # Perform bootstrap analysis
  cat("Performing bootstrap analysis with", bootstrap_n, "samples...\n")
  bootstrap_results <- perform_mediation_bootstrap(
    beta_xm, se_xm, beta_my, se_my, bootstrap_n, confidence_level
  )
  
  # Calculate proportion mediated
  prop_mediated <- calculate_proportion_mediated(effects, bootstrap_results, confidence_level)
  
  # Compile results
  results <- list(
    mediation_effects = effects,
    sobel_test = sobel_results,
    bootstrap_results = bootstrap_results,
    proportion_mediated = prop_mediated,
    input_parameters = list(
      bootstrap_n = bootstrap_n,
      confidence_level = confidence_level,
      seed = set_seed
    )
  )
  
  cat("Mediation analysis completed!\n")
  return(results)
}

#' Extract IVW estimate from MR results
#'
#' @description Helper function to extract beta or SE from IVW method
#'
#' @param mr_results MR analysis results
#' @param parameter Parameter to extract ("beta", "se", or "pval")
#'
#' @return Numeric value
extract_ivw_estimate <- function(mr_results, parameter) {
  
  # Handle both TwoSampleMR and custom result formats
  if ("method" %in% names(mr_results)) {
    ivw_row <- mr_results[mr_results$method == "Inverse variance weighted", ]
    if (nrow(ivw_row) == 0) {
      ivw_row <- mr_results[grepl("IVW", mr_results$method), ]
    }
  } else if (is.list(mr_results) && "mr_results" %in% names(mr_results)) {
    mr_data <- mr_results$mr_results
    ivw_row <- mr_data[mr_data$method == "Inverse variance weighted", ]
  } else {
    stop("Unrecognized MR results format")
  }
  
  if (nrow(ivw_row) == 0) {
    stop("IVW results not found in MR results")
  }
  
  # Map parameter to column name
  param_map <- c(
    "beta" = c("b", "beta", "effect"),
    "se" = c("se", "std_error"),
    "pval" = c("pval", "p", "p_value")
  )
  
  possible_cols <- param_map[[parameter]]
  col_name <- intersect(possible_cols, names(ivw_row))[1]
  
  if (is.na(col_name)) {
    stop("Could not find column for parameter: ", parameter)
  }
  
  return(as.numeric(ivw_row[[col_name]]))
}

#' Calculate direct, indirect, and total effects
#'
#' @description Calculates mediation effect decomposition
#'
#' @param beta_xm Exposure → mediator effect
#' @param se_xm Exposure → mediator standard error  
#' @param beta_my Mediator → outcome effect
#' @param se_my Mediator → outcome standard error
#' @param beta_total Total exposure → outcome effect
#'
#' @return Data frame with effect estimates
calculate_effects <- function(beta_xm, se_xm, beta_my, se_my, beta_total) {
  
  # Indirect effect (mediated effect)
  indirect_effect <- beta_xm * beta_my
  
  # Standard error of indirect effect (delta method)
  se_indirect <- sqrt((beta_xm^2 * se_my^2) + (beta_my^2 * se_xm^2))
  
  # Direct effect
  direct_effect <- beta_total - indirect_effect
  
  # Confidence intervals for indirect effect
  ci_lower_indirect <- indirect_effect - 1.96 * se_indirect
  ci_upper_indirect <- indirect_effect + 1.96 * se_indirect
  
  effects <- data.frame(
    effect_type = c("Total", "Direct", "Indirect"),
    beta = c(beta_total, direct_effect, indirect_effect),
    se = c(NA, NA, se_indirect),  # SE only calculated for indirect effect
    ci_lower = c(NA, NA, ci_lower_indirect),
    ci_upper = c(NA, NA, ci_upper_indirect),
    stringsAsFactors = FALSE
  )
  
  return(effects)
}

#' Perform Sobel test for mediation
#'
#' @description Tests significance of indirect effect using Sobel test
#'
#' @param beta_xm Exposure → mediator effect
#' @param se_xm Exposure → mediator standard error
#' @param beta_my Mediator → outcome effect  
#' @param se_my Mediator → outcome standard error
#'
#' @return Data frame with Sobel test results
perform_sobel_test <- function(beta_xm, se_xm, beta_my, se_my) {
  
  # Indirect effect
  indirect_effect <- beta_xm * beta_my
  
  # Standard error (Sobel formula)
  se_sobel <- sqrt((beta_xm^2 * se_my^2) + (beta_my^2 * se_xm^2))
  
  # Z-statistic
  z_statistic <- indirect_effect / se_sobel
  
  # Two-tailed p-value
  p_value <- 2 * (1 - pnorm(abs(z_statistic)))
  
  sobel_results <- data.frame(
    test = "Sobel",
    indirect_effect = indirect_effect,
    se = se_sobel,
    z_statistic = z_statistic,
    p_value = p_value,
    significant = p_value < 0.05,
    stringsAsFactors = FALSE
  )
  
  return(sobel_results)
}

#' Perform bootstrap mediation analysis
#'
#' @description Bootstrap confidence intervals for mediation effects
#'
#' @param beta_xm Exposure → mediator effect
#' @param se_xm Exposure → mediator standard error
#' @param beta_my Mediator → outcome effect
#' @param se_my Mediator → outcome standard error  
#' @param bootstrap_n Number of bootstrap samples
#' @param confidence_level Confidence level
#'
#' @return List with bootstrap results and confidence intervals
perform_mediation_bootstrap <- function(beta_xm, se_xm, beta_my, se_my, 
                                       bootstrap_n, confidence_level) {
  
  # Bootstrap function for mediation effect
  mediation_boot_func <- function(data, indices) {
    # Sample effect estimates from normal distributions
    beta_xm_boot <- rnorm(1, beta_xm, se_xm)
    beta_my_boot <- rnorm(1, beta_my, se_my)
    
    # Calculate indirect effect
    indirect <- beta_xm_boot * beta_my_boot
    return(indirect)
  }
  
  # Create dummy data for boot function
  dummy_data <- data.frame(x = 1:10)
  
  # Perform bootstrap
  boot_results <- boot(
    data = dummy_data,
    statistic = function(data, indices) {
      beta_xm_boot <- rnorm(1, beta_xm, se_xm)
      beta_my_boot <- rnorm(1, beta_my, se_my)
      return(beta_xm_boot * beta_my_boot)
    },
    R = bootstrap_n
  )
  
  # Calculate confidence intervals using different methods
  alpha <- 1 - confidence_level
  
  # Percentile method
  percentile_ci <- quantile(boot_results$t, c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  
  # Bias-corrected and accelerated (BCa) method
  bca_ci <- tryCatch({
    boot.ci(boot_results, type = "bca", conf = confidence_level)$bca[4:5]
  }, error = function(e) {
    warning("BCa CI calculation failed, using percentile method")
    percentile_ci
  })
  
  # Normal approximation
  normal_ci <- tryCatch({
    boot.ci(boot_results, type = "norm", conf = confidence_level)$normal[2:3]
  }, error = function(e) {
    warning("Normal CI calculation failed")
    c(NA, NA)
  })
  
  # Basic bootstrap
  basic_ci <- tryCatch({
    boot.ci(boot_results, type = "basic", conf = confidence_level)$basic[4:5]
  }, error = function(e) {
    warning("Basic CI calculation failed")
    c(NA, NA)
  })
  
  # Bootstrap p-value (proportion of bootstrap samples with opposite sign)
  original_effect <- beta_xm * beta_my
  boot_p_value <- 2 * min(
    mean(boot_results$t >= 0, na.rm = TRUE),
    mean(boot_results$t < 0, na.rm = TRUE)
  )
  
  bootstrap_summary <- data.frame(
    method = c("Percentile", "BCa", "Normal", "Basic"),
    ci_lower = c(percentile_ci[1], bca_ci[1], normal_ci[1], basic_ci[1]),
    ci_upper = c(percentile_ci[2], bca_ci[2], normal_ci[2], basic_ci[2]),
    stringsAsFactors = FALSE
  )
  
  results <- list(
    boot_object = boot_results,
    confidence_intervals = bootstrap_summary,
    bootstrap_p_value = boot_p_value,
    bootstrap_mean = mean(boot_results$t, na.rm = TRUE),
    bootstrap_se = sd(boot_results$t, na.rm = TRUE)
  )
  
  return(results)
}

#' Calculate proportion mediated with confidence intervals
#'
#' @description Calculates what proportion of the total effect is mediated
#'
#' @param effects Effect estimates data frame
#' @param bootstrap_results Bootstrap results
#' @param confidence_level Confidence level
#'
#' @return Data frame with proportion mediated estimates
calculate_proportion_mediated <- function(effects, bootstrap_results, confidence_level) {
  
  # Extract effects
  total_effect <- effects$beta[effects$effect_type == "Total"]
  indirect_effect <- effects$beta[effects$effect_type == "Indirect"]
  
  # Calculate proportion mediated
  prop_mediated <- ifelse(total_effect != 0, indirect_effect / total_effect, NA)
  
  # Use bootstrap samples to calculate CI for proportion
  if (!is.null(bootstrap_results$boot_object)) {
    # Calculate proportion for each bootstrap sample
    # Note: We would need the total effect bootstrap samples too for proper calculation
    # For now, use the point estimate approach
    
    alpha <- 1 - confidence_level
    
    # Simple confidence interval based on bootstrap SE
    if (!is.na(prop_mediated) && !is.null(bootstrap_results$bootstrap_se)) {
      prop_se_approx <- abs(bootstrap_results$bootstrap_se / total_effect)
      prop_ci_lower <- prop_mediated - qnorm(1 - alpha/2) * prop_se_approx
      prop_ci_upper <- prop_mediated + qnorm(1 - alpha/2) * prop_se_approx
    } else {
      prop_ci_lower <- NA
      prop_ci_upper <- NA
    }
  } else {
    prop_ci_lower <- NA
    prop_ci_upper <- NA
  }
  
  prop_results <- data.frame(
    proportion_mediated = prop_mediated,
    ci_lower = prop_ci_lower,
    ci_upper = prop_ci_upper,
    interpretation = case_when(
      is.na(prop_mediated) ~ "Cannot calculate (total effect = 0)",
      abs(prop_mediated) < 0.1 ~ "Minimal mediation (<10%)",
      abs(prop_mediated) < 0.3 ~ "Small mediation (10-30%)",
      abs(prop_mediated) < 0.7 ~ "Medium mediation (30-70%)",
      TRUE ~ "Large mediation (>70%)"
    ),
    stringsAsFactors = FALSE
  )
  
  return(prop_results)
}

#' Batch mediation analysis for multiple mediators
#'
#' @description Performs mediation analysis for multiple potential mediators
#'
#' @param exposure_data Exposure GWAS data
#' @param mediator_list List of mediator GWAS data or file paths
#' @param outcome_data Outcome GWAS data
#' @param output_dir Directory to save results
#' @param bootstrap_n Number of bootstrap samples (default: 5000 for batch)
#' @param parallel Use parallel processing (default: FALSE)
#'
#' @return List of mediation results for each mediator
#'
#' @export
run_batch_mediation_analysis <- function(exposure_data,
                                        mediator_list,
                                        outcome_data,
                                        output_dir,
                                        bootstrap_n = 5000,
                                        parallel = FALSE) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Starting batch mediation analysis...\n")
  cat("Number of mediators:", length(mediator_list), "\n")
  
  # Load required libraries for MR analysis
  library(TwoSampleMR)
  
  # Function to process single mediator
  process_mediator <- function(i, mediator_item) {
    
    cat("Processing mediator", i, "of", length(mediator_list), "\n")
    
    # Load mediator data if it's a file path
    if (is.character(mediator_item) && length(mediator_item) == 1) {
      mediator_data <- load_gwas_data(mediator_item, data_type = "outcome")
    } else {
      mediator_data <- mediator_item
    }
    
    tryCatch({
      # Step 1: Exposure → Mediator
      harmonised_xm <- harmonise_data(exposure_data, mediator_data, action = 2)
      if (nrow(harmonised_xm) < 3) {
        warning("Not enough SNPs for exposure-mediator analysis for mediator ", i)
        return(NULL)
      }
      mr_results_xm <- mr(harmonised_xm)
      
      # Step 2: Mediator → Outcome  
      harmonised_my <- harmonise_data(mediator_data, outcome_data, action = 2)
      if (nrow(harmonised_my) < 3) {
        warning("Not enough SNPs for mediator-outcome analysis for mediator ", i)
        return(NULL)
      }
      mr_results_my <- mr(harmonised_my)
      
      # Step 3: Exposure → Outcome (total effect)
      harmonised_xy <- harmonise_data(exposure_data, outcome_data, action = 2)
      if (nrow(harmonised_xy) < 3) {
        warning("Not enough SNPs for exposure-outcome analysis for mediator ", i)
        return(NULL)
      }
      mr_results_xy <- mr(harmonised_xy)
      
      # Perform mediation analysis
      mediation_results <- calculate_mediation_effect(
        exposure_mediator_results = mr_results_xm,
        mediator_outcome_results = mr_results_my,
        exposure_outcome_results = mr_results_xy,
        bootstrap_n = bootstrap_n
      )
      
      # Add mediator identifier
      mediation_results$mediator_id <- i
      mediation_results$mediator_name <- ifelse(
        is.character(mediator_item), basename(mediator_item), paste("mediator", i)
      )
      
      # Save individual results
      mediator_output_dir <- file.path(output_dir, paste0("mediator_", i))
      save_mediation_results(mediation_results, mediator_output_dir)
      
      return(mediation_results)
      
    }, error = function(e) {
      warning("Mediation analysis failed for mediator ", i, ": ", e$message)
      return(NULL)
    })
  }
  
  # Run analysis
  if (parallel && requireNamespace("foreach", quietly = TRUE)) {
    results <- foreach::foreach(i = seq_along(mediator_list)) %do% {
      process_mediator(i, mediator_list[[i]])
    }
  } else {
    results <- lapply(seq_along(mediator_list), function(i) {
      process_mediator(i, mediator_list[[i]])
    })
  }
  
  # Remove NULL results
  results <- results[!sapply(results, is.null)]
  
  # Create combined summary
  create_mediation_summary(results, output_dir)
  
  cat("Batch mediation analysis completed!\n")
  cat("Successful analyses:", length(results), "/", length(mediator_list), "\n")
  
  return(results)
}

#' Save mediation analysis results
#'
#' @description Saves mediation results to specified directory
#'
#' @param results Mediation results list
#' @param output_dir Output directory path
save_mediation_results <- function(results, output_dir) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save main effects
  write.table(results$mediation_effects,
              file.path(output_dir, "mediation_effects.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save Sobel test results
  write.table(results$sobel_test,
              file.path(output_dir, "sobel_test.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save bootstrap confidence intervals
  write.table(results$bootstrap_results$confidence_intervals,
              file.path(output_dir, "bootstrap_confidence_intervals.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save proportion mediated
  write.table(results$proportion_mediated,
              file.path(output_dir, "proportion_mediated.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save analysis summary
  summary_text <- paste(
    "=== Mediation Analysis Summary ===",
    paste("Bootstrap samples:", results$input_parameters$bootstrap_n),
    paste("Confidence level:", results$input_parameters$confidence_level),
    paste("Random seed:", results$input_parameters$seed),
    "",
    "=== Key Results ===",
    paste("Indirect effect:", round(results$mediation_effects$beta[3], 4)),
    paste("Sobel test p-value:", format(results$sobel_test$p_value, scientific = TRUE)),
    paste("Bootstrap p-value:", format(results$bootstrap_results$bootstrap_p_value, scientific = TRUE)),
    paste("Proportion mediated:", round(results$proportion_mediated$proportion_mediated, 3)),
    sep = "\n"
  )
  
  writeLines(summary_text, file.path(output_dir, "mediation_summary.txt"))
  
  cat("Mediation results saved to:", output_dir, "\n")
}

#' Create combined summary for batch mediation analysis
#'
#' @description Creates summary table for multiple mediation analyses
#'
#' @param results_list List of mediation results
#' @param output_dir Output directory
create_mediation_summary <- function(results_list, output_dir) {
  
  # Extract key results from each analysis
  summary_data <- do.call(rbind, lapply(results_list, function(result) {
    
    indirect_effect <- result$mediation_effects$beta[result$mediation_effects$effect_type == "Indirect"]
    sobel_p <- result$sobel_test$p_value
    bootstrap_p <- result$bootstrap_results$bootstrap_p_value
    prop_mediated <- result$proportion_mediated$proportion_mediated
    
    # Get best bootstrap CI (prefer BCa, fall back to percentile)
    ci_data <- result$bootstrap_results$confidence_intervals
    bca_row <- ci_data[ci_data$method == "BCa", ]
    if (nrow(bca_row) > 0 && !is.na(bca_row$ci_lower)) {
      ci_lower <- bca_row$ci_lower
      ci_upper <- bca_row$ci_upper
      ci_method <- "BCa"
    } else {
      perc_row <- ci_data[ci_data$method == "Percentile", ]
      ci_lower <- perc_row$ci_lower
      ci_upper <- perc_row$ci_upper  
      ci_method <- "Percentile"
    }
    
    data.frame(
      mediator_id = result$mediator_id,
      mediator_name = result$mediator_name,
      indirect_effect = indirect_effect,
      sobel_p_value = sobel_p,
      bootstrap_p_value = bootstrap_p,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      ci_method = ci_method,
      proportion_mediated = prop_mediated,
      significant_sobel = sobel_p < 0.05,
      significant_bootstrap = bootstrap_p < 0.05,
      stringsAsFactors = FALSE
    )
  }))
  
  # Sort by bootstrap p-value
  summary_data <- summary_data[order(summary_data$bootstrap_p_value), ]
  
  # Save summary
  write.table(summary_data,
              file.path(output_dir, "mediation_summary_combined.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Create significant results subset
  significant_results <- summary_data[summary_data$significant_bootstrap, ]
  if (nrow(significant_results) > 0) {
    write.table(significant_results,
                file.path(output_dir, "significant_mediation_results.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  cat("Combined mediation summary saved to:", output_dir, "\n")
}
