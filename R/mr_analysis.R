#' Core Mendelian Randomization Analysis Functions
#'
#' This module contains the main functions for performing two-sample MR analysis
#' including data processing, harmonization, and multiple MR methods.
#'
#' @author Xuejiao Hou
#' @export

# Load required libraries
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

#' Perform comprehensive two-sample Mendelian randomization analysis
#'
#' @description Conducts MR analysis using multiple methods including IVW, 
#'              MR Egger, weighted median, and mode-based estimates with 
#'              comprehensive sensitivity testing
#'
#' @param exposure_data Data frame with SNP-exposure associations
#' @param outcome_data Data frame with SNP-outcome associations  
#' @param clump_r2 LD clumping R² threshold (default: 0.001)
#' @param clump_kb LD clumping distance in kb (default: 10000)
#' @param pval_threshold P-value threshold for instrument selection (default: 5e-8)
#' @param f_stat_threshold F-statistic threshold for weak instrument filtering (default: 10)
#' @param harmonise_action Action for strand ambiguous SNPs (default: 2)
#' @param perform_steiger Perform Steiger directionality test (default: TRUE)
#' @param output_dir Directory to save results (default: NULL)
#' @param create_plots Generate diagnostic plots (default: TRUE)
#'
#' @return List containing:
#'   - mr_results: Main MR analysis results with OR and 95% CI
#'   - harmonised_data: Harmonized exposure-outcome dataset
#'   - sensitivity_tests: Pleiotropy and heterogeneity test results
#'   - diagnostic_plots: ggplot objects for visualization
#'   - instrument_stats: Instrument strength statistics
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(example_exposure_data)
#' data(example_outcome_data)
#' 
#' # Run MR analysis
#' results <- perform_mr_analysis(
#'   exposure_data = example_exposure_data,
#'   outcome_data = example_outcome_data
#' )
#' 
#' # View results
#' print(results$mr_results)
#' }
#'
#' @export
perform_mr_analysis <- function(exposure_data, 
                                outcome_data,
                                clump_r2 = 0.001,
                                clump_kb = 10000,
                                pval_threshold = 5e-8,
                                f_stat_threshold = 10,
                                harmonise_action = 2,
                                perform_steiger = TRUE,
                                output_dir = NULL,
                                create_plots = TRUE) {
  
  # Input validation
  if (!is.data.frame(exposure_data) || !is.data.frame(outcome_data)) {
    stop("exposure_data and outcome_data must be data frames")
  }
  
  # Check required columns
  required_cols <- c("SNP", "beta", "se", "pval", "effect_allele", "other_allele")
  if (!all(required_cols %in% names(exposure_data))) {
    stop("exposure_data missing required columns: ", 
         paste(setdiff(required_cols, names(exposure_data)), collapse = ", "))
  }
  
  cat("Starting MR analysis...\n")
  cat("Initial exposure SNPs:", nrow(exposure_data), "\n")
  cat("Initial outcome SNPs:", nrow(outcome_data), "\n")
  
  # Filter exposure data by p-value threshold
  exposure_filtered <- exposure_data[exposure_data$pval < pval_threshold, ]
  cat("Exposure SNPs after p-value filtering:", nrow(exposure_filtered), "\n")
  
  # Harmonise exposure and outcome data
  cat("Harmonizing exposure and outcome data...\n")
  harmonised_data <- harmonise_data(
    exposure_dat = exposure_filtered,
    outcome_dat = outcome_data,
    action = harmonise_action
  )
  
  cat("SNPs after harmonization:", nrow(harmonised_data), "\n")
  
  if (nrow(harmonised_data) == 0) {
    stop("No SNPs remain after harmonization")
  }
  
  # Calculate instrument strength statistics
  cat("Calculating instrument strength statistics...\n")
  harmonised_data <- calculate_instrument_strength(harmonised_data)
  
  # Filter weak instruments
  strong_instruments <- harmonised_data[harmonised_data$f_statistic > f_stat_threshold, ]
  cat("Strong instruments (F >", f_stat_threshold, "):", nrow(strong_instruments), "\n")
  
  if (nrow(strong_instruments) < 3) {
    warning("Fewer than 3 strong instruments available. Results may be unreliable.")
    final_data <- harmonised_data
  } else {
    final_data <- strong_instruments
  }
  
  # Perform Steiger directionality test if requested
  if (perform_steiger) {
    cat("Performing Steiger directionality test...\n")
    final_data <- steiger_filtering(final_data)
    cat("SNPs passing Steiger test:", sum(final_data$steiger_dir), "\n")
  }
  
  # Perform main MR analysis
  cat("Conducting MR analysis...\n")
  mr_results <- mr(final_data)
  
  # Generate odds ratios
  mr_results_or <- generate_odds_ratios(mr_results)
  
  # Perform sensitivity tests
  cat("Conducting sensitivity tests...\n")
  sensitivity_results <- perform_sensitivity_tests(final_data)
  
  # Create diagnostic plots if requested
  plots <- NULL
  if (create_plots) {
    cat("Generating diagnostic plots...\n")
    plots <- create_mr_plots(final_data, mr_results)
  }
  
  # Compile results
  results <- list(
    mr_results = mr_results_or,
    harmonised_data = final_data,
    sensitivity_tests = sensitivity_results,
    diagnostic_plots = plots,
    instrument_stats = list(
      n_total_snps = nrow(exposure_data),
      n_significant_snps = nrow(exposure_filtered),
      n_harmonised_snps = nrow(harmonised_data),
      n_strong_instruments = nrow(strong_instruments),
      n_final_instruments = nrow(final_data),
      mean_f_statistic = mean(final_data$f_statistic, na.rm = TRUE),
      total_variance_explained = sum(final_data$r2_exposure, na.rm = TRUE)
    )
  )
  
  # Save results if output directory specified
  if (!is.null(output_dir)) {
    save_mr_results(results, output_dir)
  }
  
  cat("MR analysis completed successfully!\n")
  return(results)
}

#' Calculate instrument strength statistics
#'
#' @description Calculates F-statistics and R² values for MR instruments
#'
#' @param harmonised_data Harmonized exposure-outcome data
#' 
#' @return Input data with added columns for F-statistics and R²
calculate_instrument_strength <- function(harmonised_data) {
  
  # Check if EAF is available for more accurate R² calculation
  if (all(is.na(harmonised_data$eaf.exposure))) {
    # Calculate R² without EAF (less accurate)
    harmonised_data$r2_exposure <- (2 * (harmonised_data$beta.exposure^2)) /
      (2 * (harmonised_data$beta.exposure^2) + 
       2 * harmonised_data$samplesize.exposure * harmonised_data$se.exposure^2)
  } else {
    # Calculate R² with EAF (more accurate)
    harmonised_data$r2_exposure <- (2 * (harmonised_data$beta.exposure^2) * 
                                   harmonised_data$eaf.exposure * 
                                   (1 - harmonised_data$eaf.exposure)) /
      (2 * (harmonised_data$beta.exposure^2) * 
       harmonised_data$eaf.exposure * (1 - harmonised_data$eaf.exposure) + 
       2 * harmonised_data$samplesize.exposure * harmonised_data$eaf.exposure * 
       (1 - harmonised_data$eaf.exposure) * harmonised_data$se.exposure^2)
  }
  
  # Calculate F-statistics
  harmonised_data$f_statistic <- harmonised_data$r2_exposure * 
    (harmonised_data$samplesize.exposure - 2) / (1 - harmonised_data$r2_exposure)
  
  # Calculate mean F-statistic for the instrument set
  harmonised_data$mean_f_statistic <- mean(harmonised_data$f_statistic, na.rm = TRUE)
  
  return(harmonised_data)
}

#' Perform comprehensive sensitivity tests
#'
#' @description Conducts pleiotropy, heterogeneity, and outlier tests
#'
#' @param harmonised_data Harmonized exposure-outcome data
#'
#' @return List containing sensitivity test results
perform_sensitivity_tests <- function(harmonised_data) {
  
  results <- list()
  
  # MR-Egger intercept test (pleiotropy)
  tryCatch({
    results$pleiotropy <- mr_pleiotropy_test(harmonised_data)
  }, error = function(e) {
    warning("Pleiotropy test failed: ", e$message)
    results$pleiotropy <- NULL
  })
  
  # Heterogeneity tests
  tryCatch({
    heterogeneity <- mr_heterogeneity(harmonised_data)
    # Calculate I² statistic
    heterogeneity$I2 <- pmax(0, (heterogeneity$Q - heterogeneity$Q_df) / heterogeneity$Q)
    results$heterogeneity <- heterogeneity
  }, error = function(e) {
    warning("Heterogeneity test failed: ", e$message)
    results$heterogeneity <- NULL
  })
  
  # Leave-one-out analysis
  tryCatch({
    results$leave_one_out <- mr_leaveoneout(harmonised_data)
  }, error = function(e) {
    warning("Leave-one-out analysis failed: ", e$message)
    results$leave_one_out <- NULL
  })
  
  # Single SNP analysis
  tryCatch({
    results$single_snp <- mr_singlesnp(harmonised_data)
  }, error = function(e) {
    warning("Single SNP analysis failed: ", e$message)
    results$single_snp <- NULL
  })
  
  # MR-PRESSO (if enough SNPs)
  if (nrow(harmonised_data) >= 4) {
    tryCatch({
      results$mr_presso <- run_mr_presso(harmonised_data, NbDistribution = 1000)
    }, error = function(e) {
      warning("MR-PRESSO failed: ", e$message)
      results$mr_presso <- NULL
    })
  }
  
  return(results)
}

#' Create diagnostic plots for MR analysis
#'
#' @description Generates scatter plots, forest plots, and funnel plots
#'
#' @param harmonised_data Harmonized exposure-outcome data
#' @param mr_results MR analysis results
#'
#' @return List of ggplot objects
create_mr_plots <- function(harmonised_data, mr_results) {
  
  plots <- list()
  
  # Scatter plot
  tryCatch({
    plots$scatter <- mr_scatter_plot(mr_results, harmonised_data)[[1]]
  }, error = function(e) {
    warning("Scatter plot creation failed: ", e$message)
  })
  
  # Forest plot
  tryCatch({
    single_snp_results <- mr_singlesnp(harmonised_data)
    plots$forest <- mr_forest_plot(single_snp_results)[[1]]
  }, error = function(e) {
    warning("Forest plot creation failed: ", e$message)
  })
  
  # Leave-one-out plot
  tryCatch({
    loo_results <- mr_leaveoneout(harmonised_data)
    plots$leave_one_out <- mr_leaveoneout_plot(loo_results)[[1]]
  }, error = function(e) {
    warning("Leave-one-out plot creation failed: ", e$message)
  })
  
  # Funnel plot
  tryCatch({
    single_snp_results <- mr_singlesnp(harmonised_data)
    plots$funnel <- mr_funnel_plot(single_snp_results)[[1]]
  }, error = function(e) {
    warning("Funnel plot creation failed: ", e$message)
  })
  
  return(plots)
}

#' Save MR analysis results
#'
#' @description Saves all MR results to specified directory
#'
#' @param results MR analysis results list
#' @param output_dir Output directory path
save_mr_results <- function(results, output_dir) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save main results
  write.table(results$mr_results, 
              file.path(output_dir, "mr_results.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save harmonised data
  write.table(results$harmonised_data,
              file.path(output_dir, "harmonised_data.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save sensitivity tests
  if (!is.null(results$sensitivity_tests$pleiotropy)) {
    write.table(results$sensitivity_tests$pleiotropy,
                file.path(output_dir, "pleiotropy_test.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  if (!is.null(results$sensitivity_tests$heterogeneity)) {
    write.table(results$sensitivity_tests$heterogeneity,
                file.path(output_dir, "heterogeneity_test.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Save plots
  if (!is.null(results$diagnostic_plots)) {
    for (plot_name in names(results$diagnostic_plots)) {
      if (!is.null(results$diagnostic_plots[[plot_name]])) {
        ggsave(file.path(output_dir, paste0(plot_name, "_plot.pdf")),
               results$diagnostic_plots[[plot_name]],
               width = 8, height = 8)
      }
    }
  }
  
  # Save instrument statistics
  cat("=== Instrument Strength Statistics ===\n",
      file = file.path(output_dir, "instrument_stats.txt"))
  cat(paste(names(results$instrument_stats), results$instrument_stats, 
            sep = ": ", collapse = "\n"),
      file = file.path(output_dir, "instrument_stats.txt"), append = TRUE)
  
  cat("Results saved to:", output_dir, "\n")
}

#' Load and format GWAS summary statistics
#'
#' @description Standardized function to load GWAS data with consistent column naming
#'
#' @param file_path Path to GWAS summary statistics file
#' @param data_type Type of data ("exposure" or "outcome")
#' @param snp_col Column name for SNP identifier
#' @param beta_col Column name for effect size
#' @param se_col Column name for standard error
#' @param pval_col Column name for p-value
#' @param effect_allele_col Column name for effect allele
#' @param other_allele_col Column name for other allele
#' @param eaf_col Column name for effect allele frequency (optional)
#' @param samplesize_col Column name for sample size (optional)
#' @param separator File separator (default: "\t")
#'
#' @return Formatted data frame ready for MR analysis
#'
#' @export
load_gwas_data <- function(file_path,
                          data_type = "exposure",
                          snp_col = "SNP",
                          beta_col = "beta",
                          se_col = "se", 
                          pval_col = "pval",
                          effect_allele_col = "effect_allele",
                          other_allele_col = "other_allele",
                          eaf_col = "eaf",
                          samplesize_col = "samplesize",
                          separator = "\t") {
  
  # Read data
  if (tools::file_ext(file_path) %in% c("gz", "bz2")) {
    data <- fread(file_path, sep = separator)
  } else {
    data <- fread(file_path, sep = separator)
  }
  
  # Format data based on type
  if (data_type == "exposure") {
    formatted_data <- format_data(data,
                                 type = "exposure",
                                 snp_col = snp_col,
                                 beta_col = beta_col,
                                 se_col = se_col,
                                 pval_col = pval_col,
                                 effect_allele_col = effect_allele_col,
                                 other_allele_col = other_allele_col,
                                 eaf_col = eaf_col,
                                 samplesize_col = samplesize_col)
  } else {
    formatted_data <- format_data(data,
                                 type = "outcome",
                                 snp_col = snp_col,
                                 beta_col = beta_col,
                                 se_col = se_col,
                                 pval_col = pval_col,
                                 effect_allele_col = effect_allele_col,
                                 other_allele_col = other_allele_col,
                                 eaf_col = eaf_col,
                                 samplesize_col = samplesize_col)
  }
  
  cat("Loaded", nrow(formatted_data), "SNPs from", basename(file_path), "\n")
  return(formatted_data)
}

#' Validate input data format
#'
#' @description Checks that input data has required columns and valid values
#'
#' @param data Input data frame
#' @param data_type Type of data ("exposure" or "outcome")
#'
#' @return TRUE if valid, otherwise stops with error message
#' @export
check_data_format <- function(data, data_type = "exposure") {
  
  required_cols <- c("SNP", "beta", "se", "pval", "effect_allele", "other_allele")
  
  # Check required columns exist
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for missing values in critical columns
  critical_cols <- c("SNP", "beta", "se", "pval")
  for (col in critical_cols) {
    if (any(is.na(data[[col]]))) {
      warning("Missing values found in ", col, " column")
    }
  }
  
  # Check numeric columns are numeric
  numeric_cols <- c("beta", "se", "pval")
  for (col in numeric_cols) {
    if (!is.numeric(data[[col]])) {
      stop("Column ", col, " must be numeric")
    }
  }
  
  # Check p-values are valid
  if (any(data$pval < 0 | data$pval > 1, na.rm = TRUE)) {
    stop("P-values must be between 0 and 1")
  }
  
  # Check standard errors are positive
  if (any(data$se <= 0, na.rm = TRUE)) {
    stop("Standard errors must be positive")
  }
  
  cat("Data format validation passed for", nrow(data), data_type, "SNPs\n")
  return(TRUE)
}
