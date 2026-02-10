#' GSMR (Generalized Summary-data-based Mendelian Randomization) Analysis Functions
#'
#' This module contains functions for performing GSMR analysis with HEIDI 
#' outlier detection and visualization.
#'
#' @author Xuejiao Hou
#' @export

# Load required libraries
suppressPackageStartupMessages({
  library(gsmr2)
  library(data.table)
  library(dplyr)
})

#' Perform GSMR analysis with HEIDI outlier detection
#'
#' @description Conducts GSMR analysis which is more robust to horizontal 
#'              pleiotropy compared to standard MR methods
#'
#' @param harmonised_data Harmonized exposure-outcome data from TwoSampleMR
#' @param gwas_threshold P-value threshold for instrument selection (default: 5e-8)
#' @param heidi_outlier_flag Enable HEIDI outlier detection (default: TRUE)  
#' @param multi_snps_heidi_thresh HEIDI test threshold for multi-SNP analysis (default: 0.01)
#' @param nsnps_thresh Minimum number of SNPs required (default: 3)
#' @param ld_r2_thresh LD R² threshold for outlier detection (default: 0.05)
#' @param ld_fdr_thresh FDR threshold for LD pruning (default: 0.05)
#' @param n_ref Reference sample size for LD matrix (default: 7703)
#' @param gsmr2_beta GSMR2 beta parameter (default: 1)
#' @param create_plot Generate GSMR scatter plot (default: TRUE)
#' @param output_dir Directory to save results (default: NULL)
#'
#' @return List containing:
#'   - gsmr_results: GSMR analysis results
#'   - used_snps: SNPs used in final analysis after outlier detection
#'   - outlier_snps: SNPs identified as outliers
#'   - gsmr_plot: ggplot object of GSMR scatter plot (if create_plot = TRUE)
#'   - effect_estimates: Effect size estimates with confidence intervals
#'
#' @examples
#' \dontrun{
#' # Run GSMR analysis
#' gsmr_results <- run_gsmr_analysis(harmonised_data)
#' 
#' # View results
#' print(gsmr_results$effect_estimates)
#' }
#'
#' @export
run_gsmr_analysis <- function(harmonised_data,
                              gwas_threshold = 5e-8,
                              heidi_outlier_flag = TRUE,
                              multi_snps_heidi_thresh = 0.01,
                              nsnps_thresh = 3,
                              ld_r2_thresh = 0.05,
                              ld_fdr_thresh = 0.05,
                              n_ref = 7703,
                              gsmr2_beta = 1,
                              create_plot = TRUE,
                              output_dir = NULL) {
  
  # Input validation
  if (!is.data.frame(harmonised_data)) {
    stop("harmonised_data must be a data frame")
  }
  
  if (nrow(harmonised_data) < nsnps_thresh) {
    stop("Not enough SNPs for GSMR analysis. Need at least ", nsnps_thresh, " SNPs")
  }
  
  cat("Starting GSMR analysis...\n")
  cat("Input SNPs:", nrow(harmonised_data), "\n")
  
  # Prepare GSMR data format
  gsmr_data <- prepare_gsmr_data(harmonised_data)
  
  # Create identity LD matrix (assumes independent SNPs after clumping)
  ldrho <- diag(nrow(gsmr_data))
  colnames(ldrho) <- rownames(ldrho) <- as.character(gsmr_data$SNP)
  snp_coeff_id <- snpid <- as.character(gsmr_data$SNP)
  
  # Extract effect sizes and standard errors
  bzx <- gsmr_data$bzx
  bzx_se <- gsmr_data$bzx_se
  bzx_pval <- gsmr_data$bzx_pval
  bzy <- gsmr_data$bzy
  bzy_se <- gsmr_data$bzy_se
  bzy_pval <- gsmr_data$bzy_pval
  
  cat("Running GSMR with the following parameters:\n")
  cat("  - GWAS threshold:", gwas_threshold, "\n")
  cat("  - HEIDI outlier detection:", heidi_outlier_flag, "\n")
  cat("  - Minimum SNPs:", nsnps_thresh, "\n")
  cat("  - LD R² threshold:", ld_r2_thresh, "\n")
  
  # Perform GSMR analysis
  tryCatch({
    gsmr_results <- gsmr(bzx = bzx,
                        bzx_se = bzx_se,
                        bzx_pval = bzx_pval,
                        bzy = bzy,
                        bzy_se = bzy_se,
                        bzy_pval = bzy_pval,
                        ldrho = ldrho,
                        snp_coeff_id = snp_coeff_id,
                        n_ref = n_ref,
                        heidi_outlier_flag = heidi_outlier_flag,
                        gwas_thresh = gwas_threshold,
                        multi_snps_heidi_thresh = multi_snps_heidi_thresh,
                        nsnps_thresh = nsnps_thresh,
                        ld_r2_thresh = ld_r2_thresh,
                        ld_fdr_thresh = ld_fdr_thresh,
                        gsmr2_beta = gsmr2_beta)
    
    cat("GSMR analysis completed successfully!\n")
    
  }, error = function(e) {
    stop("GSMR analysis failed: ", e$message)
  })
  
  # Extract results and create summary
  filtered_index <- gsmr_results$used_index
  outlier_index <- setdiff(1:length(bzx), filtered_index)
  
  cat("SNPs used in final analysis:", length(filtered_index), "\n")
  cat("SNPs identified as outliers:", length(outlier_index), "\n")
  
  # Calculate effect estimates with confidence intervals
  effect_estimates <- calculate_gsmr_estimates(gsmr_results)
  
  # Create GSMR plot if requested
  gsmr_plot <- NULL
  if (create_plot) {
    cat("Creating GSMR plot...\n")
    gsmr_plot <- create_gsmr_plot(bzx, bzx_se, bzy, bzy_se, 
                                 gsmr_results, 
                                 harmonised_data$exposure[1], 
                                 harmonised_data$outcome[1])
  }
  
  # Compile results
  results <- list(
    gsmr_results = gsmr_results,
    used_snps = gsmr_data[filtered_index, ],
    outlier_snps = if(length(outlier_index) > 0) gsmr_data[outlier_index, ] else NULL,
    gsmr_plot = gsmr_plot,
    effect_estimates = effect_estimates,
    analysis_info = list(
      n_input_snps = nrow(harmonised_data),
      n_used_snps = length(filtered_index),
      n_outlier_snps = length(outlier_index),
      gwas_threshold = gwas_threshold,
      heidi_enabled = heidi_outlier_flag
    )
  )
  
  # Save results if output directory specified
  if (!is.null(output_dir)) {
    save_gsmr_results(results, output_dir)
  }
  
  return(results)
}

#' Prepare data for GSMR analysis
#'
#' @description Converts harmonized MR data to GSMR format
#'
#' @param harmonised_data Harmonized exposure-outcome data
#'
#' @return Data frame in GSMR format
prepare_gsmr_data <- function(harmonised_data) {
  
  gsmr_data <- data.frame(
    SNP = harmonised_data$SNP,
    a1 = harmonised_data$effect_allele.exposure,
    a2 = harmonised_data$other_allele.exposure,
    a1_freq = harmonised_data$eaf.exposure,
    bzx_n = harmonised_data$samplesize.exposure,
    bzx = harmonised_data$beta.exposure,
    bzx_se = harmonised_data$se.exposure,
    bzx_pval = harmonised_data$pval.exposure,
    bzy = harmonised_data$beta.outcome,
    bzy_n = harmonised_data$samplesize.outcome,
    bzy_se = harmonised_data$se.outcome,
    bzy_pval = harmonised_data$pval.outcome,
    stringsAsFactors = FALSE
  )
  
  # Handle missing values
  gsmr_data$a1_freq[is.na(gsmr_data$a1_freq)] <- 0.5  # Default MAF if missing
  
  return(gsmr_data)
}

#' Calculate GSMR effect estimates with confidence intervals
#'
#' @description Calculates effect sizes, ORs, and confidence intervals from GSMR results
#'
#' @param gsmr_results GSMR analysis results object
#'
#' @return Data frame with effect estimates and confidence intervals
calculate_gsmr_estimates <- function(gsmr_results) {
  
  # Extract effect estimate and standard error
  beta <- gsmr_results$bxy
  beta_se <- gsmr_results$bxy_se
  pvalue <- gsmr_results$bxy_pval
  
  # Calculate confidence intervals
  beta_lci <- beta - 1.96 * beta_se
  beta_uci <- beta + 1.96 * beta_se
  
  # Calculate odds ratios and CIs
  or <- exp(beta)
  or_lci <- exp(beta_lci)  
  or_uci <- exp(beta_uci)
  
  # Create results data frame
  estimates <- data.frame(
    method = "GSMR",
    nsnp = length(gsmr_results$used_index),
    beta = beta,
    se = beta_se,
    pval = pvalue,
    beta_lci = beta_lci,
    beta_uci = beta_uci,
    or = or,
    or_lci = or_lci,
    or_uci = or_uci,
    stringsAsFactors = FALSE
  )
  
  return(estimates)
}

#' Create GSMR scatter plot
#'
#' @description Creates a scatter plot showing the GSMR analysis results
#'
#' @param bzx Exposure effect sizes
#' @param bzx_se Exposure standard errors
#' @param bzy Outcome effect sizes  
#' @param bzy_se Outcome standard errors
#' @param gsmr_results GSMR results object
#' @param exposure_name Name of exposure variable
#' @param outcome_name Name of outcome variable
#' 
#' @return ggplot object
create_gsmr_plot <- function(bzx, bzx_se, bzy, bzy_se, gsmr_results, 
                            exposure_name, outcome_name) {
  
  # Get indices of SNPs used in analysis
  filtered_index <- gsmr_results$used_index
  
  # Calculate plot ranges
  vals_x <- c(bzx[filtered_index] - bzx_se[filtered_index], 
              bzx[filtered_index] + bzx_se[filtered_index])
  xmin <- min(vals_x)
  xmax <- max(vals_x)
  
  vals_y <- c(bzy[filtered_index] - bzy_se[filtered_index], 
              bzy[filtered_index] + bzy_se[filtered_index])
  ymin <- min(vals_y)
  ymax <- max(vals_y)
  
  # Create data frame for plotting
  plot_data <- data.frame(
    bzx = bzx[filtered_index],
    bzy = bzy[filtered_index],
    bzx_se = bzx_se[filtered_index],
    bzy_se = bzy_se[filtered_index]
  )
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = bzx, y = bzy)) +
    geom_point(color = "#4189C8", size = 2, alpha = 0.8) +
    geom_errorbar(aes(ymin = bzy - bzy_se, ymax = bzy + bzy_se), 
                  color = "#4189C8", width = 0, alpha = 0.6) +
    geom_errorbarh(aes(xmin = bzx - bzx_se, xmax = bzx + bzx_se), 
                   color = "#4189C8", height = 0, alpha = 0.6) +
    geom_abline(intercept = 0, slope = gsmr_results$bxy, 
                linetype = "dashed", color = "grey40", size = 1) +
    xlim(xmin, xmax) +
    ylim(ymin, ymax) +
    labs(
      x = bquote(.(exposure_name) ~ "(" * italic(b[zx]) * ")"),
      y = bquote(.(outcome_name) ~ "(" * italic(b[zy]) * ")"),
      title = paste("GSMR Analysis:", exposure_name, "→", outcome_name),
      subtitle = paste("β =", round(gsmr_results$bxy, 4), 
                      ", P =", format(gsmr_results$bxy_pval, scientific = TRUE, digits = 3),
                      ", N SNPs =", length(filtered_index))
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      panel.border = element_rect(fill = NA, color = "black", size = 0.5)
    )
  
  return(p)
}

#' Run batch GSMR analysis on multiple exposures
#'
#' @description Performs GSMR analysis on multiple exposure-outcome pairs
#'
#' @param exposure_list List of exposure data frames or file paths
#' @param outcome_data Outcome data frame  
#' @param output_dir Directory to save results
#' @param parallel Use parallel processing (default: FALSE)
#' @param n_cores Number of cores for parallel processing (default: 2)
#' @param ... Additional parameters passed to run_gsmr_analysis
#'
#' @return List of GSMR results for each exposure
#'
#' @export
run_batch_gsmr <- function(exposure_list, 
                          outcome_data, 
                          output_dir,
                          parallel = FALSE,
                          n_cores = 2,
                          ...) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("Starting batch GSMR analysis...\n")
  cat("Number of exposures:", length(exposure_list), "\n")
  
  # Function to process single exposure
  process_exposure <- function(i, exposure_item) {
    
    cat("Processing exposure", i, "of", length(exposure_list), "\n")
    
    # Load exposure data if it's a file path
    if (is.character(exposure_item) && length(exposure_item) == 1) {
      exposure_data <- load_gwas_data(exposure_item, data_type = "exposure")
    } else {
      exposure_data <- exposure_item
    }
    
    # Create exposure-specific output directory
    exposure_output_dir <- file.path(output_dir, paste0("exposure_", i))
    
    # Harmonize data
    harmonised_data <- harmonise_data(exposure_data, outcome_data, action = 2)
    
    if (nrow(harmonised_data) < 3) {
      warning("Not enough SNPs for exposure ", i, ". Skipping...")
      return(NULL)
    }
    
    # Run GSMR analysis
    tryCatch({
      gsmr_result <- run_gsmr_analysis(harmonised_data,
                                      output_dir = exposure_output_dir,
                                      ...)
      return(gsmr_result)
    }, error = function(e) {
      warning("GSMR analysis failed for exposure ", i, ": ", e$message)
      return(NULL)
    })
  }
  
  # Run analysis (parallel or sequential)
  if (parallel && requireNamespace("foreach", quietly = TRUE)) {
    results <- foreach::foreach(i = seq_along(exposure_list)) %do% {
      process_exposure(i, exposure_list[[i]])
    }
  } else {
    results <- lapply(seq_along(exposure_list), function(i) {
      process_exposure(i, exposure_list[[i]])
    })
  }
  
  # Remove NULL results
  results <- results[!sapply(results, is.null)]
  
  # Combine effect estimates
  combined_estimates <- do.call(rbind, lapply(results, function(x) x$effect_estimates))
  combined_estimates$exposure_id <- rep(seq_along(results), 
                                       times = sapply(results, function(x) nrow(x$effect_estimates)))
  
  # Save combined results
  write.table(combined_estimates,
              file.path(output_dir, "combined_gsmr_results.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("Batch GSMR analysis completed!\n")
  cat("Successful analyses:", length(results), "/", length(exposure_list), "\n")
  
  return(results)
}

#' Save GSMR analysis results
#'
#' @description Saves GSMR results to specified directory
#'
#' @param results GSMR results list
#' @param output_dir Output directory path
save_gsmr_results <- function(results, output_dir) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save main GSMR results
  write.table(results$effect_estimates,
              file.path(output_dir, "gsmr_results.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save used SNPs
  write.table(results$used_snps,
              file.path(output_dir, "gsmr_used_snps.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Save outlier SNPs if any
  if (!is.null(results$outlier_snps)) {
    write.table(results$outlier_snps,
                file.path(output_dir, "gsmr_outlier_snps.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  # Save detailed GSMR output
  capture.output(results$gsmr_results, 
                file = file.path(output_dir, "gsmr_detailed_output.txt"))
  
  # Save plot if available
  if (!is.null(results$gsmr_plot)) {
    ggsave(file.path(output_dir, "gsmr_plot.pdf"),
           results$gsmr_plot,
           width = 8, height = 8)
  }
  
  # Save analysis information
  info_text <- paste(
    "=== GSMR Analysis Information ===",
    paste("Input SNPs:", results$analysis_info$n_input_snps),
    paste("SNPs used in analysis:", results$analysis_info$n_used_snps),
    paste("Outlier SNPs:", results$analysis_info$n_outlier_snps),
    paste("GWAS threshold:", results$analysis_info$gwas_threshold),
    paste("HEIDI outlier detection:", results$analysis_info$heidi_enabled),
    sep = "\n"
  )
  
  writeLines(info_text, file.path(output_dir, "gsmr_analysis_info.txt"))
  
  cat("GSMR results saved to:", output_dir, "\n")
}
