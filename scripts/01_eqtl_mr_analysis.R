#' eQTL-MR Analysis Script
#'
#' This script performs two-sample Mendelian randomization analysis using
#' eQTLGen expression QTL data and MDD GWAS data to identify causal genes
#' for Major Depressive Disorder. 
#'
#' @title Bulk Tissue eQTL-MR Analysis for Major Depressive Disorder
#' @description 
#' Conducts comprehensive MR analysis including:
#' - Data loading and harmonization
#' - Multiple MR methods (IVW, MR-Egger, Weighted Median, GSMR)
#' - Sensitivity analyses (pleiotropy, heterogeneity, MR-PRESSO)
#' - Visualization (scatter plots, forest plots, funnel plots)
#'
#' @author Xuejiao Hou
#' @date 2025-02-10
#' @version 1.0.0

# ==============================================================================
# SETUP AND CONFIGURATION
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(gsmr2)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(foreach)
  library(progress)
  library(logging)
})

# Source custom functions
source(here::here("R", "mr_analysis.R"))
source(here::here("R", "gsmr_analysis.R"))
source(here::here("R", "visualization.R"))

# Setup logging
basicConfig(level = "INFO")
log_file <- file.path("logs", paste0("eqtl_mr_analysis_", Sys.Date(), ".log"))
if (!dir.exists("logs")) dir.create("logs")
addHandler(writeToFile, file = log_file)

# Configuration parameters
config <- list(
  # Input data paths (modify these to match your data location)
  eqtl_data_path = "data/example_eqtl_data.txt",  # Your eQTL data path
  mdd_gwas_path = "data/example_mdd_gwas.txt",    # Your MDD GWAS path
  
  # Analysis parameters
  pval_threshold = 5e-8,          # P-value threshold for instrument selection
  clump_r2 = 0.1,              # LD clumping R² threshold  
  clump_kb = 10000,              # LD clumping distance in kb
  f_stat_threshold = 10,          # F-statistic threshold for weak instruments
  
  # Output directories
  output_base_dir = "results/eqtl_mr",
  figures_dir = "results/figures/eqtl_mr",
  
  # Analysis options
  perform_gsmr = TRUE,
  perform_sensitivity = TRUE,
  create_plots = TRUE,
  save_detailed_results = TRUE
)

# Create output directories
output_dirs <- c(config$output_base_dir, config$figures_dir, 
                paste0(config$output_base_dir, "/significant_results"),
                paste0(config$output_base_dir, "/all_results"))

for (dir in output_dirs) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    loginfo(paste("Created directory:", dir))
  }
}

# ==============================================================================
# DATA LOADING AND PREPARATION  
# ==============================================================================

loginfo("Starting eQTL-MR analysis for Major Depressive Disorder")
loginfo(paste("Configuration: p-value threshold =", config$pval_threshold))
loginfo(paste("F-statistic threshold =", config$f_stat_threshold))

# Load MDD GWAS outcome data
loginfo("Loading MDD GWAS outcome data...")
mdd_data <- load_gwas_data(
  file_path = config$mdd_gwas_path,
  data_type = "outcome",
  snp_col = "ID",              # Adjust these column names to match your data
  beta_col = "BETA",
  se_col = "SE", 
  pval_col = "PVAL",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  eaf_col = "eaf",
  samplesize_col = "NEFF"
)

# Validate MDD GWAS data
check_data_format(mdd_data, "outcome")
loginfo(paste("Loaded", nrow(mdd_data), "SNPs from MDD GWAS"))

# Load eQTL exposure data (this could be a single file or batch processing)
loginfo("Loading eQTL exposure data...")
eqtl_data <- load_gwas_data(
  file_path = config$eqtl_data_path,
  data_type = "exposure",
  snp_col = "SNP",             # Adjust these column names to match your data
  beta_col = "beta",
  se_col = "se",
  pval_col = "pval", 
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  eaf_col = "eaf",
  samplesize_col = "samplesize"
)

# Add exposure identifier (gene name)
if (!"exposure" %in% names(eqtl_data)) {
  # Extract gene name from filename or add manually
  gene_name <- tools::file_path_sans_ext(basename(config$eqtl_data_path))
  eqtl_data$exposure <- gene_name
  loginfo(paste("Set exposure name to:", gene_name))
}

# ==============================================================================
# SINGLE GENE MR ANALYSIS EXAMPLE
# ==============================================================================

# Example: Analyze CKAP2 gene (modify this section for your specific gene)
gene_name <- "CKAP2"  # Replace with your gene of interest
loginfo(paste("Performing MR analysis for gene:", gene_name))

# Filter eQTL data for the specific gene (if analyzing multiple genes)
if ("gene" %in% names(eqtl_data)) {
  gene_eqtl_data <- eqtl_data[eqtl_data$gene == gene_name, ]
} else {
  gene_eqtl_data <- eqtl_data
  gene_eqtl_data$exposure <- gene_name
}

loginfo(paste("Found", nrow(gene_eqtl_data), "SNPs for", gene_name))

# Perform comprehensive MR analysis
mr_results <- perform_mr_analysis(
  exposure_data = gene_eqtl_data,
  outcome_data = mdd_data,
  clump_r2 = config$clump_r2,
  clump_kb = config$clump_kb,
  pval_threshold = config$pval_threshold,
  f_stat_threshold = config$f_stat_threshold,
  output_dir = file.path(config$output_base_dir, gene_name),
  create_plots = config$create_plots
)

loginfo(paste("MR analysis completed for", gene_name))

# Extract key results for summary
ivw_result <- mr_results$mr_results[mr_results$mr_results$method == "Inverse variance weighted", ]
loginfo(paste("IVW result - OR:", round(ivw_result$or, 4), 
              "95% CI: [", round(ivw_result$or_lci95, 4), ",", 
              round(ivw_result$or_uci95, 4), "]",
              "P-value:", format(ivw_result$pval, scientific = TRUE)))

# GSMR analysis if requested
if (config$perform_gsmr) {
  loginfo("Performing GSMR analysis...")
  
  gsmr_results <- run_gsmr_analysis(
    harmonised_data = mr_results$harmonised_data,
    output_dir = file.path(config$output_base_dir, gene_name, "gsmr")
  )
  
  loginfo(paste("GSMR result - Beta:", round(gsmr_results$effect_estimates$beta, 4),
                "P-value:", format(gsmr_results$effect_estimates$pval, scientific = TRUE)))
}

# ==============================================================================
# BATCH ANALYSIS FOR MULTIPLE GENES (Optional)
# ==============================================================================

# If you have multiple genes/eQTL files to analyze
if (FALSE) {  # Set to TRUE to enable batch processing
  
  loginfo("Starting batch eQTL-MR analysis...")
  
  # List of eQTL files (modify to match your data structure)
  eqtl_files <- list.files("data/eqtl_data/", pattern = "\\.txt$", full.names = TRUE)
  
  # Initialize results storage
  batch_results <- list()
  n_files <- length(eqtl_files)
  
  # Progress bar
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent in :elapsed, ETA: :eta",
    total = n_files, clear = FALSE, width = 60
  )
  
  # Process each gene
  for (i in 1:n_files) {
    file_path <- eqtl_files[i]
    gene_name <- tools::file_path_sans_ext(basename(file_path))
    
    pb$tick()
    loginfo(paste("Processing gene", i, "of", n_files, ":", gene_name))
    
    tryCatch({
      # Load gene-specific eQTL data
      gene_eqtl <- load_gwas_data(file_path, data_type = "exposure")
      gene_eqtl$exposure <- gene_name
      
      # Skip if insufficient SNPs
      if (nrow(gene_eqtl) < 10) {
        logwarn(paste("Skipping", gene_name, "- insufficient SNPs"))
        next
      }
      
      # Run MR analysis
      gene_mr_results <- perform_mr_analysis(
        exposure_data = gene_eqtl,
        outcome_data = mdd_data,
        pval_threshold = config$pval_threshold,
        f_stat_threshold = config$f_stat_threshold,
        output_dir = file.path(config$output_base_dir, "all_results", gene_name),
        create_plots = FALSE  # Save time in batch processing
      )
      
      # Store results
      batch_results[[i]] <- list(
        gene = gene_name,
        mr_results = gene_mr_results$mr_results,
        n_snps = nrow(gene_mr_results$harmonised_data),
        instrument_stats = gene_mr_results$instrument_stats
      )
      
      # Save significant results separately
      ivw_result <- gene_mr_results$mr_results[
        gene_mr_results$mr_results$method == "Inverse variance weighted", ]
      
      if (ivw_result$pval < 0.05) {
        loginfo(paste("Significant result for", gene_name, 
                     "- saving detailed results"))
        
        # Save detailed results for significant findings
        gene_output_dir <- file.path(config$output_base_dir, "significant_results", gene_name)
        
        # Re-run with full diagnostics for significant results
        detailed_results <- perform_mr_analysis(
          exposure_data = gene_eqtl,
          outcome_data = mdd_data,
          output_dir = gene_output_dir,
          create_plots = TRUE
        )
        
        # GSMR for significant results
        if (config$perform_gsmr) {
          run_gsmr_analysis(
            harmonised_data = detailed_results$harmonised_data,
            output_dir = file.path(gene_output_dir, "gsmr")
          )
        }
      }
      
    }, error = function(e) {
      logerror(paste("Error processing", gene_name, ":", e$message))
    })
  }
  
  # Combine batch results
  loginfo("Combining batch results...")
  combined_results <- do.call(rbind, lapply(batch_results, function(x) {
    if (!is.null(x)) {
      ivw_row <- x$mr_results[x$mr_results$method == "Inverse variance weighted", ]
      data.frame(
        gene = x$gene,
        nsnp = x$n_snps,
        beta = ivw_row$b,
        se = ivw_row$se,
        pval = ivw_row$pval,
        or = ivw_row$or,
        or_lci95 = ivw_row$or_lci95,
        or_uci95 = ivw_row$or_uci95,
        total_snps = x$instrument_stats$n_total_snps,
        strong_instruments = x$instrument_stats$n_strong_instruments,
        mean_f_stat = x$instrument_stats$mean_f_statistic,
        stringsAsFactors = FALSE
      )
    }
  }))
  
  # Apply multiple testing correction
  combined_results$fdr <- p.adjust(combined_results$pval, method = "fdr")
  combined_results$bonferroni <- p.adjust(combined_results$pval, method = "bonferroni")
  
  # Sort by p-value
  combined_results <- combined_results[order(combined_results$pval), ]
  
  # Save combined results
  write.table(combined_results,
              file.path(config$output_base_dir, "eqtl_mdd_combined_results.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  loginfo(paste("Batch analysis completed.", nrow(combined_results), "genes analyzed"))
  loginfo(paste("Significant results (p < 0.05):", 
                sum(combined_results$pval < 0.05, na.rm = TRUE)))
  loginfo(paste("FDR significant (q < 0.05):",
                sum(combined_results$fdr < 0.05, na.rm = TRUE)))
}

# ==============================================================================
# VISUALIZATION AND SUMMARY
# ==============================================================================

# Create summary visualizations if results exist
if (exists("combined_results") && nrow(combined_results) > 0) {
  
  loginfo("Creating summary visualizations...")
  
  # Volcano plot
  volcano_plot <- create_volcano_plot(
    mr_results = combined_results,
    pval_threshold = 0.05,
    bonferroni_correction = TRUE,
    title = "eQTL-MR Analysis: Gene Expression → MDD Risk"
  )
  
  save_publication_plot(
    volcano_plot,
    file.path(config$figures_dir, "eqtl_mdd_volcano_plot"),
    format = "pdf", width = 10, height = 8
  )
  
  # Forest plot for top results
  top_results <- head(combined_results[combined_results$pval < 0.05, ], 20)
  
  if (nrow(top_results) > 0) {
    # Add database column for forest plot
    top_results$database <- "eQTLGen"
    
    forest_plot <- create_forest_plot(
      mr_results = top_results,
      group_var = "database",
      title = "Forest Plot: Top eQTL-MR Results for MDD"
    )
    
    save_publication_plot(
      forest_plot,
      file.path(config$figures_dir, "eqtl_mdd_forest_plot"),
      format = "pdf", width = 12, height = 8
    )
  }
}

# ==============================================================================
# GENERATE ANALYSIS REPORT
# ==============================================================================

# Create analysis summary report
report_file <- file.path(config$output_base_dir, "analysis_summary_report.txt")
report_content <- c(
  "=================================================================",
  "eQTL-MR Analysis Summary Report",
  "=================================================================",
  paste("Analysis Date:", Sys.Date()),
  paste("Analysis Time:", Sys.time()),
  "",
  "CONFIGURATION:",
  paste("  P-value threshold:", config$pval_threshold),
  paste("  F-statistic threshold:", config$f_stat_threshold),
  paste("  LD clumping R²:", config$clump_r2),
  paste("  LD clumping distance:", config$clump_kb, "kb"),
  "",
  "INPUT DATA:",
  paste("  eQTL data path:", config$eqtl_data_path),
  paste("  MDD GWAS path:", config$mdd_gwas_path),
  ""
)

# Add single gene results if available
if (exists("mr_results")) {
  ivw_res <- mr_results$mr_results[mr_results$mr_results$method == "Inverse variance weighted", ]
  report_content <- c(report_content,
    "SINGLE GENE ANALYSIS RESULTS:",
    paste("  Gene analyzed:", gene_name),
    paste("  Number of instruments:", nrow(mr_results$harmonised_data)),
    paste("  IVW Odds Ratio:", round(ivw_res$or, 4)),
    paste("  IVW 95% CI: [", round(ivw_res$or_lci95, 4), ",", round(ivw_res$or_uci95, 4), "]"),
    paste("  IVW P-value:", format(ivw_res$pval, scientific = TRUE)),
    ""
  )
}

# Add batch results summary if available
if (exists("combined_results")) {
  n_total <- nrow(combined_results)  
  n_significant <- sum(combined_results$pval < 0.05, na.rm = TRUE)
  n_fdr_sig <- sum(combined_results$fdr < 0.05, na.rm = TRUE)
  
  report_content <- c(report_content,
    "BATCH ANALYSIS RESULTS:",
    paste("  Total genes analyzed:", n_total),
    paste("  Nominally significant (p < 0.05):", n_significant),
    paste("  FDR significant (q < 0.05):", n_fdr_sig),
    paste("  Proportion significant:", round(n_significant/n_total * 100, 2), "%"),
    ""
  )
  
  # Top results
  if (n_significant > 0) {
    top_genes <- head(combined_results$gene[combined_results$pval < 0.05], 10)
    report_content <- c(report_content,
      "TOP SIGNIFICANT GENES:",
      paste("  ", top_genes),
      ""
    )
  }
}

report_content <- c(report_content,
  "OUTPUT LOCATIONS:",
  paste("  Main results:", config$output_base_dir),
  paste("  Figures:", config$figures_dir),
  paste("  Log file:", log_file),
  "",
  "=================================================================",
  "Analysis completed successfully!",
  "================================================================="
)

writeLines(report_content, report_file)
loginfo("Analysis summary report saved to: " + report_file)

# ==============================================================================
# CLEANUP AND SESSION INFO
# ==============================================================================

loginfo("eQTL-MR analysis pipeline completed successfully!")
loginfo(paste("Results saved to:", config$output_base_dir))
loginfo(paste("Figures saved to:", config$figures_dir))

# Save session information for reproducibility
session_info_file <- file.path(config$output_base_dir, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_info_file)

# Print final summary to console
cat("\n=================================================================\n")
cat("eQTL-MR Analysis Completed Successfully!\n")
cat("=================================================================\n")
if (exists("combined_results")) {
  cat("Genes analyzed:", nrow(combined_results), "\n")
  cat("Significant results:", sum(combined_results$pval < 0.05, na.rm = TRUE), "\n")
}
cat("Results directory:", config$output_base_dir, "\n")
cat("Log file:", log_file, "\n")
cat("=================================================================\n")

# End logging
removeHandler("basic.stdout")
