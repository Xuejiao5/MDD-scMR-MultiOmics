#' Visualization Functions for Multi-omics Mendelian Randomization
#'
#' This module contains functions for creating publication-quality plots including
#' volcano plots, forest plots, Venn diagrams, and composite figures.
#'
#' @author Xuejiao Hou
#' @export

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(patchwork)
  library(VennDiagram)
  library(RColorBrewer)
  library(magick)
  library(data.table)
})

#' Create volcano plot for MR analysis results
#'
#' @description Creates a volcano plot showing effect sizes vs -log10(p-values)
#'              with customizable significance thresholds and gene labeling
#'
#' @param mr_results Data frame with MR results containing columns: id, or, pval
#' @param pval_threshold P-value threshold for significance line (default: 0.05)
#' @param bonferroni_correction Apply Bonferroni correction (default: TRUE)
#' @param n_tests Number of tests for Bonferroni correction (auto-detected if NULL)
#' @param label_genes Vector of gene names to label (default: NULL for auto-labeling)
#' @param label_top_n Number of top significant genes to label (default: 10)
#' @param colors Named vector of colors for different point types
#' @param title Plot title (default: "Volcano Plot: MR Analysis Results")
#' @param x_limits X-axis limits (default: NULL for auto)
#' @param y_limits Y-axis limits (default: NULL for auto)
#'
#' @return ggplot2 object
#'
#' @examples
#' \dontrun{
#' # Create volcano plot
#' volcano_plot <- create_volcano_plot(
#'   mr_results = eqtl_mdd_results,
#'   pval_threshold = 0.05,
#'   label_top_n = 5
#' )
#' print(volcano_plot)
#' }
#'
#' @export
create_volcano_plot <- function(mr_results,
                                pval_threshold = 0.05,
                                bonferroni_correction = TRUE,
                                n_tests = NULL,
                                label_genes = NULL,
                                label_top_n = 10,
                                colors = c("grey" = "grey60", "protective" = "#14365F", "risk" = "#D64F38"),
                                title = "Volcano Plot: MR Analysis Results",
                                x_limits = NULL,
                                y_limits = NULL) {
  
  # Prepare data
  plot_data <- mr_results %>%
    mutate(
      log_or = log(or),
      neg_log_p = -log10(pval),
      gene_name = id
    ) %>%
    filter(!is.na(or) & !is.na(pval) & is.finite(log_or) & is.finite(neg_log_p))
  
  # Calculate significance threshold
  if (is.null(n_tests)) {
    n_tests <- nrow(plot_data)
  }
  
  if (bonferroni_correction) {
    sig_threshold <- -log10(pval_threshold / n_tests)
    threshold_label <- paste0("Bonferroni (p=", pval_threshold, "/", n_tests, ")")
  } else {
    sig_threshold <- -log10(pval_threshold)
    threshold_label <- paste0("p=", pval_threshold)
  }
  
  # Assign colors based on significance and direction
  plot_data <- plot_data %>%
    mutate(
      point_color = case_when(
        neg_log_p > sig_threshold & or > 1 ~ "risk",
        neg_log_p > sig_threshold & or < 1 ~ "protective", 
        TRUE ~ "grey"
      )
    )
  
  # Determine genes to label
  if (is.null(label_genes)) {
    # Auto-select top significant genes
    significant_genes <- plot_data %>%
      filter(neg_log_p > sig_threshold) %>%
      arrange(pval) %>%
      slice_head(n = label_top_n) %>%
      pull(gene_name)
    
    plot_data$label <- ifelse(plot_data$gene_name %in% significant_genes, 
                             plot_data$gene_name, NA)
  } else {
    plot_data$label <- ifelse(plot_data$gene_name %in% label_genes,
                             plot_data$gene_name, NA)
  }
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = log_or, y = neg_log_p)) +
    geom_point(aes(color = point_color, size = neg_log_p), alpha = 0.7) +
    scale_color_manual(values = colors, guide = "none") +
    scale_size_continuous(range = c(1, 4), guide = "none") +
    
    # Add significance threshold line
    geom_hline(yintercept = sig_threshold, 
               linetype = "dashed", color = "#76A2B9", size = 0.8) +
    
    # Add OR = 1 reference line
    geom_vline(xintercept = 0, 
               linetype = "dashed", color = "#76A2B9", size = 0.8) +
    
    # Add nominal significance line if using Bonferroni
    {if(bonferroni_correction) 
      geom_hline(yintercept = -log10(pval_threshold), 
                 linetype = "dashed", color = "grey50", size = 0.5)},
    
    # Add gene labels
    geom_label_repel(
      aes(label = label),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.size = 0.3,
      segment.alpha = 0.6,
      max.overlaps = 20,
      na.rm = TRUE
    ) +
    
    # Customize axes
    labs(
      x = "ln(Odds Ratio)",
      y = "-log₁₀(P-value)", 
      title = title,
      caption = paste0("Significance threshold: ", threshold_label, 
                      "\nProtective effects (OR < 1) in blue, Risk effects (OR > 1) in red")
    ) +
    
    # Set axis limits if provided
    {if(!is.null(x_limits)) xlim(x_limits)} +
    {if(!is.null(y_limits)) ylim(y_limits)} +
    
    # Theme
    theme_classic() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.caption = element_text(hjust = 0, size = 9),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.border = element_rect(fill = NA, color = "black", size = 0.5)
    )
  
  return(p)
}

#' Create forest plot for MR results
#'
#' @description Creates a publication-quality forest plot showing odds ratios
#'              with confidence intervals, organized by groups
#'
#' @param mr_results Data frame with MR results
#' @param group_var Column name for grouping (default: "database") 
#' @param sort_by Column to sort results by (default: "pval")
#' @param show_heterogeneity Include heterogeneity test results (default: TRUE)
#' @param colors Color scheme for different groups
#' @param title Plot title
#'
#' @return ggplot2 object
#'
#' @export
create_forest_plot <- function(mr_results,
                               group_var = "database",
                               sort_by = "pval",
                               show_heterogeneity = TRUE,
                               colors = c("#72B6A1", "#4189C8", "#F08E64"),
                               title = "Forest Plot: MR Analysis Results") {
  
  # Prepare data
  plot_data <- mr_results %>%
    arrange(!!sym(sort_by)) %>%
    mutate(
      row_id = rev(row_number()),
      or_ci = sprintf("%.3f (%.3f, %.3f)", or, or_lci95, or_uci95),
      p_formatted = case_when(
        pval < 0.001 ~ sprintf("%.2e", pval),
        TRUE ~ sprintf("%.3f", pval)
      ),
      significant = pval < 0.05
    )
  
  # Create alternating background
  if (!is.null(group_var) && group_var %in% names(plot_data)) {
    plot_data <- plot_data %>%
      mutate(group_numeric = as.numeric(factor(!!sym(group_var)))) %>%
      mutate(alternating_bg = group_numeric %% 2 == 0)
  } else {
    plot_data <- plot_data %>%
      mutate(alternating_bg = row_id %% 2 == 0)
  }
  
  # Position parameters
  gene_x <- -1.2
  group_x <- -0.6
  nsnp_x <- -0.1
  or_ci_x <- 1.8
  pval_x <- 2.6
  plot_max_x <- 3.2
  
  # Create the plot
  p <- ggplot(plot_data, aes(y = factor(row_id, levels = row_id))) +
    
    # Add alternating background
    geom_rect(aes(xmin = -1.4, xmax = plot_max_x,
                  ymin = row_id - 0.5, ymax = row_id + 0.5,
                  fill = alternating_bg),
              alpha = 0.3) +
    scale_fill_manual(values = c("TRUE" = "#E8E8E8", "FALSE" = "white"), guide = "none") +
    
    # Add reference line at OR = 1
    geom_vline(xintercept = 1, linetype = "dashed", color = "#76A2B9", size = 0.8) +
    
    # Add confidence intervals
    geom_errorbarh(aes(xmin = or_lci95, xmax = or_uci95, 
                       color = factor(significant)),
                   height = 0.3, size = 1.2) +
    
    # Add point estimates
    geom_point(aes(x = or, color = factor(significant)), size = 3) +
    
    # Color scheme
    scale_color_manual(values = c("TRUE" = colors[2], "FALSE" = colors[1]), guide = "none") +
    
    # Add text labels
    geom_text(aes(x = gene_x, label = id), hjust = 0, size = 3.5) +
    {if (!is.null(group_var) && group_var %in% names(plot_data))
      geom_text(aes(x = group_x, label = !!sym(group_var)), hjust = 0, size = 3.5)} +
    geom_text(aes(x = nsnp_x, label = nsnp), hjust = 0.5, size = 3.5) +
    geom_text(aes(x = or_ci_x, label = or_ci), hjust = 0, size = 3.5) +
    geom_text(aes(x = pval_x, label = p_formatted), hjust = 0, size = 3.5) +
    
    # Add column headers
    annotate("text", x = gene_x, y = max(plot_data$row_id) + 1,
             label = "Gene", size = 4, fontface = "bold", hjust = 0) +
    {if (!is.null(group_var) && group_var %in% names(plot_data))
      annotate("text", x = group_x, y = max(plot_data$row_id) + 1,
               label = str_to_title(group_var), size = 4, fontface = "bold", hjust = 0)} +
    annotate("text", x = nsnp_x, y = max(plot_data$row_id) + 1,
             label = "N SNPs", size = 4, fontface = "bold", hjust = 0.5) +
    annotate("text", x = or_ci_x, y = max(plot_data$row_id) + 1,
             label = "OR (95% CI)", size = 4, fontface = "bold", hjust = 0) +
    annotate("text", x = pval_x, y = max(plot_data$row_id) + 1,
             label = "P-value", size = 4, fontface = "bold", hjust = 0) +
    
    # Add separator lines
    geom_segment(x = -1.4, xend = plot_max_x,
                y = max(plot_data$row_id) + 0.5, yend = max(plot_data$row_id) + 0.5,
                color = "black", size = 0.5) +
    geom_segment(x = -1.4, xend = plot_max_x,
                y = 0.5, yend = 0.5,
                color = "black", size = 0.5) +
    
    # Customize axes and theme
    scale_x_continuous(limits = c(-1.4, plot_max_x),
                      breaks = c(0, 0.5, 1, 1.5, 2),
                      labels = c("0", "0.5", "1", "1.5", "2")) +
    scale_y_discrete(expand = expansion(add = 1)) +
    labs(x = "Odds Ratio", y = "", title = title) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title.x = element_text(size = 12),
      axis.text.x = element_text(size = 10),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  return(p)
}

#' Create Venn diagram for multi-dataset comparison
#'
#' @description Creates publication-quality Venn diagrams to show overlap
#'              between different datasets or analyses
#'
#' @param gene_lists Named list of character vectors (gene sets)
#' @param colors Vector of colors for each set
#' @param output_file Output file path (with extension: .png, .pdf, .tiff)
#' @param width Image width in inches (default: 8)
#' @param height Image height in inches (default: 8) 
#' @param resolution Image resolution in DPI (default: 300)
#' @param title Plot title
#'
#' @return File path of saved plot
#'
#' @examples
#' \dontrun{
#' # Create Venn diagram
#' gene_lists <- list(
#'   eQTL = c("GENE1", "GENE2", "GENE3"),
#'   pQTL = c("GENE2", "GENE3", "GENE4"),
#'   scRNA = c("GENE3", "GENE4", "GENE5")
#' )
#' 
#' venn_plot <- create_venn_diagram(
#'   gene_lists = gene_lists,
#'   output_file = "gene_overlap.png"
#' )
#' }
#'
#' @export
create_venn_diagram <- function(gene_lists,
                                colors = c("#3182bd", "#e6550d", "#31a354", "#756bb1"),
                                output_file = "venn_diagram.png",
                                width = 8,
                                height = 8,
                                resolution = 300,
                                title = NULL) {
  
  # Validate input
  if (!is.list(gene_lists) || is.null(names(gene_lists))) {
    stop("gene_lists must be a named list")
  }
  
  n_sets <- length(gene_lists)
  if (n_sets < 2 || n_sets > 4) {
    stop("Venn diagrams support 2-4 sets only")
  }
  
  # Create temporary file if no extension provided
  if (!grepl("\\.(png|pdf|tiff)$", output_file)) {
    output_file <- paste0(output_file, ".png")
  }
  
  # Prepare colors
  if (length(colors) < n_sets) {
    colors <- rainbow(n_sets)
  }
  colors <- colors[1:n_sets]
  
  # Create Venn diagram based on number of sets
  if (n_sets == 2) {
    venn_plot <- venn.diagram(
      x = gene_lists,
      filename = output_file,
      imagetype = tools::file_ext(output_file),
      height = height * 100,
      width = width * 100,
      resolution = resolution,
      compression = "lzw",
      units = "px",
      lwd = 2,
      lty = "solid",
      col = colors,
      fill = alpha(colors, 0.4),
      cex = 1.5,
      fontface = "bold",
      cat.cex = 1.8,
      cat.fontface = "bold",
      cat.dist = c(0.05, 0.05),
      cat.col = colors,
      margin = 0.1,
      main = title,
      main.cex = 2
    )
  } else if (n_sets == 3) {
    venn_plot <- venn.diagram(
      x = gene_lists,
      filename = output_file,
      imagetype = tools::file_ext(output_file),
      height = height * 100,
      width = width * 100,
      resolution = resolution,
      compression = "lzw",
      units = "px",
      lwd = 2,
      lty = "solid",
      col = colors,
      fill = alpha(colors, 0.4),
      cex = 1.5,
      fontface = "bold",
      cat.cex = 1.8,
      cat.fontface = "bold",
      cat.dist = c(0.08, 0.08, 0.08),
      cat.col = colors,
      margin = 0.1,
      main = title,
      main.cex = 2,
      rotation.degree = 0
    )
  } else if (n_sets == 4) {
    venn_plot <- venn.diagram(
      x = gene_lists,
      filename = output_file,
      imagetype = tools::file_ext(output_file),
      height = height * 100,
      width = width * 100,
      resolution = resolution,
      compression = "lzw",
      units = "px",
      lwd = 2,
      lty = "solid",
      col = colors,
      fill = alpha(colors, 0.35),
      cex = 1.2,
      fontface = "bold",
      cat.cex = 1.5,
      cat.fontface = "bold",
      cat.dist = c(0.22, 0.22, 0.22, 0.22),
      cat.col = colors,
      margin = 0.1,
      main = title,
      main.cex = 2
    )
  }
  
  cat("Venn diagram saved to:", output_file, "\n")
  return(output_file)
}

#' Create composite figure with multiple panels
#'
#' @description Combines multiple ggplot objects into a publication-ready figure
#'
#' @param plot_list Named list of ggplot objects
#' @param layout Layout specification (e.g., "AB/CD" for 2x2 grid)
#' @param labels Panel labels (default: uppercase letters)
#' @param label_size Size of panel labels (default: 16)
#' @param title Overall figure title (optional)
#' @param output_file Output file path (optional)
#' @param width Figure width in inches (default: 12)
#' @param height Figure height in inches (default: 10)
#' @param dpi Resolution in DPI (default: 300)
#'
#' @return patchwork object
#'
#' @examples
#' \dontrun{
#' # Create composite figure
#' plot_list <- list(
#'   volcano = volcano_plot,
#'   forest = forest_plot
#' )
#' 
#' composite <- create_composite_figure(
#'   plot_list = plot_list,
#'   layout = "AB",
#'   output_file = "figure1.pdf"
#' )
#' }
#'
#' @export
create_composite_figure <- function(plot_list,
                                   layout = NULL,
                                   labels = NULL,
                                   label_size = 16,
                                   title = NULL,
                                   output_file = NULL,
                                   width = 12,
                                   height = 10,
                                   dpi = 300) {
  
  # Generate default layout if not provided
  if (is.null(layout)) {
    n_plots <- length(plot_list)
    if (n_plots <= 2) {
      layout <- paste(toupper(letters[1:n_plots]), collapse = "")
    } else if (n_plots <= 4) {
      layout <- paste0(paste(toupper(letters[1:2]), collapse = ""), "/",
                      paste(toupper(letters[3:min(4, n_plots)]), collapse = ""))
    } else {
      # Default to grid layout
      n_cols <- ceiling(sqrt(n_plots))
      layout <- NULL  # Let patchwork handle it
    }
  }
  
  # Generate default labels if not provided
  if (is.null(labels)) {
    labels <- toupper(letters[1:length(plot_list)])
  }
  
  # Create composite plot
  if (!is.null(layout)) {
    composite_plot <- wrap_plots(plot_list, design = layout) +
      plot_annotation(
        tag_levels = list(labels),
        tag_prefix = "",
        tag_suffix = "",
        theme = theme(plot.tag = element_text(size = label_size, face = "bold")),
        title = title
      )
  } else {
    composite_plot <- wrap_plots(plot_list) +
      plot_annotation(
        tag_levels = list(labels),
        tag_prefix = "",
        tag_suffix = "",
        theme = theme(plot.tag = element_text(size = label_size, face = "bold")),
        title = title
      )
  }
  
  # Save if output file specified
  if (!is.null(output_file)) {
    ggsave(output_file, composite_plot,
           width = width, height = height, dpi = dpi,
           bg = "white")
    cat("Composite figure saved to:", output_file, "\n")
  }
  
  return(composite_plot)
}

#' Create advanced composite figure using magick
#'
#' @description Creates complex multi-panel figures with precise positioning
#'
#' @param image_files Vector of image file paths
#' @param layout_matrix Matrix specifying image positions
#' @param canvas_width Canvas width in pixels (default: 3000)
#' @param canvas_height Canvas height in pixels (default: 4000)
#' @param labels Panel labels to add
#' @param label_size Label font size (default: 120)
#' @param label_positions List of label positions (x, y coordinates)
#' @param output_file Output file path
#' @param quality Image quality 1-100 (default: 95)
#'
#' @return magick image object
#'
#' @export
create_advanced_composite <- function(image_files,
                                     layout_matrix = NULL,
                                     canvas_width = 3000,
                                     canvas_height = 4000,
                                     labels = NULL,
                                     label_size = 120,
                                     label_positions = NULL,
                                     output_file = NULL,
                                     quality = 95) {
  
  # Read images
  images <- lapply(image_files, image_read)
  
  # Create white canvas
  canvas <- image_blank(canvas_width, canvas_height, color = "white")
  
  # Add images based on layout matrix or positions
  if (!is.null(layout_matrix)) {
    # Use matrix layout
    n_rows <- nrow(layout_matrix)
    n_cols <- ncol(layout_matrix)
    
    img_width <- floor(canvas_width / n_cols)
    img_height <- floor(canvas_height / n_rows)
    
    for (i in 1:length(images)) {
      # Find position in matrix
      pos <- which(layout_matrix == i, arr.ind = TRUE)
      if (nrow(pos) > 0) {
        x_pos <- (pos[1, "col"] - 1) * img_width
        y_pos <- (pos[1, "row"] - 1) * img_height
        
        # Resize and composite image
        img_resized <- image_scale(images[[i]], paste0(img_width, "x", img_height))
        canvas <- image_composite(canvas, img_resized, 
                                offset = paste0("+", x_pos, "+", y_pos))
      }
    }
  } else {
    # Manual positioning (implement as needed)
    warning("Manual positioning not implemented yet. Using default layout.")
  }
  
  # Add labels if specified
  if (!is.null(labels)) {
    if (is.null(label_positions)) {
      # Generate default positions
      label_positions <- list()
      for (i in 1:length(labels)) {
        label_positions[[i]] <- list(x = 100, y = 100 + (i-1) * 200)
      }
    }
    
    for (i in 1:length(labels)) {
      if (i <= length(label_positions)) {
        pos <- label_positions[[i]]
        canvas <- image_annotate(
          canvas, labels[i],
          size = label_size,
          font = "Arial Bold",
          color = "black",
          location = paste0("+", pos$x, "+", pos$y)
        )
      }
    }
  }
  
  # Save if output file specified
  if (!is.null(output_file)) {
    image_write(canvas, output_file, quality = quality)
    cat("Advanced composite saved to:", output_file, "\n")
  }
  
  return(canvas)
}

#' Save plot with publication specifications
#'
#' @description Standardized function to save plots with publication settings
#'
#' @param plot ggplot object
#' @param filename Output filename
#' @param format File format ("pdf", "png", "tiff", "eps")
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution for raster formats
#' @param compression Compression for TIFF ("lzw", "none")
#'
#' @return TRUE if successful
#'
#' @export
save_publication_plot <- function(plot,
                                 filename,
                                 format = "pdf",
                                 width = 8,
                                 height = 8,
                                 dpi = 300,
                                 compression = "lzw") {
  
  # Determine file extension
  if (!grepl(paste0("\\.", format, "$"), filename)) {
    filename <- paste0(filename, ".", format)
  }
  
  # Save based on format
  tryCatch({
    if (format == "tiff") {
      ggsave(filename, plot, 
             width = width, height = height, dpi = dpi,
             compression = compression, bg = "white")
    } else {
      ggsave(filename, plot,
             width = width, height = height, dpi = dpi,
             bg = "white")
    }
    
    cat("Plot saved to:", filename, "\n")
    return(TRUE)
    
  }, error = function(e) {
    warning("Failed to save plot: ", e$message)
    return(FALSE)
  })
}
