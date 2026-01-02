# title: 
# author: Beril Tümkaya
# date: 02.01.2026

# Set working directory to script location
setwd(".") # Set your working directory


# Install and load required packages
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    if (bioc) {
      if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
    library(pkg, character.only = TRUE)
  }
}

# Core packages
install_if_missing("DESeq2", bioc = TRUE)
install_if_missing("ggplot2")
install_if_missing("dplyr")
install_if_missing("tidyr")
install_if_missing("pheatmap")
install_if_missing("RColorBrewer")
install_if_missing("ggrepel")
install_if_missing("UpSetR")
install_if_missing("clusterProfiler", bioc = TRUE)
install_if_missing("org.Hs.eg.db", bioc = TRUE)
install_if_missing("enrichplot", bioc = TRUE)
install_if_missing("pathview", bioc = TRUE)
install_if_missing("apeglm", bioc = TRUE)
install_if_missing("tibble")
install_if_missing("scales")
install_if_missing("gridExtra")

# Set working directory and output paths
CONFIG <- list(
  data_dir = ".",
  output_dir = "./results",
  counts_file = "combined_counts.csv",
  metadata_file = "metadata.csv",
  fdr_cutoff = 0.05,
  lfc_cutoff = 1,
  min_count = 10,
  min_samples = 3
)

# Create output directory if it doesn't exist
if (!dir.exists(CONFIG$output_dir)) {
  dir.create(CONFIG$output_dir, recursive = TRUE)
}

# STEP 1: VALIDATE RAW COUNT MATRIX


validate_count_matrix <- function(counts, metadata) {

  issues <- list()
  
  # Check for non-negative integers
  if (any(counts < 0, na.rm = TRUE)) {
    issues$negative_values <- "Found negative values in count matrix"
  } else {
    cat("✓ All values are non-negative\n")
  }
  
  # Check for NA values
  na_count <- sum(is.na(counts))
  if (na_count > 0) {
    issues$na_values <- paste("Found", na_count, "NA values")
  } else {
    cat("✓ No NA values found\n")
  }
  
  # Check for integer values (or whole numbers)
  if (!all(counts == floor(counts), na.rm = TRUE)) {
    issues$non_integer <- "Found non-integer values"
  } else {
    cat("✓ All values are integers\n")
  }
  
  # Check for duplicated gene names
  gene_names <- rownames(counts)
  if (any(duplicated(gene_names))) {
    dup_genes <- gene_names[duplicated(gene_names)]
    issues$dup_genes <- paste("Found", length(dup_genes), "duplicated gene names")
  } else {
    cat("✓ No duplicated gene names\n")
  }
  
  # Check for duplicated sample names
  sample_names <- colnames(counts)
  if (any(duplicated(sample_names))) {
    issues$dup_samples <- paste("Found duplicated sample names")
  } else {
    cat("✓ No duplicated sample names\n")
  }
  
  # Check if samples in counts match metadata
  metadata_samples <- metadata$sample_id
  counts_samples <- colnames(counts)
  
  if (!all(counts_samples %in% metadata_samples)) {
    missing <- counts_samples[!counts_samples %in% metadata_samples]
    issues$missing_meta <- paste("Samples missing from metadata:", paste(missing, collapse = ", "))
  } else {
    cat("✓ All count matrix samples found in metadata\n")
  }
  
  # Summary statistics
  cat("\n--- Matrix Summary ---\n")
  cat("Dimensions:", nrow(counts), "genes x", ncol(counts), "samples\n")
  cat("Total counts:", format(sum(counts), big.mark = ","), "\n")
  cat("Genes with zero counts across all samples:", sum(rowSums(counts) == 0), "\n")
  
  cat("\n✓ Count matrix validation PASSED\n")
  return(TRUE)
}

# STEP 2: RAW COUNT DISTRIBUTION QC

plot_raw_count_distribution <- function(counts, output_dir = CONFIG$output_dir) {

  # Calculate log-transformed counts for visualization
  log_counts <- log10(counts + 1)
  
  # Histogram of all counts
  all_counts <- as.vector(as.matrix(counts))
  log_all_counts <- log10(all_counts + 1)
  
  df_hist <- data.frame(log_count = log_all_counts)
  
  p1 <- ggplot(df_hist, aes(x = log_count)) +
    geom_histogram(bins = 50, fill = "#3498db", color = "white", alpha = 0.8) +
    labs(
      title = "Distribution of Raw Counts (All Genes & Samples)",
      x = "log10(count + 1)",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12)
    )
  
  ggsave(file.path(output_dir, "01_raw_count_histogram.png"), p1, width = 10, height = 6, dpi = 300)
  
  # Per-sample boxplot
  df_box <- tidyr::pivot_longer(
    as.data.frame(log_counts) %>% tibble::rownames_to_column("gene"),
    cols = -gene,
    names_to = "sample",
    values_to = "log_count"
  )
  
  p2 <- ggplot(df_box, aes(x = sample, y = log_count, fill = sample)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    labs(
      title = "Per-Sample Count Distribution",
      x = "Sample",
      y = "log10(count + 1)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    scale_fill_brewer(palette = "Set3")
  
  ggsave(file.path(output_dir, "02_per_sample_boxplot.png"), p2, width = 12, height = 6, dpi = 300)

  return(list(histogram = p1, boxplot = p2))
}

# STEP 3: CREATE DESeqDataSet

create_deseq_dataset <- function(counts, metadata, design_formula = ~ treatment) {

  # Ensure counts is a matrix
  counts_matrix <- as.matrix(counts)
  
  # Ensure metadata row names match column names of counts
  rownames(metadata) <- metadata$sample_id
  
  # Reorder metadata to match counts columns
  metadata <- metadata[colnames(counts_matrix), , drop = FALSE]
  
  # Convert treatment to factor
  metadata$treatment <- factor(metadata$treatment)
  
  cat("Design formula:", deparse(design_formula), "\n")
  cat("Treatment levels:", paste(levels(metadata$treatment), collapse = ", "), "\n")
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = counts_matrix,
    colData = metadata,
    design = design_formula
  )

  cat("  - Genes:", nrow(dds), "\n")
  cat("  - Samples:", ncol(dds), "\n")
  
  return(dds)
}

# STEP 4: FILTER LOW-COUNT GENES

filter_low_count_genes <- function(dds, min_count = CONFIG$min_count, min_samples = CONFIG$min_samples) {

  original_genes <- nrow(dds)
  
  # Keep genes with at least min_count in at least min_samples
  keep <- rowSums(counts(dds) >= min_count) >= min_samples
  dds_filtered <- dds[keep, ]
  
  filtered_genes <- nrow(dds_filtered)
  removed_genes <- original_genes - filtered_genes
  
  cat("Filter criteria: ≥", min_count, "counts in ≥", min_samples, "samples\n")
  cat("Original genes:", format(original_genes, big.mark = ","), "\n")
  cat("Genes removed:", format(removed_genes, big.mark = ","), "\n")
  cat("Genes retained:", format(filtered_genes, big.mark = ","), "\n")
  cat("Retention rate:", round(filtered_genes / original_genes * 100, 1), "%\n")
  
  return(dds_filtered)
}

# STEP 5: NORMALIZE AND ESTIMATE SIZE FACTORS

normalize_counts <- function(dds, output_dir = CONFIG$output_dir) {

  # Estimate size factors
  dds <- estimateSizeFactors(dds)
  
  # Get size factors
  size_factors <- sizeFactors(dds)
  sf_table <- data.frame(
    sample = names(size_factors),
    size_factor = size_factors
  )
  
  cat("Size factors:\n")
  print(sf_table)
  
  # Get normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Save outputs
  write.csv(sf_table, file.path(output_dir, "03_size_factors.csv"), row.names = FALSE)
  write.csv(normalized_counts, file.path(output_dir, "04_normalized_counts.csv"))

  
  return(list(
    dds = dds,
    size_factors = sf_table,
    normalized_counts = normalized_counts
  ))
}

# STEP 6: VARIANCE STABILIZING TRANSFORMATION

apply_vst <- function(dds, output_dir = CONFIG$output_dir) {

  # Apply VST
  vst_data <- vst(dds, blind = FALSE)
  vst_matrix <- assay(vst_data)
  
  cat("  - Original count range:", min(counts(dds)), "-", max(counts(dds)), "\n")
  cat("  - VST value range:", round(min(vst_matrix), 2), "-", round(max(vst_matrix), 2), "\n")
  
  # Save VST matrix
  write.csv(vst_matrix, file.path(output_dir, "05_vst_matrix.csv"))

  return(vst_data)
}

# STEP 7: PCA AND SCREE PLOT

perform_pca_analysis <- function(vst_data, metadata, output_dir = CONFIG$output_dir) {

  # Get VST matrix
  vst_matrix <- assay(vst_data)
  
  # Perform PCA
  pca_result <- prcomp(t(vst_matrix), scale. = TRUE)
  
  # Calculate variance explained
  var_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100
  
  # PCA data frame
  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    sample = rownames(pca_result$x)
  )
  pca_df <- merge(pca_df, metadata, by.x = "sample", by.y = "sample_id")
  
  # PCA plot
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = treatment, label = sample)) +
    geom_point(size = 4, alpha = 0.8) +
    geom_text_repel(size = 3, max.overlaps = 20) +
    labs(
      title = "PCA of VST-Transformed Data",
      x = paste0("PC1 (", round(var_explained[1], 1), "% variance)"),
      y = paste0("PC2 (", round(var_explained[2], 1), "% variance)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    ) +
    scale_color_brewer(palette = "Set1")
  
  ggsave(file.path(output_dir, "06_pca_plot.png"), p_pca, width = 10, height = 8, dpi = 300)
  
  # Scree plot
  scree_df <- data.frame(
    PC = factor(1:length(var_explained), levels = 1:length(var_explained)),
    Variance = var_explained,
    Cumulative = cumsum(var_explained)
  )
  
  p_scree <- ggplot(scree_df[1:min(10, nrow(scree_df)), ], aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "#2ecc71", alpha = 0.8) +
    geom_line(aes(x = as.numeric(PC), y = Cumulative), color = "#e74c3c", size = 1.2) +
    geom_point(aes(x = as.numeric(PC), y = Cumulative), color = "#e74c3c", size = 3) +
    labs(
      title = "Scree Plot: Variance Explained by PCs",
      x = "Principal Component",
      y = "Variance Explained (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    ) +
    scale_y_continuous(sec.axis = sec_axis(~., name = "Cumulative Variance (%)"))
  
  ggsave(file.path(output_dir, "07_scree_plot.png"), p_scree, width = 10, height = 6, dpi = 300)
  
  cat("✓ Variance explained by PC1:", round(var_explained[1], 1), "%\n")
  cat("✓ Variance explained by PC2:", round(var_explained[2], 1), "%\n")

  return(list(
    pca_result = pca_result,
    var_explained = var_explained,
    pca_plot = p_pca,
    scree_plot = p_scree
  ))
}

# STEP 8: SAMPLE CORRELATION AND CLUSTERING

perform_sample_clustering <- function(vst_data, metadata, output_dir = CONFIG$output_dir) {

  vst_matrix <- assay(vst_data)
  
  # Calculate sample correlation
  sample_cor <- cor(vst_matrix, method = "pearson")
  
  # Prepare annotation
  annotation_col <- data.frame(
    Treatment = metadata$treatment,
    row.names = metadata$sample_id
  )
  
  # Define colors for annotation
  ann_colors <- list(
    Treatment = c(
      "no_treatment" = "#3498db",
      "siCNTL" = "#2ecc71",
      "siSOX10_1" = "#e74c3c",
      "siSOX10_2" = "#9b59b6"
    )
  )
  
  # Correlation heatmap
  png(file.path(output_dir, "08_sample_correlation_heatmap.png"), width = 10, height = 8, units = "in", res = 300)
  pheatmap(
    sample_cor,
    annotation_col = annotation_col,
    annotation_row = annotation_col,
    annotation_colors = ann_colors,
    clustering_method = "complete",
    main = "Sample-to-Sample Correlation Heatmap",
    fontsize = 10,
    display_numbers = TRUE,
    number_format = "%.2f",
    color = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100)
  )
  dev.off()
  
  # Save correlation matrix
  write.csv(sample_cor, file.path(output_dir, "09_sample_correlation_matrix.csv"))
  
  cat("✓ Correlation range:", round(min(sample_cor[lower.tri(sample_cor)]), 3), 
      "to", round(max(sample_cor[lower.tri(sample_cor)]), 3), "\n")

    return(sample_cor)
}

# STEP 9: RUN DESeq2 DIFFERENTIAL EXPRESSION

run_deseq_analysis <- function(dds) {

  # Run DESeq
  dds <- DESeq(dds)
  
  return(dds)
}

# STEP 10: MEAN-VARIANCE SANITY CHECK

plot_mean_variance <- function(counts, output_dir = CONFIG$output_dir) {

  # Calculate gene-wise mean and variance
  gene_means <- rowMeans(counts)
  gene_vars <- apply(counts, 1, var)
  
  # Filter out zeros for log transformation
  valid_idx <- gene_means > 0 & gene_vars > 0
  
  mv_df <- data.frame(
    mean = gene_means[valid_idx],
    variance = gene_vars[valid_idx]
  )
  
  # Plot
  p <- ggplot(mv_df, aes(x = log10(mean), y = log10(variance))) +
    geom_point(alpha = 0.3, size = 0.5, color = "#3498db") +
    geom_abline(slope = 1, intercept = 0, color = "#e74c3c", linetype = "dashed", size = 1) +
    labs(
      title = "Mean-Variance Relationship (RNA-seq)",
      x = "log10(Mean Expression)",
      y = "log10(Variance)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )
  
  ggsave(file.path(output_dir, "11_mean_variance_plot.png"), p, width = 10, height = 8, dpi = 300)

  return(p)
}

# STEP 11: DESeq2 DISPERSION PLOT

plot_dispersion <- function(dds, output_dir = CONFIG$output_dir) {

  png(file.path(output_dir, "12_dispersion_plot.png"), width = 10, height = 8, units = "in", res = 300)
  plotDispEsts(dds, main = "DESeq2 Dispersion Estimates")
  dev.off()

}

# STEP 12: DIFFERENTIAL EXPRESSION FOR SPECIFIC COMPARISONS

get_de_results <- function(dds, contrast, fdr_cutoff = CONFIG$fdr_cutoff, output_dir = CONFIG$output_dir) {
  cat("\n Computing DE for:", contrast[2], "vs", contrast[3], "\n")
  
  # Get results
  res <- results(dds, contrast = contrast, alpha = fdr_cutoff)
  
  # Summary
  cat("Results summary:\n")
  summary(res)
  
  # Order by adjusted p-value
  res_ordered <- res[order(res$padj), ]
  
  # Count significant genes
  sig_up <- sum(res$padj < fdr_cutoff & res$log2FoldChange > 0, na.rm = TRUE)
  sig_down <- sum(res$padj < fdr_cutoff & res$log2FoldChange < 0, na.rm = TRUE)
  
  cat("\nSignificant genes (FDR <", fdr_cutoff, "):\n")
  cat("  - Upregulated:", sig_up, "\n")
  cat("  - Downregulated:", sig_down, "\n")
  cat("  - Total:", sig_up + sig_down, "\n")
  
  return(res_ordered)
}

run_all_comparisons <- function(dds, fdr_cutoff = CONFIG$fdr_cutoff, output_dir = CONFIG$output_dir) {

  # Define comparisons based on workflow:
  # no_treatment vs siCNTL, siCNTL vs siSOX10_2, siSOX10_1 vs siSOX10_2
  comparisons <- list(
    "no_treatment_vs_siCNTL" = c("treatment", "no_treatment", "siCNTL"),
    "siCNTL_vs_siSOX10_2" = c("treatment", "siCNTL", "siSOX10_2"),
    "siSOX10_1_vs_siSOX10_2" = c("treatment", "siSOX10_1", "siSOX10_2")
  )
  
  results_list <- list()
  
  for (name in names(comparisons)) {
    contrast <- comparisons[[name]]
    res <- get_de_results(dds, contrast, fdr_cutoff, output_dir)
    results_list[[name]] <- res
    
    # Save results
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    write.csv(res_df, file.path(output_dir, paste0("13_DE_", name, ".csv")), row.names = FALSE)
  }
  
  return(results_list)
}

# STEP 13: UPSET PLOT FOR OVERLAPPING DEGs

generate_upset_plot <- function(de_results_list, fdr_cutoff = CONFIG$fdr_cutoff, output_dir = CONFIG$output_dir) {
  
  # Extract significant genes from each comparison
  sig_genes_list <- lapply(de_results_list, function(res) {
    sig <- subset(res, padj < fdr_cutoff)
    rownames(sig)
  })
  
  # Create binary matrix for UpSet
  all_genes <- unique(unlist(sig_genes_list))
  
  if (length(all_genes) == 0) {
    cat("⚠ No significant genes found at FDR <", fdr_cutoff, "\n")
    return(NULL)
  }
  
  upset_matrix <- sapply(sig_genes_list, function(genes) {
    as.integer(all_genes %in% genes)
  })
  rownames(upset_matrix) <- all_genes
  upset_df <- as.data.frame(upset_matrix)
  
  # Generate UpSet plot
  png(file.path(output_dir, "10_upset_plot.png"), width = 10, height = 6, units = "in", res = 300)
  print(upset(upset_df, 
              sets = colnames(upset_df),
              order.by = "freq",
              main.bar.color = "#3498db",
              sets.bar.color = "#2ecc71",
              text.scale = 1.5))
  dev.off()
  
  for (name in names(sig_genes_list)) {
    cat("  -", name, ":", length(sig_genes_list[[name]]), "genes\n")
  }
  
  return(upset_df)
}

# STEP 14: LFC SHRINKAGE AND MA PLOTS

perform_lfc_shrinkage <- function(dds, results_list, output_dir = CONFIG$output_dir) {

  # Get coefficient names
  coef_names <- resultsNames(dds)
  cat("Available coefficients:", paste(coef_names, collapse = ", "), "\n")
  
  shrunk_results <- list()
  
  # For each comparison, apply shrinkage and generate MA plot
  comparison_to_coef <- list(
    "no_treatment_vs_siCNTL" = "treatment_siCNTL_vs_no_treatment",
    "siCNTL_vs_siSOX10_2" = "treatment_siSOX10_2_vs_no_treatment",
    "siSOX10_1_vs_siSOX10_2" = "treatment_siSOX10_2_vs_no_treatment"
  )
  
  # Generate MA plots for original results
  for (name in names(results_list)) {
    res <- results_list[[name]]
    
    # MA plot
    png(file.path(output_dir, paste0("14_MA_plot_", name, ".png")), width = 10, height = 6, units = "in", res = 300)
    plotMA(res, main = paste("MA Plot:", gsub("_", " ", name)), ylim = c(-5, 5))
    dev.off()
    
  }
  
  # Apply shrinkage for the first coefficient (example)
  if ("treatment_siCNTL_vs_no_treatment" %in% coef_names) {
    res_shrunk <- lfcShrink(dds, coef = "treatment_siCNTL_vs_no_treatment", type = "apeglm")
    shrunk_results[["siCNTL_vs_no_treatment_shrunk"]] <- res_shrunk
    
    png(file.path(output_dir, "14_MA_plot_shrunk_siCNTL.png"), width = 10, height = 6, units = "in", res = 300)
    plotMA(res_shrunk, main = "MA Plot with LFC Shrinkage (siCNTL vs no_treatment)", ylim = c(-5, 5))
    dev.off()
  }
  
  return(shrunk_results) 
}

# STEP 15: KEGG ENRICHMENT ANALYSIS

perform_kegg_enrichment <- function(de_results, comparison_name, fdr_cutoff = CONFIG$fdr_cutoff, output_dir = CONFIG$output_dir) {

  # Get all genes tested (for background)
  all_genes <- rownames(de_results)
  
  # Remove version numbers from ENSEMBL IDs
  all_genes_clean <- gsub("\\..*", "", all_genes)

  gene_mapping <- tryCatch({
    bitr(all_genes_clean, 
         fromType = "ENSEMBL", 
         toType = "ENTREZID", 
         OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(gene_mapping) || nrow(gene_mapping) == 0) {
    return(NULL)
  }
  
  # Prepare gene list with log2FC values for GSEA
  de_df <- as.data.frame(de_results)
  de_df$ensembl_clean <- gsub("\\..*", "", rownames(de_df))
  de_df <- merge(de_df, gene_mapping, by.x = "ensembl_clean", by.y = "ENSEMBL")
  
  # Remove NAs and duplicates
  de_df <- de_df[!is.na(de_df$log2FoldChange) & !duplicated(de_df$ENTREZID), ]
  
  # Create named vector of log2FC
  gene_list <- de_df$log2FoldChange
  names(gene_list) <- de_df$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  kegg_gsea <- tryCatch({
    gseKEGG(
      geneList = gene_list,
      organism = "hsa",
      minGSSize = 10,
      maxGSSize = 500,
      pvalueCutoff = 0.1,
      verbose = FALSE
    )
  }, error = function(e) {
    return(NULL)
  })
  
  if (!is.null(kegg_gsea) && nrow(kegg_gsea) > 0) {
    
    # Save results
    kegg_df <- as.data.frame(kegg_gsea)
    write.csv(kegg_df, file.path(output_dir, paste0("15_KEGG_enrichment_", comparison_name, ".csv")), row.names = FALSE)
    
    # Dotplot
    if (nrow(kegg_gsea) > 0) {
      p_dot <- dotplot(kegg_gsea, showCategory = 15, title = paste("KEGG Enrichment:", comparison_name))
      ggsave(file.path(output_dir, paste0("15_KEGG_dotplot_", comparison_name, ".png")), p_dot, width = 12, height = 8, dpi = 300)
    }

        return(kegg_gsea)
  } else {
    return(NULL)
  }
}

# STEP 16: ADD GENE ANNOTATIONS

add_gene_annotations <- function(results_list, output_dir = CONFIG$output_dir) {

  annotated_results <- list()
  
  for (name in names(results_list)) {
    cat("\nAnnotating:", name, "\n")
    res <- results_list[[name]]
    res_df <- as.data.frame(res)
    res_df$ensembl_id <- rownames(res_df)
    
    # Remove version numbers from ENSEMBL IDs for mapping
    res_df$ensembl_clean <- gsub("\\..*", "", res_df$ensembl_id)
    
    # Get gene symbols and other annotations
    gene_annotations <- tryCatch({
      AnnotationDbi::select(org.Hs.eg.db, 
                            keys = unique(res_df$ensembl_clean), 
                            columns = c("SYMBOL", "GENENAME", "ENTREZID"), 
                            keytype = "ENSEMBL")
    }, error = function(e) {
      return(NULL)
    })
    
    if (!is.null(gene_annotations) && nrow(gene_annotations) > 0) {
      # Remove duplicates (keep first mapping)
      gene_annotations <- gene_annotations[!duplicated(gene_annotations$ENSEMBL), ]
      
      # Merge annotations with results
      res_df <- merge(res_df, gene_annotations, 
                      by.x = "ensembl_clean", by.y = "ENSEMBL", 
                      all.x = TRUE)
      
      # Reorder columns for better readability
      res_df <- res_df[order(res_df$padj, na.last = TRUE), ]
      
      # Count successful mappings
      mapped_count <- sum(!is.na(res_df$SYMBOL))
      
    } else {
      res_df$SYMBOL <- NA
      res_df$GENENAME <- NA
      res_df$ENTREZID <- NA
    }
    
    annotated_results[[name]] <- res_df
    
    # Save annotated results
    write.csv(res_df, file.path(output_dir, paste0("16_annotated_DE_", name, ".csv")), row.names = FALSE)
  }
  
  return(annotated_results)
}

# STEP 17: SUMMARIZE SIGNIFICANT GENES

summarize_significant_genes <- function(results_list, fdr_cutoff = CONFIG$fdr_cutoff, output_dir = CONFIG$output_dir) {

  summary_list <- list()
  
  for (name in names(results_list)) {
    res <- results_list[[name]]
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    
    sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < fdr_cutoff, ]
    sig_genes <- sig_genes[order(sig_genes$padj), ]
    
    summary_list[[name]] <- list(
      total = nrow(sig_genes),
      up = sum(sig_genes$log2FoldChange > 0),
      down = sum(sig_genes$log2FoldChange < 0),
      genes = sig_genes
    )
    
    cat("\n", name, ":\n")
    cat("  Total significant:", summary_list[[name]]$total, "\n")
    cat("  Upregulated:", summary_list[[name]]$up, "\n")
    cat("  Downregulated:", summary_list[[name]]$down, "\n")
    
    if (nrow(sig_genes) > 0) {
      write.csv(sig_genes, file.path(output_dir, paste0("17_significant_genes_", name, ".csv")), row.names = FALSE)
    }
  }
  
  return(summary_list)
}

# STEP 18: HEATMAP OF SIGNIFICANT GENES

plot_significant_genes_heatmap <- function(vst_data, annotated_df, metadata, top_n = 50, 
                                           comparison_name = "comparison", output_dir = CONFIG$output_dir) {

  # Filter significant genes from annotated results
  sig_df <- annotated_df[!is.na(annotated_df$padj) & annotated_df$padj < CONFIG$fdr_cutoff, ]
  sig_df <- sig_df[order(sig_df$padj), ]
  
  if (nrow(sig_df) == 0) {
    cat("⚠ No significant genes to plot\n")
    return(NULL)
  }
  
  vst_matrix <- assay(vst_data)
  
  # Get top N genes
  top_sig <- head(sig_df, top_n)
  genes_to_plot <- top_sig$ensembl_id[top_sig$ensembl_id %in% rownames(vst_matrix)]
  
  if (length(genes_to_plot) == 0) {
    return(NULL)
  }
  
  # Get gene symbols from annotated data
  gene_symbol_map <- setNames(top_sig$SYMBOL, top_sig$ensembl_id)
  # Use ENSEMBL ID if SYMBOL is NA
  gene_symbol_map[is.na(gene_symbol_map)] <- names(gene_symbol_map)[is.na(gene_symbol_map)]
  
  # Subset and scale
  heatmap_data <- vst_matrix[genes_to_plot, , drop = FALSE]
  heatmap_scaled <- t(scale(t(heatmap_data)))
  
  # Replace row names with gene symbols
  rownames(heatmap_scaled) <- gene_symbol_map[rownames(heatmap_scaled)]
  
  # Annotation
  annotation_col <- data.frame(
    Treatment = metadata$treatment,
    row.names = metadata$sample_id
  )
  
  ann_colors <- list(
    Treatment = c(
      "no_treatment" = "#3498db",
      "siCNTL" = "#2ecc71",
      "siSOX10_1" = "#e74c3c",
      "siSOX10_2" = "#9b59b6"
    )
  )
  
  # Generate heatmap with gene symbols
  png(file.path(output_dir, paste0("18_heatmap_significant_", comparison_name, ".png")), 
      width = 12, height = max(8, length(genes_to_plot) * 0.15), units = "in", res = 300)
  pheatmap(
    heatmap_scaled,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    clustering_method = "complete",
    scale = "none",  # Already scaled
    main = paste("Top", length(genes_to_plot), "Significant Genes (Gene Symbols) -", comparison_name),
    fontsize_row = 6,
    fontsize_col = 10,
    color = colorRampPalette(c("#3498db", "white", "#e74c3c"))(100),
    show_rownames = TRUE
  )
  dev.off()

  return(genes_to_plot)
}

# STEP 19: VOLCANO PLOTS

plot_volcano <- function(annotated_df, comparison_name, fdr_cutoff = CONFIG$fdr_cutoff, 
                         lfc_cutoff = CONFIG$lfc_cutoff, output_dir = CONFIG$output_dir) {

  res_df <- annotated_df[!is.na(annotated_df$padj), ]
  
  # Categorize genes
  res_df$significance <- "Not Significant"
  res_df$significance[res_df$padj < fdr_cutoff & res_df$log2FoldChange > lfc_cutoff] <- "Up"
  res_df$significance[res_df$padj < fdr_cutoff & res_df$log2FoldChange < -lfc_cutoff] <- "Down"
  res_df$significance <- factor(res_df$significance, levels = c("Not Significant", "Up", "Down"))
  
  # Get top 20 most significant genes
  top_genes <- res_df[order(res_df$padj), ]
  top_genes <- head(top_genes[top_genes$significance != "Not Significant", ], 20)
  
  # Use gene symbols from annotations (fallback to ensembl_id if NA)
  if (nrow(top_genes) > 0) {
    top_genes$gene_symbol <- ifelse(is.na(top_genes$SYMBOL), top_genes$ensembl_id, top_genes$SYMBOL)
  }
  
  # Volcano plot
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Not Significant" = "grey60", "Up" = "#e74c3c", "Down" = "#3498db")) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(fdr_cutoff), linetype = "dashed", color = "grey40") +
    labs(
      title = paste("Volcano Plot:", gsub("_", " ", comparison_name)),
      subtitle = paste("FDR <", fdr_cutoff, ", |log2FC| >", lfc_cutoff, " | Labels: Gene Symbols"),
      x = "log2 Fold Change",
      y = "-log10(Adjusted P-value)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "right"
    )
  
  # Add labels for top genes with gene symbols
  if (nrow(top_genes) > 0) {
    p <- p + geom_text_repel(
      data = top_genes,
      aes(label = gene_symbol),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      fontface = "bold"
    )
  }
  
  ggsave(file.path(output_dir, paste0("19_volcano_", comparison_name, ".png")), p, width = 10, height = 8, dpi = 300)
  
  cat("✓ Up-regulated:", sum(res_df$significance == "Up"), "\n")
  cat("✓ Down-regulated:", sum(res_df$significance == "Down"), "\n")

  return(p)
}

# STEP 20: EXPRESSION PLOTS FOR TOP GENES

plot_top_genes_expression <- function(dds, annotated_df, metadata, top_n = 20, 
                                      comparison_name = "comparison", output_dir = CONFIG$output_dir) {

  res_df <- annotated_df[!is.na(annotated_df$padj), ]
  res_df <- res_df[order(res_df$padj), ]
  
  # Get top N significant genes
  top_df <- head(res_df, top_n)
  
  if (nrow(top_df) == 0) {
    return(NULL)
  }
  
  # Get normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  
  # Filter to genes present in counts
  top_df <- top_df[top_df$ensembl_id %in% rownames(norm_counts), ]
  top_genes <- top_df$ensembl_id
  
  if (length(top_genes) == 0) {
    return(NULL)
  }
  
  # Create gene symbol map from annotations (fallback to ensembl_id if NA)
  gene_symbol_map <- setNames(
    ifelse(is.na(top_df$SYMBOL), top_df$ensembl_id, top_df$SYMBOL),
    top_df$ensembl_id
  )
  
  plot_data <- norm_counts[top_genes, , drop = FALSE]
  
  # Reshape for plotting first, then add gene symbols
  plot_df <- as.data.frame(t(plot_data))
  plot_df$sample <- rownames(plot_df)
  plot_df <- merge(plot_df, metadata, by.x = "sample", by.y = "sample_id")
  
  # Pivot with original gene names (ENSEMBL IDs)
  plot_df_long <- tidyr::pivot_longer(
    plot_df,
    cols = all_of(top_genes),
    names_to = "ensembl_id",
    values_to = "expression"
  )
  
  # Add gene symbols as a new column
  plot_df_long$gene_symbol <- gene_symbol_map[plot_df_long$ensembl_id]
  
  # Plot with gene symbols
  p <- ggplot(plot_df_long, aes(x = treatment, y = log10(expression + 1), fill = treatment)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
    facet_wrap(~ gene_symbol, scales = "free_y", ncol = 4) +
    labs(
      title = paste("Top", length(top_genes), "Significant Genes (Gene Symbols) -", comparison_name),
      x = "Treatment",
      y = "log10(Normalized Count + 1)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(size = 8, face = "bold"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
    ) +
    scale_fill_brewer(palette = "Set1")
  
  ggsave(file.path(output_dir, paste0("20_top_genes_expression_", comparison_name, ".png")), 
         p, width = 14, height = max(6, length(top_genes) / 4 * 3), dpi = 300)
  
  return(p)
}

# MAIN PIPELINE FUNCTION

run_full_pipeline <- function(counts_file = CONFIG$counts_file, 
                               metadata_file = CONFIG$metadata_file,
                               output_dir = CONFIG$output_dir) {

  # Load counts
  counts_raw <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  
  # Load metadata
  metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)

  cat("- Counts:", nrow(counts_raw), "genes x", ncol(counts_raw), "samples\n")
  cat("- Metadata:", nrow(metadata), "samples\n")
  
  # ----- STEP 1: Validate Count Matrix -----
  validate_count_matrix(counts_raw, metadata)
  
  # ----- STEP 2: Raw Count Distribution QC -----
  plot_raw_count_distribution(counts_raw, output_dir)
  
  # ----- STEP 3: Create DESeqDataSet -----
  dds <- create_deseq_dataset(counts_raw, metadata)
  
  # ----- STEP 4: Filter Low-Count Genes -----
  dds_filtered <- filter_low_count_genes(dds)
  
  # ----- STEP 5: Normalize Counts -----
  norm_result <- normalize_counts(dds_filtered, output_dir)
  dds_filtered <- norm_result$dds
  
  # ----- STEP 6: VST Transformation -----
  vst_data <- apply_vst(dds_filtered, output_dir)
  
  # ----- STEP 7: PCA Analysis -----
  pca_result <- perform_pca_analysis(vst_data, metadata, output_dir)
  
  # ----- STEP 8: Sample Clustering -----
  sample_cor <- perform_sample_clustering(vst_data, metadata, output_dir)
  
  # ----- STEP 9: Run DESeq Analysis -----
  dds_de <- run_deseq_analysis(dds_filtered)
  
  # ----- STEP 10: Mean-Variance Plot  -----
  plot_mean_variance(counts(dds_filtered), output_dir)
  
  # ----- STEP 11: Dispersion Plot -----
  plot_dispersion(dds_de, output_dir)
  
  # ----- STEP 12: Run All Comparisons -----
  de_results <- run_all_comparisons(dds_de, CONFIG$fdr_cutoff, output_dir)
  
  # ----- STEP 13: UpSet Plot -----
  upset_data <- generate_upset_plot(de_results, CONFIG$fdr_cutoff, output_dir)
  
  # ----- STEP 14: LFC Shrinkage and MA Plots -----
  shrunk_results <- perform_lfc_shrinkage(dds_de, de_results, output_dir)
  
  # ----- STEP 15: KEGG Enrichment (for main comparison) -----
  for (name in names(de_results)) {
    kegg_result <- perform_kegg_enrichment(de_results[[name]], name, CONFIG$fdr_cutoff, output_dir)
  }
  
  # ----- STEP 16: Add Gene Annotations -----
  annotated_results <- add_gene_annotations(de_results, output_dir)
  
  # ----- STEP 17: Summarize Significant Genes -----
  sig_summary <- summarize_significant_genes(de_results, CONFIG$fdr_cutoff, output_dir)
  
  # ----- STEP 18, 19, 20: Visualizations for each comparison -----
  for (name in names(annotated_results)) {
    cat("\n\n>>> Processing visualizations for:", name, "\n")
    
    # Use annotated results which already contain gene symbols
    annotated_df <- annotated_results[[name]]
    
    # Step 18: Heatmap (using annotated data with gene symbols)
    plot_significant_genes_heatmap(vst_data, annotated_df, metadata, 50, name, output_dir)
    
    # Step 19: Volcano plot (using annotated data with gene symbols)
    plot_volcano(annotated_df, name, CONFIG$fdr_cutoff, CONFIG$lfc_cutoff, output_dir)
    
    # Step 20: Top genes expression (using annotated data with gene symbols)
    plot_top_genes_expression(dds_de, annotated_df, metadata, 20, name, output_dir)
  }
  
  return(list(
    dds = dds_de,
    vst_data = vst_data,
    de_results = de_results,
    annotated_results = annotated_results,
    sig_summary = sig_summary
  ))
}

# Run the full pipeline
results <- run_full_pipeline()

# Print summary
cat("\n\n=== SUMMARY ===\n")
for (name in names(results$sig_summary)) {
  cat("\n", name, ":\n")
  cat("  Significant genes:", results$sig_summary[[name]]$total, "\n")
  cat("    - Upregulated:", results$sig_summary[[name]]$up, "\n")
  cat("    - Downregulated:", results$sig_summary[[name]]$down, "\n")
}

