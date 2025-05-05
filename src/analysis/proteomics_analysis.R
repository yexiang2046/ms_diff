#!/usr/bin/env Rscript

# Load required libraries
library(DEP)
library(PerformanceAnalytics)
library(vsn)
library(limma)
library(tidyverse)

#' Process MaxQuant proteinGroups data
#' @param protein_groups_file Path to the proteinGroups.txt file
#' @param fdr_threshold FDR threshold for filtering (default: 0.01)
#' @return Processed data frame
process_protein_groups <- function(protein_groups_file, fdr_threshold = 0.01) {
    # Read proteinGroups file
    data <- read.delim(protein_groups_file, stringsAsFactors = FALSE, sep = "\t")
    
    # Filter contaminants and reverse hits
    data <- data %>%
        filter(Potential.contaminant != "+") %>%
        filter(Reverse != "+")
    
    # Filter by FDR (using Q-value column)
    data <- data %>%
        filter(Q.value <= fdr_threshold)
    
    # Extract intensity columns
    intensity_cols <- grep("^Intensity\\.", names(data), value = TRUE)
    
    # Log2 transform intensities
    data[, intensity_cols] <- log2(data[, intensity_cols])
    
    # Add protein information columns
    data <- data %>%
        mutate(
            Protein.IDs = `Protein.IDs`,
            Gene.names = `Gene.names`,
            Protein.names = `Protein.names`,
            Sequence.coverage = `Sequence.coverage....`,
            Mol..weight = `Mol..weight..kDa.`,
            Sequence.length = `Sequence.length`
        )
    
    # Select relevant columns
    selected_cols <- c(
        "Protein.IDs", "Gene.names", "Protein.names",
        "Sequence.coverage", "Mol..weight", "Sequence.length",
        intensity_cols
    )
    
    data <- data[, selected_cols]
    
    return(data)
}

#' Perform quality control analysis
#' @param data Processed protein data
#' @param output_dir Directory for QC plots
perform_qc_analysis <- function(data, output_dir) {
    # Create output directory if it doesn't exist
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Extract intensity columns
    intensity_cols <- grep("^Intensity\\.", names(data), value = TRUE)
    
    # Correlation analysis
    cor_matrix <- cor(data[, intensity_cols], use = "pairwise.complete.obs")
    
    # Save correlation plot
    pdf(file.path(output_dir, "correlation_plot.pdf"))
    chart.Correlation(data[, intensity_cols], histogram = TRUE, pch = 19)
    dev.off()
    
    # Missing value analysis
    missing_values <- colSums(is.na(data[, intensity_cols]))
    missing_df <- data.frame(
        Sample = names(missing_values),
        Missing_Values = missing_values
    )
    
    # Save missing values plot
    pdf(file.path(output_dir, "missing_values.pdf"))
    print(ggplot(missing_df, aes(x = Sample, y = Missing_Values)) +
        geom_bar(stat = "identity") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = "Missing Values per Sample",
             y = "Number of Missing Values"))
    dev.off()
    
    # Distribution of protein intensities
    pdf(file.path(output_dir, "intensity_distribution.pdf"))
    data_long <- data %>%
        select(all_of(intensity_cols)) %>%
        pivot_longer(everything(), names_to = "Sample", values_to = "Intensity")
    
    print(ggplot(data_long, aes(x = Intensity, fill = Sample)) +
        geom_density(alpha = 0.5) +
        theme_minimal() +
        labs(title = "Distribution of Protein Intensities",
             x = "Log2 Intensity",
             y = "Density"))
    dev.off()
    
    # Log QC statistics
    cat(sprintf("Total proteins: %d\n", nrow(data)))
    cat(sprintf("Proteins with missing values: %d\n", sum(rowSums(is.na(data[, intensity_cols])) > 0)))
    cat(sprintf("Average missing values per sample: %.2f\n", mean(missing_values)))
}

#' Normalize protein intensities using VSN
#' @param data Processed protein data
#' @return Normalized data
normalize_intensities <- function(data) {
    # Extract intensity columns
    intensity_cols <- grep("^Intensity\\.", names(data), value = TRUE)
    
    # Replace infinite values with NA
    data[, intensity_cols] <- lapply(data[, intensity_cols], function(x) {
        x[is.infinite(x)] <- NA
        return(x)
    })
    
    # Remove rows with all NA values
    valid_rows <- rowSums(!is.na(data[, intensity_cols])) > 0
    data <- data[valid_rows, ]
    
    # Apply VSN normalization
    tryCatch({
        vsn_fit <- vsn2(as.matrix(data[, intensity_cols]))
        normalized_data <- predict(vsn_fit, newdata = as.matrix(data[, intensity_cols]))
        
        # Replace original intensities with normalized values
        data[, intensity_cols] <- normalized_data
        
        # Log the number of rows used for normalization
        cat(sprintf("Normalized %d proteins with valid intensity values\n", sum(valid_rows)))
    }, error = function(e) {
        cat("Error in VSN normalization, falling back to log2 normalization\n")
        # Fallback to log2 normalization
        data[, intensity_cols] <- log2(data[, intensity_cols])
    })
    
    return(data)
}

#' Impute missing values using Mindet method
#' @param data Normalized protein data
#' @return Data with imputed values
impute_missing_values <- function(data) {
    # Extract intensity columns
    intensity_cols <- grep("^Intensity\\.", names(data), value = TRUE)
    
    # Find minimum value for each protein
    min_values <- apply(data[, intensity_cols], 1, min, na.rm = TRUE)
    
    # Impute missing values with minimum value
    for (i in 1:nrow(data)) {
        data[i, intensity_cols][is.na(data[i, intensity_cols])] <- min_values[i]
    }
    
    # Log the number of imputed values
    total_na <- sum(is.na(data[, intensity_cols]))
    if (total_na > 0) {
        cat(sprintf("Imputed %d missing values\n", total_na))
    }
    
    return(data)
}

#' Perform differential expression analysis
#' @param data Normalized and imputed protein data
#' @param design_matrix Design matrix for the experiment
#' @param contrast_matrix Contrast matrix for the comparisons
#' @return Results of differential expression analysis
perform_de_analysis <- function(data, design_matrix, contrast_matrix) {
    # Extract intensity columns
    intensity_cols <- grep("^Intensity\\.", names(data), value = TRUE)
    
    # Create expression matrix
    expr_matrix <- as.matrix(data[, intensity_cols])
    
    # Fit linear model
    fit <- lmFit(expr_matrix, design_matrix)
    
    # Apply contrasts
    fit2 <- contrasts.fit(fit, contrast_matrix)
    
    # Apply empirical Bayes moderation
    fit2 <- eBayes(fit2)
    
    # Get results for each contrast
    results_list <- list()
    for (i in 1:ncol(contrast_matrix)) {
        contrast_name <- colnames(contrast_matrix)[i]
        results <- topTable(fit2, coef = i, adjust.method = "BH", number = Inf)
        
        # Add protein information
        results$Protein.IDs <- data$Protein.IDs
        results$Gene.names <- data$Gene.names
        results$Protein.names <- data$Protein.names
        
        # Rename columns to include contrast name
        colnames(results)[1:6] <- paste0(colnames(results)[1:6], "_", contrast_name)
        
        results_list[[contrast_name]] <- results
    }
    
    return(results_list)
}

#' Create volcano plots for differential expression results
#' @param results_list List of differential expression results
#' @param output_dir Directory for output files
create_volcano_plots <- function(results_list, output_dir) {
    for (contrast_name in names(results_list)) {
        results <- results_list[[contrast_name]]
        
        # Get the logFC and adj.P.Val columns for this contrast
        logfc_col <- grep(paste0("^logFC_", contrast_name), names(results), value = TRUE)
        pval_col <- grep(paste0("^adj.P.Val_", contrast_name), names(results), value = TRUE)
        
        if (length(logfc_col) > 0 && length(pval_col) > 0) {
            # Create volcano plot
            pdf(file.path(output_dir, paste0("volcano_plot_", contrast_name, ".pdf")))
            print(ggplot(results, aes_string(x = logfc_col, y = paste0("-log10(", pval_col, ")"))) +
                geom_point(aes_string(color = paste0(pval_col, " < 0.05"))) +
                theme_minimal() +
                labs(title = paste("Volcano Plot -", contrast_name),
                     x = "Log2 Fold Change",
                     y = "-log10 Adjusted P-value"))
            dev.off()
        }
    }
}

#' Main analysis pipeline
#' @param protein_groups_file Path to proteinGroups.txt file
#' @param output_dir Directory for output files
#' @param design_matrix Design matrix for the experiment
#' @param contrast_matrix Contrast matrix for the comparisons
main <- function(protein_groups_file, output_dir, design_matrix, contrast_matrix) {
    # Process protein groups
    cat("Processing protein groups...\n")
    data <- process_protein_groups(protein_groups_file)
    
    # Perform QC analysis
    cat("Performing quality control analysis...\n")
    perform_qc_analysis(data, file.path(output_dir, "qc"))
    
    # Normalize intensities
    cat("Normalizing intensities...\n")
    normalized_data <- normalize_intensities(data)
    
    # Impute missing values
    cat("Imputing missing values...\n")
    imputed_data <- impute_missing_values(normalized_data)
    
    # Perform differential expression analysis
    cat("Performing differential expression analysis...\n")
    de_results <- perform_de_analysis(imputed_data, design_matrix, contrast_matrix)
    
    # Create volcano plots
    cat("Creating volcano plots...\n")
    create_volcano_plots(de_results, file.path(output_dir, "qc"))
    
    # Save results
    cat("Saving results...\n")
    for (contrast_name in names(de_results)) {
        write.csv(de_results[[contrast_name]], 
                 file.path(output_dir, paste0("differential_expression_results_", contrast_name, ".csv")))
    }
    write.csv(imputed_data, file.path(output_dir, "normalized_imputed_data.csv"))
    
    cat("Analysis complete!\n")
}

# Example usage (commented out)
# design_matrix <- model.matrix(~0 + factor(c("Control", "Control", "Treatment", "Treatment")))
# contrast_matrix <- makeContrasts(Treatment - Control, levels = design_matrix)
# main("path/to/proteinGroups.txt", "output", design_matrix, contrast_matrix) 
