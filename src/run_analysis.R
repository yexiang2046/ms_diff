#!/usr/bin/env Rscript

# Load required libraries
library(DEP)
library(PerformanceAnalytics)
library(vsn)
library(limma)
library(tidyverse)

# Source analysis script and configuration
source("src/analysis/proteomics_analysis.R")
source("config/analysis_config.R")

# Run the analysis pipeline
main(
    protein_groups_file = FILE_PATHS$protein_groups,
    output_dir = OUTPUT_DIRS$base,
    design_matrix = DESIGN_MATRIX,
    contrast_matrix = CONTRAST_MATRIX
)

# Print completion message
cat("\nAnalysis pipeline completed successfully!\n")
cat("Results can be found in:", OUTPUT_DIRS$base, "\n")
cat("QC plots are available in:", OUTPUT_DIRS$qc, "\n")
cat("Differential expression results:", FILE_PATHS$results, "\n")
cat("Normalized and imputed data:", FILE_PATHS$normalized_data, "\n") 
