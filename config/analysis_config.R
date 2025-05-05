# Analysis configuration parameters

# FDR threshold for protein identification
FDR_THRESHOLD <- 0.01

# Sample groups
SAMPLE_GROUPS <- c(
    "GFPnegRNase", "GFPnegRNase",
    "K9negRNase", "K9negRNase",
    "GFPposRNase", "GFPposRNase",
    "K9posRNase", "K9posRNase"
)

# Create design matrix
DESIGN_MATRIX <- model.matrix(~0 + factor(SAMPLE_GROUPS))
colnames(DESIGN_MATRIX) <- c("GFPnegRNase", "K9negRNase", "GFPposRNase", "K9posRNase")

# Create contrast matrix for all pairwise comparisons
CONTRAST_MATRIX <- makeContrasts(
    K9negRNase_vs_GFPnegRNase = K9negRNase - GFPnegRNase,
    K9posRNase_vs_GFPposRNase = K9posRNase - GFPposRNase,
    K9posRNase_vs_K9negRNase = K9posRNase - K9negRNase,
    GFPposRNase_vs_GFPnegRNase = GFPposRNase - GFPnegRNase,
    levels = DESIGN_MATRIX
)

# Output directories
OUTPUT_DIRS <- list(
    base = "output",
    qc = "output/qc",
    results = "output/results"
)

# Create output directories
for (dir in OUTPUT_DIRS) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
}

# File paths
FILE_PATHS <- list(
    protein_groups = "proteinGroups_6361.txt",
    results = file.path(OUTPUT_DIRS$results, "differential_expression_results.csv"),
    normalized_data = file.path(OUTPUT_DIRS$results, "normalized_imputed_data.csv")
)

# Plot parameters
PLOT_PARAMS <- list(
    width = 10,
    height = 8,
    dpi = 300,
    theme = theme_minimal()
) 
