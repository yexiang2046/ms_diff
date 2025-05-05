# Proteomics Differential Expression Analysis Pipeline

This pipeline performs differential expression analysis on proteomics data from MaxQuant output, specifically analyzing protein groups data.

## Overview

The pipeline processes MaxQuant proteinGroups.txt files and performs:
- Quality control analysis
- Data normalization using VSN
- Missing value imputation
- Differential expression analysis
- Visualization of results

## Requirements

The pipeline is containerized using Docker. You only need:
- Docker installed on your system
- Input data files (proteinGroups.txt from MaxQuant)

## Quick Start

1. Pull the Docker image:
```bash
docker pull xiang2019/ms_diff:v1.0.0
```

2. modify config/analysis_config.R to add your input data file
   protein_groups = "your maxquant data file proteinGroups.txt"

2. Run the analysis:
```bash
docker run -v /path/to/your/data:/app -v /path/to/ms_diff/output:/app/output xiang2019/ms_diff:v1.0.0 Rscript /app/src/run_analysis.R
```

Replace `/path/to/your/data` with the directory containing your proteinGroups.txt file and `/path/to/output` with your desired output directory.

## Input Data

The pipeline expects:
- A proteinGroups.txt file from MaxQuant output
- The file should contain intensity columns starting with "Intensity."

## Output

The pipeline generates:
- Quality control plots in the `qc` directory:
  - Correlation plot
  - Missing values plot
  - Intensity distribution plot
  - Volcano plots for each comparison
- Results in the `results` directory:
  - Differential expression results for each comparison
  - Normalized and imputed data

## Configuration

The analysis parameters can be modified in `config/analysis_config.R`:
- FDR threshold for protein identification
- Sample groups and experimental design
- Output directories
- Plot parameters

## Analysis Steps

1. **Data Processing**
   - Filter contaminants and reverse hits
   - Apply FDR threshold
   - Log2 transform intensities

2. **Quality Control**
   - Sample correlation analysis
   - Missing value analysis
   - Intensity distribution analysis

3. **Data Normalization**
   - VSN normalization
   - Fallback to log2 normalization if VSN fails

4. **Missing Value Imputation**
   - Impute missing values using minimum detection method

5. **Differential Expression Analysis**
   - Linear model fitting
   - Empirical Bayes moderation
   - Multiple testing correction

6. **Visualization**
   - Volcano plots for each comparison
   - Quality control plots

## Troubleshooting

If you encounter any issues:
1. Ensure your input file is properly formatted
2. Check that you have sufficient disk space for output
3. Verify that the Docker container has proper permissions to access input/output directories

## License

This project is licensed under the MIT License - see the LICENSE file for details. 