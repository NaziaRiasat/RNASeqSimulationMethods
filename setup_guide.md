# Setup Guide: Improved RNA-seq Simulation Method

## Overview

This guide will walk you through setting up and using the improved RNA-seq simulation method that addresses the key limitations of existing simulators like SPsimSeq and Splatter. The method focuses on preserving gene-to-gene correlations and realistic biological structure.

## Prerequisites

### Required R Packages

```r
# Install required packages
install.packages(c("MASS", "Matrix", "igraph", "cluster", "ggplot2", 
                   "gridExtra", "reshape2", "corrplot", "devtools"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("edgeR", "DESeq2", "SPsimSeq"))
```

### System Requirements

- R version 4.0 or higher
- At least 8GB RAM (16GB recommended for large datasets)
- Sufficient disk space for intermediate files

## Quick Start

### 1. Load the Improved Simulation Method

```r
# Source the improved simulation functions
source("ImprovedRNASimR.R")

# Load your existing data (replace with your actual data paths)
Counts <- read.csv("project_read_counts.csv", row.names = 1)
Annot <- read.csv("gene_annotations.csv")
```

### 2. Prepare Your Reference Data

```r
# Create group factor (adjust based on your experimental design)
Group <- factor(c(rep("Nox2KO_PBS", 6), rep("Nox2KO_LPS", 6)))

# Filter and prepare data (similar to your existing workflow)
library(edgeR)
y <- DGEList(counts = Counts, genes = Annot, group = Group)
isexpr <- filterByExpr(y)
hasannot <- rowSums(is.na(y$genes)) == 0
y <- y[isexpr & hasannot, keep.lib.sizes = FALSE]

# Use the filtered count matrix as reference
reference_data <- y$counts
```

### 3. Run the Improved Simulation

```r
# Run improved simulation
improved_result <- ImprovedRNASeqSimulate(
  reference_data = reference_data,
  groups = Group,
  n_genes = 1000,           # Number of genes to simulate
  n_samples = 20,           # Number of samples to simulate
  max_iterations = 50       # Refinement iterations
)

# Access results
simulated_counts <- improved_result$counts
sample_metadata <- improved_result$sample_metadata
gene_metadata <- improved_result$gene_metadata
```

## Working with Your Existing Data

### Using Your Project Data

Here's how to adapt the method to your existing `SPSimseq.R` workflow:

```r
# Replace the SPsimSeq call in your existing code
# OLD:
# sim.data.bulk <- SPsimSeq(n.sim = 1, s.data = y$counts, group = Group, ...)

# NEW:
improved_sim_result <- ImprovedRNASeqSimulate(
  reference_data = y$counts,
  groups = Group,
  n_genes = 12600,
  n_samples = 172,
  max_iterations = 50
)

# Extract counts in the same format as SPsimSeq
simulated_counts <- improved_sim_result$counts
```

### Integration with Your PCA Analysis

```r
# Use the improved simulation results in your existing PCA workflow
# Replace 'sim.data.bulk1$counts' with 'improved_sim_result$counts'

# Create DGEList object with improved simulation data
Datta_improved <- DGEList(
  counts = improved_sim_result$counts,
  genes = improved_sim_result$gene_metadata,
  group = improved_sim_result$sample_metadata$Group
)

# Continue with your existing filtering and analysis
isexpr_improved <- filterByExpr(Datta_improved)
Datta_improved <- Datta_improved[isexpr_improved, keep.lib.sizes = FALSE]

# Apply TMM normalization
Datta_improved <- calcNormFactors(Datta_improved)

# Perform PCA
improved_data_t <- t(Datta_improved$counts)
pca_improved <- prcomp(improved_data_t, scale. = TRUE, center = TRUE)
```

## Quality Assessment

### Compare with Your Existing Results

```r
# Load your existing SPsimSeq results
load("your_existing_spsimseq_results.RData")  # Adjust path as needed

# Assess quality of both methods
spsim_quality <- assess_simulation_quality(reference_data, 
                                         list(counts = your_spsim_counts))

improved_quality <- assess_simulation_quality(reference_data, 
                                            improved_sim_result)

# Print comparison
print("Quality Comparison:")
print(paste("SPsimSeq RV coefficient:", round(spsim_quality$correlation_quality$rv_coefficient, 3)))
print(paste("Improved RV coefficient:", round(improved_quality$correlation_quality$rv_coefficient, 3)))
```

### Visualize Improvements

```r
# Create comparison plots
source("demonstration_script.R")

# Run comparison using your actual data
plots <- create_comparison_plots(
  real_data = reference_data,
  spsim_data = your_spsim_counts,
  improved_data = improved_sim_result$counts
)

# Save plots
ggsave("my_simulation_comparison.png", 
       grid.arrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2), 
       width = 12, height = 10)
```

## Advanced Configuration

### Customizing Parameters

```r
# For better correlation preservation
improved_result <- ImprovedRNASeqSimulate(
  reference_data = reference_data,
  groups = Group,
  n_genes = 1000,
  n_samples = 20,
  max_iterations = 100,        # More iterations for better accuracy
  correlation_threshold = 0.6,  # Adjust correlation detection threshold
  min_network_size = 5         # Minimum genes per network
)
```

### Batch Effect Simulation

```r
# If you have batch information
batch_info <- factor(rep(c("Batch1", "Batch2"), each = 6))

improved_result_batch <- ImprovedRNASeqSimulate(
  reference_data = reference_data,
  groups = Group,
  n_genes = 1000,
  n_samples = 20,
  batch_info = batch_info,
  max_iterations = 50
)
```

## Troubleshooting

### Common Issues

1. **Memory Issues**: Reduce the number of genes or samples if you encounter memory problems
2. **Slow Performance**: Decrease `max_iterations` or use a smaller reference dataset
3. **Correlation Matrix Issues**: The method automatically handles non-positive definite matrices

### Performance Tips

```r
# For large datasets, use parallel processing
library(parallel)
cl <- makeCluster(detectCores() - 1)
# Use cluster for parameter estimation (implementation depends on specific needs)
stopCluster(cl)

# Monitor memory usage
pryr::mem_used()
```

## Integration with Existing Workflows

### Replacing SPsimSeq in Your Pipeline

1. **Find SPsimSeq calls** in your existing code
2. **Replace with ImprovedRNASeqSimulate** calls
3. **Update result extraction** to use the new format
4. **Maintain compatibility** with downstream analyses

### Example Replacement

```r
# OLD CODE:
# sim.data.bulk1 <- SPsimSeq(n.sim = 1, s.data = y$counts, group = Group,
#                           n.genes = 1260, batch.config = 1, ...)
# counts.SPsimseq <- sim.data.bulk1$counts

# NEW CODE:
sim.data.bulk1 <- ImprovedRNASeqSimulate(
  reference_data = y$counts,
  groups = Group,
  n_genes = 1260,
  n_samples = 80,
  max_iterations = 50
)
counts.SPsimseq <- sim.data.bulk1$counts
```

## Validation and Testing

### Recommended Validation Steps

1. **Correlation Preservation**: Check that gene-gene correlations are better preserved
2. **Distributional Similarity**: Verify that mean and variance distributions match real data
3. **Downstream Analysis**: Test that differential expression analysis works appropriately
4. **Biological Realism**: Ensure pathway structure is maintained

### Running the Full Comparison

```r
# Run the complete demonstration
source("demonstration_script.R")

# This will:
# 1. Create reference data
# 2. Run both SPsimSeq and improved method
# 3. Compare quality metrics
# 4. Generate visualization plots
# 5. Save results for further analysis
```

## Next Steps

1. **Test with your actual data**: Replace the example data with your real RNA-seq dataset
2. **Optimize parameters**: Adjust simulation parameters based on your specific needs
3. **Integrate with analysis pipeline**: Update your existing analysis workflow
4. **Validate results**: Ensure the simulated data meets your quality requirements
5. **Scale up**: Apply to larger datasets once you're satisfied with the results

## Support and Documentation

- **Main simulation function**: `ImprovedRNASeqSimulate()`
- **Quality assessment**: `assess_simulation_quality()`
- **Visualization**: Functions in `demonstration_script.R`
- **Full documentation**: See `improved_rnaseq_simulation_plan.md`

## Example Complete Workflow

```r
# Load libraries and functions
library(edgeR)
library(ggplot2)
source("ImprovedRNASimR.R")

# Load and prepare your data
Counts <- read.csv("project_read_counts.csv", row.names = 1)
Group <- factor(c(rep("Nox2KO_PBS", 6), rep("Nox2KO_LPS", 6)))

# Create DGEList and filter
y <- DGEList(counts = Counts, group = Group)
isexpr <- filterByExpr(y)
y <- y[isexpr, keep.lib.sizes = FALSE]

# Run improved simulation
improved_result <- ImprovedRNASeqSimulate(
  reference_data = y$counts,
  groups = Group,
  n_genes = 1000,
  n_samples = 20,
  max_iterations = 50
)

# Print summary
print_simulation_summary(improved_result)

# Assess quality
quality_metrics <- assess_simulation_quality(y$counts[1:1000, ], improved_result)
print(quality_metrics)

# Save results
save(improved_result, quality_metrics, file = "improved_simulation_results.RData")
```

This completes the setup guide for the improved RNA-seq simulation method. The approach addresses the key limitations you identified in existing methods while maintaining compatibility with your existing analysis workflow.