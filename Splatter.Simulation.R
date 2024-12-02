# Install and load necessary libraries
if (!requireNamespace("splatter", quietly = TRUE)) {
  install.packages("splatter")
}
library(splatter)


# Load necessary libraries (install them first if not already installed)
if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages("data.table")
}
library(data.table)


# Read the counts file (gene expression data)
Counts <- fread("project_read_counts.csv")

# Read the gene annotation file
Annot <- fread("gene_annotations.csv")

# Define the group factor for sample conditions
Group <- factor(c(
  rep("WT_PBS", 6), 
  rep("WT_LPS", 6), 
  rep("Nox2KO_PBS", 6), 
  rep("Nox2KO_LPS", 6)
))

# Inspect the data
print(head(Counts))
print(head(Annot))
print(Group)

# Estimate parameters from the Counts data
params <- splatEstimate(as.matrix(Counts1))

# Simulate RNA-seq bulk data using Splatter
simulated_bulk_splatter <- splatSimulate(
  params, 
  group.prob = c(0.25, 0.25, 0.25, 0.25),  # Proportions for each group
  method = "groups"
)

# Extract the simulated count matrix
simulated_counts_splatter <- counts(simulated_bulk_splatter)
head(simulated_counts_splatter)

# Save the simulated data
write.csv(simulated_counts_splatter, "simulated_bulk_splatter.csv", row.names = TRUE)

library(ggplot2)

# Convert real and simulated data to matrices
real_data <- as.matrix(Counts1)
simulated_data <- as.matrix(simulated_counts_splatter)

# Function to filter genes with zero variance
filter_genes_with_variance <- function(data) {
  # Calculate the variance for each gene
  gene_variance <- apply(data, 1, var)
  # Filter out genes with zero variance
  data_filtered <- data[gene_variance > 0, ]
  return(data_filtered)
}

# Filter real and simulated data
real_data_filtered <- filter_genes_with_variance(real_data)
simulated_data_filtered <- filter_genes_with_variance(simulated_data)

# Function to calculate pairwise correlations for genes
calculate_gene_correlation <- function(data) {
  # Transpose the data so genes are in columns
  data_t <- t(as.matrix(data))
  # Calculate pairwise correlations between genes
  correlations <- cor(data_t, method = "pearson")
  # Flatten the correlation matrix and remove self-correlations
  correlations <- correlations[upper.tri(correlations)]
  return(correlations)
}

# Calculate correlations for filtered real and simulated data
real_correlations <- calculate_gene_correlation(real_data_filtered)
simulated_correlations <- calculate_gene_correlation(simulated_data_filtered)

# Create density plots
library(ggplot2)

real_density_plot <- ggplot(data.frame(correlation = real_correlations), aes(x = correlation)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Gene-wise Correlation Density (Real Data)", x = "Correlation", y = "Density") +
  theme_minimal()

simulated_density_plot <- ggplot(data.frame(correlation = simulated_correlations), aes(x = correlation)) +
  geom_density(fill = "red", alpha = 0.5) +
  labs(title = "Gene-wise Correlation Density (Simulated Data)", x = "Correlation", y = "Density") +
  theme_minimal()

# Print plots
print(real_density_plot)
print(simulated_density_plot)

# Calculate library size (total counts per sample)
library_size <- colSums(simulated_data)

# Calculate mean expression of genes
mean_expression <- rowMeans(simulated_data)

# Calculate variance expression of genes
variance_expression <- apply(simulated_data, 1, var)

# Create a data frame to combine the metrics
metrics <- data.frame(
  Metric = c(
    rep("Library Size", length(library_size)),
    rep("Mean Expression", length(mean_expression)),
    rep("Variance Expression", length(variance_expression))
  ),
  Value = c(library_size, mean_expression, variance_expression)
)

# Improved boxplot with log-scale y-axis
boxplot_plot <- ggplot(metrics, aes(x = Metric, y = Value, fill = Metric)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1, outlier.colour = "red") +
  scale_y_log10() +  # Apply log10 scale to handle large range of values
  labs(
    title = "Boxplots of Library Size, Mean Expression, and Variance Expression (Log Scale)",
    x = "Metric",
    y = "Log-scaled Value"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
  )

# Print the improved plot
print(boxplot_plot)