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


setwd("D:/Users/Nazia Khan/Documents/Fall 2021/Home_Work/Project")


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

# Save the simulated data
write.csv(simulated_counts_splatter, "simulated_bulk_splatter.csv", row.names = TRUE)

print("Splatter simulation complete. Data saved to simulated_bulk_splatter.csv.")
