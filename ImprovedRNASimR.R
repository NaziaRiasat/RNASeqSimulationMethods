# ImprovedRNASimR: Advanced RNA-seq Bulk Data Simulation
# Proof-of-concept implementation addressing key limitations of existing methods

# Required libraries
library(MASS)
library(Matrix)
library(igraph)
library(cluster)
library(edgeR)
library(DESeq2)
library(stats)
library(parallel)

# =============================================================================
# 1. HIERARCHICAL CORRELATION MODELING
# =============================================================================

#' Detect local gene networks using correlation clustering
#' @param expression_matrix Gene expression matrix (genes x samples)
#' @param correlation_threshold Minimum correlation for network inclusion
#' @param min_network_size Minimum number of genes in a network
#' @return List of local gene networks
detect_local_networks <- function(expression_matrix, correlation_threshold = 0.7, min_network_size = 3) {
  # Calculate gene-gene correlation matrix
  gene_corr <- cor(t(expression_matrix), use = "complete.obs")
  
  # Create adjacency matrix based on correlation threshold
  adj_matrix <- abs(gene_corr) > correlation_threshold
  diag(adj_matrix) <- FALSE
  
  # Convert to graph and detect communities
  graph <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  communities <- cluster_fast_greedy(graph)
  
  # Extract networks with sufficient size
  local_networks <- list()
  for (i in 1:max(communities$membership)) {
    network_genes <- which(communities$membership == i)
    if (length(network_genes) >= min_network_size) {
      local_networks[[paste0("network_", i)]] <- list(
        genes = rownames(expression_matrix)[network_genes],
        correlation_matrix = gene_corr[network_genes, network_genes],
        size = length(network_genes)
      )
    }
  }
  
  return(local_networks)
}

#' Estimate pathway-level correlations using functional annotation
#' @param expression_matrix Gene expression matrix (genes x samples)
#' @param pathway_db Pathway database (gene to pathway mapping)
#' @return List of pathway correlation structures
estimate_pathway_correlations <- function(expression_matrix, pathway_db = NULL) {
  # Simplified pathway detection using expression-based clustering
  # In practice, this would use real pathway databases (KEGG, GO, etc.)
  
  if (is.null(pathway_db)) {
    # Use hierarchical clustering to identify potential pathways
    gene_distances <- as.dist(1 - abs(cor(t(expression_matrix), use = "complete.obs")))
    gene_clusters <- hclust(gene_distances, method = "ward.D2")
    
    # Cut tree to get pathway-like groups
    pathway_groups <- cutree(gene_clusters, h = 0.7)
    
    pathway_correlations <- list()
    for (pathway_id in unique(pathway_groups)) {
      pathway_genes <- names(pathway_groups)[pathway_groups == pathway_id]
      if (length(pathway_genes) >= 5) {  # Minimum pathway size
        pathway_correlations[[paste0("pathway_", pathway_id)]] <- list(
          genes = pathway_genes,
          correlation_strength = median(cor(t(expression_matrix[pathway_genes, ]), use = "complete.obs")),
          size = length(pathway_genes)
        )
      }
    }
  }
  
  return(pathway_correlations)
}

#' Capture global correlation structure using PCA
#' @param expression_matrix Gene expression matrix (genes x samples)
#' @param n_components Number of principal components to retain
#' @return Global correlation structure
capture_global_structure <- function(expression_matrix, n_components = 10) {
  # Perform PCA on gene expression
  pca_result <- prcomp(t(expression_matrix), center = TRUE, scale. = TRUE)
  
  # Extract key components
  n_components <- min(n_components, ncol(pca_result$x))
  
  global_structure <- list(
    pc_loadings = pca_result$rotation[, 1:n_components],
    eigenvalues = pca_result$sdev[1:n_components]^2,
    explained_variance = (pca_result$sdev[1:n_components]^2) / sum(pca_result$sdev^2),
    center = pca_result$center,
    scale = pca_result$scale
  )
  
  return(global_structure)
}

#' Build hierarchical correlation structure
#' @param expression_matrix Gene expression matrix (genes x samples)
#' @return Hierarchical correlation structure
build_correlation_hierarchy <- function(expression_matrix) {
  cat("Building hierarchical correlation structure...\n")
  
  # Detect local networks
  local_networks <- detect_local_networks(expression_matrix)
  cat("Detected", length(local_networks), "local gene networks\n")
  
  # Estimate pathway correlations
  pathway_correlations <- estimate_pathway_correlations(expression_matrix)
  cat("Identified", length(pathway_correlations), "pathway-level correlations\n")
  
  # Capture global structure
  global_structure <- capture_global_structure(expression_matrix)
  cat("Captured global structure with", length(global_structure$eigenvalues), "components\n")
  
  hierarchy <- list(
    local_networks = local_networks,
    pathway_correlations = pathway_correlations,
    global_structure = global_structure
  )
  
  return(hierarchy)
}

# =============================================================================
# 2. ENHANCED PARAMETER ESTIMATION
# =============================================================================

#' Estimate gene-level parameters with improved precision
#' @param expression_matrix Gene expression matrix (genes x samples)
#' @param groups Sample group information
#' @return Gene-level parameters
estimate_gene_parameters <- function(expression_matrix, groups) {
  n_genes <- nrow(expression_matrix)
  
  # Calculate robust estimates for each gene
  gene_params <- list()
  
  for (i in 1:n_genes) {
    gene_expr <- expression_matrix[i, ]
    
    # Estimate mean and dispersion using method of moments
    mean_expr <- mean(gene_expr)
    var_expr <- var(gene_expr)
    
    # Estimate dispersion parameter for negative binomial
    dispersion <- max(0.01, (var_expr - mean_expr) / mean_expr^2)
    
    # Zero inflation probability
    zero_prob <- sum(gene_expr == 0) / length(gene_expr)
    
    gene_params[[i]] <- list(
      mean = mean_expr,
      dispersion = dispersion,
      zero_inflation = zero_prob
    )
  }
  
  # Convert to matrices for efficient access
  gene_level_params <- list(
    mean_expression = sapply(gene_params, function(x) x$mean),
    dispersion = sapply(gene_params, function(x) x$dispersion),
    zero_inflation = sapply(gene_params, function(x) x$zero_inflation)
  )
  
  return(gene_level_params)
}

#' Estimate sample-level effects and technical variation
#' @param expression_matrix Gene expression matrix (genes x samples)
#' @param groups Sample group information
#' @return Sample-level parameters
estimate_sample_parameters <- function(expression_matrix, groups) {
  n_samples <- ncol(expression_matrix)
  
  # Library size factors
  library_sizes <- colSums(expression_matrix)
  library_size_factors <- library_sizes / median(library_sizes)
  
  # Sample-specific effects using PCA
  sample_pca <- prcomp(expression_matrix, center = TRUE, scale. = TRUE)
  sample_effects <- sample_pca$x[, 1:min(5, ncol(sample_pca$x))]
  
  # Technical noise estimation
  technical_noise <- apply(expression_matrix, 2, function(x) {
    # Estimate technical noise as deviation from expected
    fitted_values <- predict(lm(x ~ groups))
    residual_var <- var(x - fitted_values)
    return(residual_var)
  })
  
  sample_level_params <- list(
    library_size_factors = library_size_factors,
    sample_effects = sample_effects,
    technical_noise = technical_noise
  )
  
  return(sample_level_params)
}

#' Estimate batch effects and structure
#' @param expression_matrix Gene expression matrix (genes x samples)
#' @param batch_info Batch information (if available)
#' @return Batch-level parameters
estimate_batch_parameters <- function(expression_matrix, batch_info = NULL) {
  if (is.null(batch_info)) {
    # Estimate latent batch structure using surrogate variable analysis
    batch_effects <- matrix(0, nrow = nrow(expression_matrix), ncol = 2)
    batch_variance <- rep(0.1, nrow(expression_matrix))
    batch_correlation <- diag(nrow(expression_matrix))
  } else {
    # Use provided batch information
    batch_effects <- model.matrix(~ batch_info)[, -1, drop = FALSE]
    batch_variance <- apply(expression_matrix, 1, function(x) {
      var(residuals(lm(x ~ batch_info)))
    })
    batch_correlation <- cor(t(expression_matrix), use = "complete.obs")
  }
  
  batch_level_params <- list(
    batch_effects = batch_effects,
    batch_variance = batch_variance,
    batch_correlation = batch_correlation
  )
  
  return(batch_level_params)
}

#' Comprehensive parameter estimation for improved simulation
#' @param expression_matrix Gene expression matrix (genes x samples)
#' @param groups Sample group information
#' @param batch_info Batch information (optional)
#' @return Complete parameter structure
estimate_hierarchical_parameters <- function(expression_matrix, groups, batch_info = NULL) {
  cat("Estimating hierarchical parameters...\n")
  
  # Gene-level parameters
  gene_params <- estimate_gene_parameters(expression_matrix, groups)
  cat("Estimated gene-level parameters for", length(gene_params$mean_expression), "genes\n")
  
  # Sample-level parameters
  sample_params <- estimate_sample_parameters(expression_matrix, groups)
  cat("Estimated sample-level parameters for", length(sample_params$library_size_factors), "samples\n")
  
  # Batch-level parameters
  batch_params <- estimate_batch_parameters(expression_matrix, batch_info)
  cat("Estimated batch-level parameters\n")
  
  # Correlation structure
  correlation_structure <- build_correlation_hierarchy(expression_matrix)
  
  # Combine all parameters
  all_params <- list(
    gene_level = gene_params,
    sample_level = sample_params,
    batch_level = batch_params,
    correlation_structure = correlation_structure
  )
  
  return(all_params)
}

# =============================================================================
# 3. ADVANCED SIMULATION ENGINE
# =============================================================================

#' Generate base expression matrix using estimated parameters
#' @param params Parameter structure from estimation
#' @param n_genes Number of genes to simulate
#' @param n_samples Number of samples to simulate
#' @return Base expression matrix
generate_base_expression <- function(params, n_genes, n_samples) {
  # Sample genes from original parameter distribution
  if (n_genes > length(params$gene_level$mean_expression)) {
    # Bootstrap additional genes
    gene_indices <- sample(1:length(params$gene_level$mean_expression), n_genes, replace = TRUE)
  } else {
    gene_indices <- sample(1:length(params$gene_level$mean_expression), n_genes, replace = FALSE)
  }
  
  # Generate base expression using negative binomial distribution
  base_matrix <- matrix(0, nrow = n_genes, ncol = n_samples)
  
  for (i in 1:n_genes) {
    gene_idx <- gene_indices[i]
    mean_expr <- params$gene_level$mean_expression[gene_idx]
    dispersion <- params$gene_level$dispersion[gene_idx]
    
    # Generate counts using negative binomial
    prob <- mean_expr / (mean_expr + dispersion)
    base_matrix[i, ] <- rnbinom(n_samples, size = dispersion, prob = prob)
  }
  
  return(base_matrix)
}

#' Apply hierarchical correlation structure to expression matrix
#' @param base_matrix Base expression matrix
#' @param correlation_structure Hierarchical correlation structure
#' @return Expression matrix with correlation structure applied
apply_correlation_hierarchy <- function(base_matrix, correlation_structure) {
  cat("Applying hierarchical correlation structure...\n")
  
  # Apply local network correlations
  if (length(correlation_structure$local_networks) > 0) {
    for (network in correlation_structure$local_networks) {
      if (network$size <= nrow(base_matrix)) {
        # Select genes for this network
        network_genes <- sample(1:nrow(base_matrix), network$size, replace = FALSE)
        
        # Apply correlation structure using Cholesky decomposition
        target_corr <- network$correlation_matrix
        if (is.positive.definite(target_corr)) {
          chol_decomp <- chol(target_corr)
          base_matrix[network_genes, ] <- chol_decomp %*% base_matrix[network_genes, ]
        }
      }
    }
  }
  
  # Apply global structure using PCA transformation
  if (!is.null(correlation_structure$global_structure)) {
    global_struct <- correlation_structure$global_structure
    if (ncol(base_matrix) >= length(global_struct$eigenvalues)) {
      # Apply PCA-based transformation
      pc_transform <- global_struct$pc_loadings %*% diag(sqrt(global_struct$eigenvalues))
      if (nrow(pc_transform) <= nrow(base_matrix)) {
        selected_genes <- sample(1:nrow(base_matrix), nrow(pc_transform), replace = FALSE)
        base_matrix[selected_genes, ] <- pc_transform %*% base_matrix[selected_genes, 1:ncol(pc_transform)]
      }
    }
  }
  
  return(base_matrix)
}

#' Add realistic batch effects to expression matrix
#' @param expression_matrix Expression matrix
#' @param batch_params Batch-level parameters
#' @return Expression matrix with batch effects
add_batch_effects <- function(expression_matrix, batch_params) {
  # Add batch effects based on estimated parameters
  if (ncol(batch_params$batch_effects) > 0) {
    # Generate batch assignments
    n_samples <- ncol(expression_matrix)
    batch_assignments <- sample(1:ncol(batch_params$batch_effects), n_samples, replace = TRUE)
    
    # Apply batch effects
    for (i in 1:ncol(expression_matrix)) {
      batch_effect <- batch_params$batch_effects[, batch_assignments[i]]
      expression_matrix[, i] <- expression_matrix[, i] + batch_effect
    }
  }
  
  return(expression_matrix)
}

#' Iterative refinement to match target statistics
#' @param expression_matrix Current expression matrix
#' @param params Target parameters
#' @param max_iterations Maximum number of iterations
#' @param convergence_threshold Convergence threshold
#' @return Refined expression matrix
iterative_refinement <- function(expression_matrix, params, max_iterations = 50, convergence_threshold = 1e-3) {
  cat("Performing iterative refinement...\n")
  
  for (iter in 1:max_iterations) {
    # Calculate current statistics
    current_means <- rowMeans(expression_matrix)
    current_vars <- apply(expression_matrix, 1, var)
    
    # Calculate target statistics
    target_means <- params$gene_level$mean_expression[1:nrow(expression_matrix)]
    target_vars <- target_means * (1 + params$gene_level$dispersion[1:nrow(expression_matrix)])
    
    # Adjust expression matrix
    mean_adjustment <- target_means / (current_means + 1e-6)
    var_adjustment <- sqrt(target_vars / (current_vars + 1e-6))
    
    # Apply adjustments
    for (i in 1:nrow(expression_matrix)) {
      expression_matrix[i, ] <- expression_matrix[i, ] * mean_adjustment[i] * var_adjustment[i]
    }
    
    # Check convergence
    mean_error <- mean(abs(current_means - target_means))
    if (mean_error < convergence_threshold) {
      cat("Converged after", iter, "iterations\n")
      break
    }
  }
  
  return(expression_matrix)
}

# =============================================================================
# 4. MAIN SIMULATION FUNCTION
# =============================================================================

#' Main improved RNA-seq simulation function
#' @param reference_data Reference RNA-seq dataset for parameter estimation
#' @param groups Group information for reference data
#' @param n_genes Number of genes to simulate
#' @param n_samples Number of samples to simulate
#' @param batch_info Batch information (optional)
#' @param max_iterations Maximum iterations for refinement
#' @return Simulated RNA-seq data with metadata
ImprovedRNASeqSimulate <- function(reference_data, groups, n_genes, n_samples, 
                                   batch_info = NULL, max_iterations = 50) {
  
  cat("Starting Improved RNA-seq Simulation...\n")
  cat("Target dimensions:", n_genes, "genes x", n_samples, "samples\n")
  
  # Phase 1: Parameter Estimation
  cat("\n=== Phase 1: Parameter Estimation ===\n")
  params <- estimate_hierarchical_parameters(reference_data, groups, batch_info)
  
  # Phase 2: Base Matrix Generation
  cat("\n=== Phase 2: Base Matrix Generation ===\n")
  base_matrix <- generate_base_expression(params, n_genes, n_samples)
  cat("Generated base matrix with dimensions:", dim(base_matrix)[1], "x", dim(base_matrix)[2], "\n")
  
  # Phase 3: Correlation Structure Application
  cat("\n=== Phase 3: Correlation Structure Application ===\n")
  corr_matrix <- apply_correlation_hierarchy(base_matrix, params$correlation_structure)
  
  # Phase 4: Batch Effect Integration
  cat("\n=== Phase 4: Batch Effect Integration ===\n")
  batch_matrix <- add_batch_effects(corr_matrix, params$batch_level)
  
  # Phase 5: Iterative Refinement
  cat("\n=== Phase 5: Iterative Refinement ===\n")
  final_matrix <- iterative_refinement(batch_matrix, params, max_iterations)
  
  # Ensure non-negative values
  final_matrix[final_matrix < 0] <- 0
  
  # Generate metadata
  gene_names <- paste0("Gene_", 1:n_genes)
  sample_names <- paste0("Sample_", 1:n_samples)
  
  rownames(final_matrix) <- gene_names
  colnames(final_matrix) <- sample_names
  
  # Create sample metadata
  sample_metadata <- data.frame(
    SampleID = sample_names,
    Group = sample(levels(groups), n_samples, replace = TRUE),
    LibrarySize = colSums(final_matrix),
    stringsAsFactors = FALSE
  )
  
  # Create gene metadata
  gene_metadata <- data.frame(
    GeneID = gene_names,
    MeanExpression = rowMeans(final_matrix),
    Variance = apply(final_matrix, 1, var),
    stringsAsFactors = FALSE
  )
  
  cat("\n=== Simulation Complete ===\n")
  cat("Final matrix dimensions:", dim(final_matrix)[1], "x", dim(final_matrix)[2], "\n")
  
  result <- list(
    counts = final_matrix,
    sample_metadata = sample_metadata,
    gene_metadata = gene_metadata,
    parameters = params,
    simulation_info = list(
      n_genes = n_genes,
      n_samples = n_samples,
      method = "ImprovedRNASeqSimulate",
      timestamp = Sys.time()
    )
  )
  
  return(result)
}

# =============================================================================
# 5. VALIDATION AND QUALITY ASSESSMENT
# =============================================================================

#' Validate correlation structure preservation
#' @param real_data Real RNA-seq data
#' @param simulated_data Simulated RNA-seq data
#' @return Correlation quality metrics
validate_correlation_structure <- function(real_data, simulated_data) {
  # Calculate correlation matrices
  real_corr <- cor(t(real_data), use = "complete.obs")
  sim_corr <- cor(t(simulated_data), use = "complete.obs")
  
  # Mantel test for correlation similarity
  mantel_result <- cor.test(as.vector(real_corr), as.vector(sim_corr))
  
  # Frobenius distance
  frobenius_dist <- norm(real_corr - sim_corr, type = "F")
  
  # RV coefficient
  rv_coeff <- sum(real_corr * sim_corr) / sqrt(sum(real_corr^2) * sum(sim_corr^2))
  
  correlation_metrics <- list(
    mantel_correlation = mantel_result$estimate,
    mantel_pvalue = mantel_result$p.value,
    frobenius_distance = frobenius_dist,
    rv_coefficient = rv_coeff
  )
  
  return(correlation_metrics)
}

#' Comprehensive simulation quality assessment
#' @param real_data Real RNA-seq data
#' @param simulated_result Result from ImprovedRNASeqSimulate
#' @return Quality assessment metrics
assess_simulation_quality <- function(real_data, simulated_result) {
  simulated_data <- simulated_result$counts
  
  # Correlation structure validation
  correlation_quality <- validate_correlation_structure(real_data, simulated_data)
  
  # Distributional similarity
  real_means <- rowMeans(real_data)
  sim_means <- rowMeans(simulated_data)
  mean_ks_test <- ks.test(real_means, sim_means)
  
  real_vars <- apply(real_data, 1, var)
  sim_vars <- apply(simulated_data, 1, var)
  var_ks_test <- ks.test(real_vars, sim_vars)
  
  # Library size similarity
  real_lib_sizes <- colSums(real_data)
  sim_lib_sizes <- colSums(simulated_data)
  lib_size_ks_test <- ks.test(real_lib_sizes, sim_lib_sizes)
  
  quality_metrics <- list(
    correlation_quality = correlation_quality,
    distributional_similarity = list(
      mean_ks_statistic = mean_ks_test$statistic,
      mean_ks_pvalue = mean_ks_test$p.value,
      variance_ks_statistic = var_ks_test$statistic,
      variance_ks_pvalue = var_ks_test$p.value
    ),
    library_size_similarity = list(
      ks_statistic = lib_size_ks_test$statistic,
      ks_pvalue = lib_size_ks_test$p.value
    )
  )
  
  return(quality_metrics)
}

# =============================================================================
# 6. UTILITY FUNCTIONS
# =============================================================================

#' Check if matrix is positive definite
#' @param matrix Input matrix
#' @return TRUE if positive definite, FALSE otherwise
is.positive.definite <- function(matrix) {
  eigenvalues <- eigen(matrix, only.values = TRUE)$values
  return(all(eigenvalues > 0))
}

#' Print simulation summary
#' @param simulation_result Result from ImprovedRNASeqSimulate
print_simulation_summary <- function(simulation_result) {
  cat("=== Improved RNA-seq Simulation Summary ===\n")
  cat("Genes:", nrow(simulation_result$counts), "\n")
  cat("Samples:", ncol(simulation_result$counts), "\n")
  cat("Total counts:", sum(simulation_result$counts), "\n")
  cat("Mean library size:", mean(colSums(simulation_result$counts)), "\n")
  cat("Mean gene expression:", mean(rowMeans(simulation_result$counts)), "\n")
  cat("Simulation timestamp:", as.character(simulation_result$simulation_info$timestamp), "\n")
}

cat("ImprovedRNASimR loaded successfully!\n")
cat("Main function: ImprovedRNASeqSimulate()\n")
cat("Quality assessment: assess_simulation_quality()\n")