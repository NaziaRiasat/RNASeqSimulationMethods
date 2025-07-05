# Installation Script for Improved RNA-seq Simulation Method
# This script installs all required packages and dependencies

cat("Installing Required Packages for Improved RNA-seq Simulation Method\n")
cat("====================================================================\n")

# Function to install packages with error handling
install_package <- function(package_name, source = "CRAN") {
  if (source == "CRAN") {
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      cat("Installing", package_name, "from CRAN...\n")
      install.packages(package_name, dependencies = TRUE)
    } else {
      cat(package_name, "already installed\n")
    }
  } else if (source == "Bioconductor") {
    if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
      cat("Installing", package_name, "from Bioconductor...\n")
      BiocManager::install(package_name)
    } else {
      cat(package_name, "already installed\n")
    }
  }
}

# Install BiocManager first
if (!require("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager")
}

# Core statistical packages
cat("\n--- Installing Core Statistical Packages ---\n")
core_packages <- c("MASS", "Matrix", "stats", "parallel")
for (pkg in core_packages) {
  install_package(pkg, "CRAN")
}

# Network analysis packages
cat("\n--- Installing Network Analysis Packages ---\n")
network_packages <- c("igraph", "cluster")
for (pkg in network_packages) {
  install_package(pkg, "CRAN")
}

# Visualization packages
cat("\n--- Installing Visualization Packages ---\n")
viz_packages <- c("ggplot2", "gridExtra", "reshape2", "corrplot")
for (pkg in viz_packages) {
  install_package(pkg, "CRAN")
}

# RNA-seq specific packages
cat("\n--- Installing RNA-seq Specific Packages ---\n")
rnaseq_packages <- c("edgeR", "DESeq2", "SPsimSeq")
for (pkg in rnaseq_packages) {
  install_package(pkg, "Bioconductor")
}

# Additional utility packages
cat("\n--- Installing Additional Utility Packages ---\n")
utility_packages <- c("devtools", "testthat", "roxygen2", "data.table", "tibble", "dplyr")
for (pkg in utility_packages) {
  install_package(pkg, "CRAN")
}

# Optional packages for advanced features
cat("\n--- Installing Optional Packages ---\n")
optional_packages <- c("doParallel", "foreach", "pryr", "microbenchmark")
for (pkg in optional_packages) {
  tryCatch({
    install_package(pkg, "CRAN")
  }, error = function(e) {
    cat("Warning: Could not install optional package", pkg, "\n")
  })
}

# Check if all essential packages are loaded
cat("\n--- Verifying Installation ---\n")
essential_packages <- c("MASS", "Matrix", "igraph", "cluster", "ggplot2", 
                       "edgeR", "DESeq2", "SPsimSeq")

all_installed <- TRUE
for (pkg in essential_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("ERROR: Package", pkg, "failed to install\n")
    all_installed <- FALSE
  } else {
    cat("✓", pkg, "successfully installed\n")
  }
}

if (all_installed) {
  cat("\n🎉 All essential packages installed successfully!\n")
  cat("You can now use the Improved RNA-seq Simulation Method.\n")
  cat("\nNext steps:\n")
  cat("1. Load the simulation functions: source('ImprovedRNASimR.R')\n")
  cat("2. Follow the setup guide: see 'setup_guide.md'\n")
  cat("3. Run the demonstration: source('demonstration_script.R')\n")
} else {
  cat("\n❌ Some packages failed to install.\n")
  cat("Please check the error messages above and install missing packages manually.\n")
}

# System information
cat("\n--- System Information ---\n")
cat("R Version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("Available Memory:", round(as.numeric(gsub("[^0-9]", "", system("free -h | grep '^Mem' | awk '{print $2}'", intern = TRUE))), 1), "GB\n")
cat("Number of CPU cores:", parallel::detectCores(), "\n")

# Performance recommendations
cat("\n--- Performance Recommendations ---\n")
if (parallel::detectCores() >= 4) {
  cat("✓ You have", parallel::detectCores(), "CPU cores - good for parallel processing\n")
} else {
  cat("⚠ You have", parallel::detectCores(), "CPU cores - consider using smaller datasets\n")
}

# Check available memory (Linux/Mac)
if (Sys.info()["sysname"] %in% c("Linux", "Darwin")) {
  tryCatch({
    mem_info <- system("free -h 2>/dev/null || vm_stat 2>/dev/null", intern = TRUE)
    if (length(mem_info) > 0) {
      cat("✓ System memory information available\n")
    }
  }, error = function(e) {
    cat("ℹ Could not determine system memory\n")
  })
}

cat("\n--- Installation Complete ---\n")
cat("Date:", Sys.time(), "\n")
cat("For support, see the documentation files in this directory.\n")