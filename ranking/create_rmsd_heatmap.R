#!/usr/bin/env Rscript

# Install required packages if not already installed
# Removed tidyverse dependency, using only essential packages
packages <- c("viridis", "pheatmap", "RColorBrewer")
install_missing_packages <- function() {
  cat("Checking and installing required packages...\n")
  for (package in packages) {
    if (!requireNamespace(package, quietly = TRUE)) {
      cat("Installing package:", package, "\n")
      tryCatch({
        install.packages(package, repos = "https://cloud.r-project.org", dependencies = TRUE)
      }, error = function(e) {
        cat("Failed to install", package, ". Error:", conditionMessage(e), "\n")
      })
    }
  }
  
  # Verify all packages are available
  missing_packages <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    cat("Warning: The following packages are still missing:", paste(missing_packages, collapse=", "), "\n")
  } else {
    cat("All required packages are installed successfully.\n")
  }
}

# Run package installation
install_missing_packages()

# Load required libraries with error handling
load_packages <- function() {
  # Load required packages
  if (requireNamespace("viridis", quietly = TRUE)) {
    library(viridis)
  } else {
    cat("Warning: viridis package is missing. Heatmaps may not render correctly.\n")
  }
  
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    library(pheatmap)
  } else {
    cat("Error: pheatmap package is missing. Cannot continue.\n")
    stop("Required package pheatmap is not available")
  }
  
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    library(RColorBrewer)
  } else {
    cat("Warning: RColorBrewer package is missing. Will use default colour palette.\n")
  }
}

# Load packages
load_packages()

# Function to create a more gradual colour palette
# palette: RColorBrewer name (e.g. RdYlBu, RdYlGn, YlOrRd) or "viridis"/"plasma"
create_gradient_palette <- function(n = 100, palette = "RdYlBu") {
  if (palette %in% c("viridis", "plasma", "inferno", "magma", "cividis") && requireNamespace("viridis", quietly = TRUE)) {
    return(viridis(n, option = palette))
  }
  if (requireNamespace("RColorBrewer", quietly = TRUE) && palette %in% rownames(brewer.pal.info)) {
    n_brewer <- min(11, max(3, brewer.pal.info[palette, "maxcolors"]))
    if (palette %in% c("RdYlBu", "RdYlGn", "RdBu", "PiYG", "BrBG")) {
      colors <- colorRampPalette(rev(brewer.pal(n_brewer, palette)))(n)
    } else {
      colors <- colorRampPalette(brewer.pal(n_brewer, palette))(n)
    }
    return(colors)
  }
  if (requireNamespace("viridis", quietly = TRUE)) return(viridis(n))
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    return(colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(n))
  }
  stop("Need viridis or RColorBrewer for colour palettes")
}

# Function to process a single CSV file
process_rmsd_file <- function(csv_file, color_palette_name = "RdYlBu") {
    # Read the CSV file
    data <- read.csv(csv_file, row.names = 1, check.names = FALSE)
    
    # Convert '-' to NA
    data[data == '-'] <- NA
    
    # Convert remaining values to numeric
    data <- as.data.frame(apply(data, 2, as.numeric))
    rownames(data) <- rownames(data)
    
    # Get subdomain from filename
    subdomain <- gsub("rmsd_table_", "", gsub(".csv", "", basename(csv_file)))
    
    # Create output filename
    output_file <- file.path(dirname(csv_file), paste0("rmsd_heatmap_", subdomain, ".pdf"))
    
    # Create colour palette
    color_palette <- create_gradient_palette(100, palette = color_palette_name)
    
    # Create heatmap
    pheatmap(data,
             color = color_palette,
             na_col = "white",
             display_numbers = TRUE,
             number_format = "%.2f",
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             main = paste("RMSD Heatmap -", subdomain),
             filename = output_file,
             width = 10,
             height = 10,
             fontsize = 8,
             number_color = "black")
    
    print(paste("Created heatmap:", output_file))
}

# Function to process the combined CSV file
process_combined_file <- function(combined_file, color_palette_name = "RdYlBu") {
    # Check if the file exists
    if (!file.exists(combined_file)) {
        cat("Combined file not found:", combined_file, "\n")
        return(FALSE)
    }
    
    # Read the CSV file
    data <- read.csv(combined_file, check.names = FALSE)
    
    # Check that the required columns are present
    if (!"Subdomain" %in% colnames(data) || !"Model" %in% colnames(data)) {
        cat("Combined file does not have the expected format\n")
        return(FALSE)
    }
    
    # Define the preferred group/subdomain order (customise per CSV group labels)
    subdomain_order <- c("sd", "D1", "D2", "ID", "ad")
    
    # Extract the protein name prefix from model names
    # Model column is used as the label; Subdomain column identifies the group
    # And Subdomain column contains the subdomain (ad, sd, etc.)
    
    # Add a column for protein names (extracted from Model)
    data$Protein <- data$Model
    
    # Get all subdomains
    all_subdomains <- unique(data$Subdomain)
    
    # Order subdomains based on the preferred order
    # First check whether any of the preferred subdomains are present
    subdomains <- c()
    
    # First add subdomains in the specified order
    for (sd in subdomain_order) {
        if (sd %in% all_subdomains) {
            subdomains <- c(subdomains, sd)
        }
    }
    
    # Then add any remaining subdomains not in the specified order
    for (sd in all_subdomains) {
        if (!(sd %in% subdomains)) {
            subdomains <- c(subdomains, sd)
        }
    }
    
    # Print the order of subdomains
    cat("Using subdomain order:", paste(subdomains, collapse=", "), "\n")
    
    # Get all unique proteins
    all_proteins <- unique(data$Protein)
    cat("Found", length(all_proteins), "unique proteins\n")
    
    # Create output filename
    output_file <- file.path(dirname(combined_file), "combined_rmsd_heatmap.pdf")
    
    # Create colour palette
    color_palette <- create_gradient_palette(100, palette = color_palette_name)
    
    # Create a list to store protein names for each subdomain
    subdomain_proteins <- list()
    
    # Get unique proteins for each subdomain
    for (subdomain in subdomains) {
        subdomain_data <- data[data$Subdomain == subdomain, ]
        # Sort proteins alphabetically
        subdomain_proteins[[subdomain]] <- sort(unique(subdomain_data$Protein))
    }
    
    # Create model names in correct order by subdomain
    model_names <- c()
    for (subdomain in subdomains) {
        # Each subdomain contains alphabetically sorted proteins
        for (protein in subdomain_proteins[[subdomain]]) {
            model_names <- c(model_names, paste(protein, subdomain, sep = "_"))
        }
    }
    
    # Create a matrix to hold all RMSD values
    rmsd_matrix <- matrix(NA, nrow = length(model_names), ncol = length(model_names))
    rownames(rmsd_matrix) <- model_names
    colnames(rmsd_matrix) <- model_names
    
    # Verify the order of subdomains in the matrix row names (for debugging)
    matrix_subdomains <- sub(".*_", "", rownames(rmsd_matrix))
    cat("Matrix subdomain order:", paste(unique(matrix_subdomains), collapse=", "), "\n")
    
    # Fill the matrix with RMSD values
    for (i in 1:nrow(data)) {
        # Get the current row from the data
        row_subdomain <- data$Subdomain[i]
        row_protein <- data$Protein[i]
        row_name <- paste(row_protein, row_subdomain, sep = "_")
        
        if (row_name %in% rownames(rmsd_matrix)) {
            # Get all column proteins in the same subdomain
            for (col_protein in all_proteins) {
                col_name <- paste(col_protein, row_subdomain, sep = "_")
                
                if (col_name %in% colnames(rmsd_matrix) && col_protein != row_protein) {
                    # Get RMSD value
                    rmsd_val <- data[i, col_protein]
                    if (!is.na(rmsd_val) && rmsd_val != "-" && rmsd_val != "") {
                        rmsd_matrix[row_name, col_name] <- as.numeric(rmsd_val)
                    }
                }
            }
        }
    }
    
    # Create subdomain annotation
    annotation_row <- data.frame(
        Subdomain = factor(sub(".*_", "", rownames(rmsd_matrix)), levels = subdomains)
    )
    rownames(annotation_row) <- rownames(rmsd_matrix)
    
    annotation_col <- data.frame(
        Subdomain = factor(sub(".*_", "", colnames(rmsd_matrix)), levels = subdomains)
    )
    rownames(annotation_col) <- colnames(rmsd_matrix)
    
    # Create subdomain colours with matching order
    subdomain_colors <- colorRampPalette(brewer.pal(8, "Dark2"))(length(subdomains))
    names(subdomain_colors) <- subdomains
    annotation_colors <- list(Subdomain = subdomain_colors)
    
    # Use only protein names for display
    # These are the parts before the "_" in the model_names
    row_labels <- sub("_.*", "", rownames(rmsd_matrix))
    col_labels <- sub("_.*", "", colnames(rmsd_matrix))
    
    # Create heatmap
    pheatmap(rmsd_matrix,
             color = color_palette,
             na_col = "white",
             display_numbers = FALSE,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             main = "Combined RMSD Heatmap",
             filename = output_file,
             width = 20,
             height = 20,
             fontsize = 8,
             annotation_row = annotation_row,
             annotation_col = annotation_col,
             annotation_colors = annotation_colors,
             labels_row = row_labels,
             labels_col = col_labels)
    
    # Create a second version with clustering
    output_file_clustered <- file.path(dirname(combined_file), "combined_rmsd_heatmap_clustered.pdf")
    
    # Try to create clustered heatmap, but handle errors gracefully
    tryCatch({
        # Calculate distance matrix with NA handling
        # First create a copy of the matrix and replace NAs with a high value
        clustering_matrix <- rmsd_matrix
        # Find the maximum value and replace NAs with a value higher than max
        max_val <- max(clustering_matrix, na.rm = TRUE)
        clustering_matrix[is.na(clustering_matrix)] <- max_val * 1.5
        
        # For self-comparisons (diagonal), use 0
        diag(clustering_matrix) <- 0
        
        # Now create the clustered heatmap with the modified matrix
        # Cluster within subdomains by using gaps
        
        # Determine where the subdomain boundaries are
        subdomain_sizes <- sapply(subdomains, function(sd) length(subdomain_proteins[[sd]]))
        gap_positions <- cumsum(subdomain_sizes)
        gap_positions <- gap_positions[-length(gap_positions)]  # Remove the last one
        
        pheatmap(rmsd_matrix,  # Use original matrix for display
                 color = color_palette,
                 na_col = "white",
                 display_numbers = FALSE,
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 clustering_distance_rows = dist(clustering_matrix),  # Pre-calculated distance
                 clustering_distance_cols = dist(clustering_matrix),  # Pre-calculated distance
                 clustering_method = "complete",
                 gaps_row = gap_positions,
                 gaps_col = gap_positions,
                 main = "Combined RMSD Heatmap (Clustered)",
                 filename = output_file_clustered,
                 width = 25,
                 height = 25,
                 fontsize = 8,
                 annotation_row = annotation_row,
                 annotation_col = annotation_col,
                 annotation_colors = annotation_colors,
                 labels_row = row_labels,
                 labels_col = col_labels)
        
        cat("Created clustered combined heatmap:", output_file_clustered, "\n")
    }, error = function(e) {
        cat("Could not create clustered heatmap due to error:", conditionMessage(e), "\n")
        cat("This is usually due to missing values or NA patterns in the data.\n")
        cat("Only the non-clustered version will be available.\n")
    })
    
    return(TRUE)
}

# RStudio standalone mode
# If running in RStudio, ask for directory or use working directory
color_palette_name <- "RdYlBu"
if (interactive()) {
  # Use file dialog to select directory
  cat("Select the directory containing RMSD CSV files\n")
  base_dir <- readline(prompt = "Enter directory path (press Enter to use current working directory): ")
  
  if (base_dir == "") {
    base_dir <- getwd()
    cat("Using current working directory:", base_dir, "\n")
  }
} else {
  # Get command line arguments when running as script
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) {
    stop("Usage: Rscript create_rmsd_heatmap.R <base_directory> [palette]\n  palette: RdYlBu (default), RdYlGn, YlOrRd, viridis, plasma, etc.")
  }
  base_dir <- args[1]
  color_palette_name <- if (length(args) >= 2) args[2] else "RdYlBu"
  cat("Using colour palette:", color_palette_name, "\n")
}

# Find all rmsd_table_*.csv files
csv_files <- list.files(base_dir, 
                       pattern = "rmsd_table_.*\\.csv$",
                       recursive = TRUE,
                       full.names = TRUE)

if (length(csv_files) == 0) {
  cat("No RMSD CSV files found in", base_dir, "\n")
} else {
  cat("Found", length(csv_files), "RMSD CSV files to process\n")
  
  # Process each CSV file
  for (csv_file in csv_files) {
    cat("Processing:", csv_file, "\n")
    process_rmsd_file(csv_file, color_palette_name = color_palette_name)
  }
}

# Check for combined table
combined_file <- file.path(base_dir, "combined_rmsd_table.csv")
if (file.exists(combined_file)) {
  cat("Found combined RMSD table, creating combined heatmap...\n")
  process_combined_file(combined_file, color_palette_name = color_palette_name)
} else {
  cat("No combined RMSD table found at", combined_file, "\n")
  cat("Run rmsd_to_csv.py --scan-dir <base> first to generate the combined table.\n")
} 