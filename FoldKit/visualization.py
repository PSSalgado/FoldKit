#!/usr/bin/env python3
"""
R-based Visualization Module
===========================

Generate R scripts for creating plots and visualizations for crystal packing analysis results.
This approach integrates with existing R-based workflows and provides publication-quality plots.
"""

import numpy as np
import pandas as pd
from pathlib import Path
import json

class PackingVisualizer:
    """Visualizer for crystal packing analysis results using R."""
    
    def __init__(self):
        """Initialize the visualizer."""
        self.r_template_header = '''
# Crystal Packing Analysis Visualization
# Generated automatically by the crystal packing analysis pipeline
# Date: {date}

# Load required libraries
if (!require("ggplot2")) install.packages("ggplot2", dependencies=TRUE)
if (!require("viridis")) install.packages("viridis", dependencies=TRUE)
if (!require("pheatmap")) install.packages("pheatmap", dependencies=TRUE)
if (!require("corrplot")) install.packages("corrplot", dependencies=TRUE)
if (!require("RColorBrewer")) install.packages("RColorBrewer", dependencies=TRUE)
if (!require("reshape2")) install.packages("reshape2", dependencies=TRUE)
if (!require("gridExtra")) install.packages("gridExtra", dependencies=TRUE)

library(ggplot2)
library(viridis)
library(pheatmap)
library(corrplot)
library(RColorBrewer)
library(reshape2)
library(gridExtra)

# Set theme for consistent plotting
theme_set(theme_minimal() + theme(
  plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "bottom"
))
'''
    
    def create_comparison_plots(self, comparison_results, output_dir):
        """
        Create R scripts for comparison plots of multiple structures.
        
        Parameters:
        -----------
        comparison_results : dict
            Results from comparative analysis
        output_dir : Path or str
            Directory to save R scripts and data
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            # Save data as CSV files for R
            self._save_data_for_r(comparison_results, output_dir)
            
            # Generate R scripts
            self._generate_summary_stats_script(comparison_results, output_dir)
            self._generate_correlations_script(comparison_results, output_dir)
            self._generate_pca_script(comparison_results, output_dir)
            self._generate_clustering_script(comparison_results, output_dir)
            
            # Create master script to run all analyses
            self._create_master_script(output_dir)
            
            print(f"R scripts and data saved to {output_dir}")
            print("To generate plots, run: Rscript master_visualization.R")
            
        except Exception as e:
            print(f"Error creating R scripts: {e}")
    
    def _save_data_for_r(self, comparison_results, output_dir):
        """Save analysis data as CSV files for R processing."""
        
        # Summary statistics
        summary_stats = comparison_results.get('summary_stats', {})
        if summary_stats:
            summary_df = pd.DataFrame(summary_stats).T
            summary_df.to_csv(output_dir / 'summary_statistics.csv')
        
        # Correlations
        correlations = comparison_results.get('correlations', {})
        if correlations:
            corr_data = []
            for pair, correlation in correlations.items():
                corr_data.append({'Metric_Pair': pair, 'Correlation': correlation})
            corr_df = pd.DataFrame(corr_data)
            corr_df.to_csv(output_dir / 'correlations.csv', index=False)
        
        # PCA results
        pca_results = comparison_results.get('pca_results', {})
        if 'explained_variance_ratio' in pca_results:
            pca_df = pd.DataFrame({
                'Component': range(1, len(pca_results['explained_variance_ratio']) + 1),
                'Explained_Variance': pca_results['explained_variance_ratio'],
                'Cumulative_Variance': pca_results.get('cumulative_variance', [])
            })
            pca_df.to_csv(output_dir / 'pca_results.csv', index=False)
        
        # Clustering results
        clustering = comparison_results.get('clustering', {})
        if 'cluster_assignments' in clustering:
            cluster_data = []
            for cluster_id, structures in clustering['cluster_assignments'].items():
                for structure in structures:
                    cluster_data.append({'Structure': structure, 'Cluster': f'Cluster {cluster_id}'})
            
            if cluster_data:
                cluster_df = pd.DataFrame(cluster_data)
                cluster_df.to_csv(output_dir / 'clustering_results.csv', index=False)
    
    def _generate_summary_stats_script(self, comparison_results, output_dir):
        """Generate R script for summary statistics plots."""
        
        summary_stats = comparison_results.get('summary_stats', {})
        if not summary_stats:
            return
        
        script_content = f'''
{self.r_template_header.format(date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"))}

# Summary Statistics Visualization
cat("Generating summary statistics plots...\\n")

# Read data
summary_data <- read.csv("summary_statistics.csv", row.names = 1)

# Create plots
pdf("summary_statistics.pdf", width = 15, height = 12)

# Prepare data for plotting
metrics <- rownames(summary_data)
means <- summary_data$mean
stds <- summary_data$std
cvs <- ifelse(means != 0, stds / means, 0)
ranges <- summary_data$max - summary_data$min

# Create a multi-panel plot
par(mfrow = c(2, 2), mar = c(10, 4, 3, 2))

# Mean values
barplot(means, names.arg = metrics, las = 2, 
        main = "Mean Values", ylab = "Value",
        col = viridis(length(metrics), alpha = 0.7))

# Standard deviations
barplot(stds, names.arg = metrics, las = 2,
        main = "Standard Deviations", ylab = "Standard Deviation",
        col = plasma(length(metrics), alpha = 0.7))

# Coefficient of variation
barplot(cvs, names.arg = metrics, las = 2,
        main = "Coefficient of Variation", ylab = "CV (std/mean)",
        col = inferno(length(metrics), alpha = 0.7))

# Ranges
barplot(ranges, names.arg = metrics, las = 2,
        main = "Ranges (Max - Min)", ylab = "Range",
        col = cividis(length(metrics), alpha = 0.7))

dev.off()

# Alternative ggplot2 version
library(reshape2)
plot_data <- data.frame(
  Metric = metrics,
  Mean = means,
  Std = stds,
  CV = cvs,
  Range = ranges
)

# Melt data for ggplot
plot_data_long <- melt(plot_data, id.vars = "Metric")

# Create faceted plot
p <- ggplot(plot_data_long, aes(x = reorder(Metric, value), y = value, fill = variable)) +
  geom_col(alpha = 0.7) +
  facet_wrap(~variable, scales = "free") +
  coord_flip() +
  scale_fill_viridis_d() +
  labs(title = "Crystal Packing Metrics Summary",
       x = "Metrics", y = "Value") +
  theme(strip.text = element_text(size = 12, face = "bold"))

ggsave("summary_statistics_ggplot.pdf", p, width = 15, height = 12)

cat("Summary statistics plots saved to summary_statistics.pdf and summary_statistics_ggplot.pdf\\n")
'''
        
        with open(output_dir / 'summary_statistics.R', 'w') as f:
            f.write(script_content)
    
    def _generate_correlations_script(self, comparison_results, output_dir):
        """Generate R script for correlation plots."""
        
        correlations = comparison_results.get('correlations', {})
        if not correlations:
            return
        
        script_content = f'''
{self.r_template_header.format(date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"))}

# Correlation Analysis Visualization
cat("Generating correlation plots...\\n")

# Read correlation data
corr_data <- read.csv("correlations.csv")

# Sort by absolute correlation value
corr_data$abs_corr <- abs(corr_data$Correlation)
corr_data <- corr_data[order(-corr_data$abs_corr), ]

# Take top 20 correlations for readability
top_corr <- head(corr_data, 20)

# Create correlation plot
pdf("correlations.pdf", width = 12, height = 10)

# Horizontal bar plot
par(mar = c(5, 15, 4, 2))
colors <- ifelse(top_corr$Correlation >= 0, "steelblue", "red")

barplot(top_corr$Correlation, names.arg = gsub("_vs_", "\\nvs\\n", top_corr$Metric_Pair),
        horiz = TRUE, las = 1, col = colors, alpha = 0.7,
        main = "Top Pairwise Correlations Between Metrics",
        xlab = "Correlation Coefficient")

# Add reference line at zero
abline(v = 0, col = "black", lty = 2, alpha = 0.5)

# Add value labels
text(top_corr$Correlation + ifelse(top_corr$Correlation >= 0, 0.02, -0.02),
     1:nrow(top_corr), sprintf("%.3f", top_corr$Correlation),
     pos = ifelse(top_corr$Correlation >= 0, 4, 2), cex = 0.8)

dev.off()

# Alternative ggplot2 version
p <- ggplot(top_corr, aes(x = reorder(Metric_Pair, Correlation), y = Correlation)) +
  geom_col(aes(fill = Correlation > 0), alpha = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "red")) +
  scale_x_discrete(labels = function(x) gsub("_vs_", "\\nvs\\n", x)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Correlation)),
            hjust = ifelse(top_corr$Correlation >= 0, -0.1, 1.1), size = 3) +
  labs(title = "Top Pairwise Correlations Between Metrics",
       x = "Metric Pairs", y = "Correlation Coefficient") +
  theme(legend.position = "none")

ggsave("correlations_ggplot.pdf", p, width = 12, height = 10)

cat("Correlation plots saved to correlations.pdf and correlations_ggplot.pdf\\n")
'''
        
        with open(output_dir / 'correlations.R', 'w') as f:
            f.write(script_content)
    
    def _generate_pca_script(self, comparison_results, output_dir):
        """Generate R script for PCA plots."""
        
        pca_results = comparison_results.get('pca_results', {})
        if 'explained_variance_ratio' in pca_results:
            
            script_content = f'''
{self.r_template_header.format(date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"))}

# PCA Analysis Visualization
cat("Generating PCA plots...\\n")

# Read PCA data
pca_data <- read.csv("pca_results.csv")

# Create PCA plots
pdf("pca_analysis.pdf", width = 15, height = 6)

par(mfrow = c(1, 2))

# Explained variance plot
barplot(pca_data$Explained_Variance, names.arg = paste("PC", pca_data$Component),
        main = "Explained Variance by Component",
        ylab = "Explained Variance Ratio",
        xlab = "Principal Component",
        col = viridis(nrow(pca_data), alpha = 0.7))
grid()

# Cumulative variance plot
plot(pca_data$Component, pca_data$Cumulative_Variance, type = "o",
     pch = 16, col = "steelblue", lwd = 2,
     main = "Cumulative Explained Variance",
     xlab = "Number of Components",
     ylab = "Cumulative Explained Variance",
     ylim = c(0, 1))
abline(h = 0.95, col = "red", lty = 2, lwd = 2)
text(max(pca_data$Component) * 0.7, 0.97, "95% threshold", col = "red")
grid()

dev.off()

# ggplot2 version
p1 <- ggplot(pca_data, aes(x = factor(Component), y = Explained_Variance)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  labs(title = "Explained Variance by Component",
       x = "Principal Component", y = "Explained Variance Ratio") +
  scale_x_discrete(labels = paste("PC", pca_data$Component))

p2 <- ggplot(pca_data, aes(x = Component, y = Cumulative_Variance)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "steelblue", size = 2) +
  geom_hline(yintercept = 0.95, color = "red", linetype = "dashed", size = 1) +
  annotate("text", x = max(pca_data$Component) * 0.7, y = 0.97, 
           label = "95% threshold", color = "red") +
  labs(title = "Cumulative Explained Variance",
       x = "Number of Components", y = "Cumulative Explained Variance") +
  ylim(0, 1)

combined_plot <- grid.arrange(p1, p2, ncol = 2)
ggsave("pca_analysis_ggplot.pdf", combined_plot, width = 15, height = 6)

cat("PCA plots saved to pca_analysis.pdf and pca_analysis_ggplot.pdf\\n")
'''
            
            with open(output_dir / 'pca_analysis.R', 'w') as f:
                f.write(script_content)
    
    def _generate_clustering_script(self, comparison_results, output_dir):
        """Generate R script for clustering plots."""
        
        clustering = comparison_results.get('clustering', {})
        if 'cluster_assignments' in clustering:
            
            script_content = f'''
{self.r_template_header.format(date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"))}

# Clustering Analysis Visualization
cat("Generating clustering plots...\\n")

# Read clustering data
cluster_data <- read.csv("clustering_results.csv")

# Count structures per cluster
cluster_counts <- table(cluster_data$Cluster)

# Create clustering plots
pdf("clustering_results.pdf", width = 10, height = 8)

# Pie chart of cluster sizes
pie_colors <- brewer.pal(length(cluster_counts), "Set3")
pie(cluster_counts, labels = paste(names(cluster_counts), "\\n(", cluster_counts, ")"),
    col = pie_colors, main = paste("Structure Distribution Across", length(cluster_counts), "Clusters"))

dev.off()

# ggplot2 version with better layout
cluster_summary <- data.frame(
  Cluster = names(cluster_counts),
  Count = as.numeric(cluster_counts),
  Percentage = as.numeric(cluster_counts) / sum(cluster_counts) * 100
)

p <- ggplot(cluster_summary, aes(x = "", y = Count, fill = Cluster)) +
  geom_col(width = 1) +
  coord_polar("y", start = 0) +
  scale_fill_brewer(type = "qual", palette = "Set3") +
  labs(title = paste("Structure Distribution Across", length(cluster_counts), "Clusters")) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave("clustering_results_ggplot.pdf", p, width = 10, height = 8)

# Create a detailed table
cat("\\nCluster Assignments:\\n")
for (cluster in unique(cluster_data$Cluster)) {{
  structures <- cluster_data$Structure[cluster_data$Cluster == cluster]
  cat(paste(cluster, ":", paste(structures, collapse = ", "), "\\n"))
}}

# Save cluster assignments to text file
sink("cluster_assignments.txt")
cat("Cluster Assignments\\n")
cat("==================\\n\\n")
for (cluster in unique(cluster_data$Cluster)) {{
  structures <- cluster_data$Structure[cluster_data$Cluster == cluster]
  cat(paste(cluster, ":", paste(structures, collapse = ", "), "\\n"))
}}
sink()

cat("Clustering plots saved to clustering_results.pdf and clustering_results_ggplot.pdf\\n")
cat("Cluster assignments saved to cluster_assignments.txt\\n")
'''
            
            with open(output_dir / 'clustering_analysis.R', 'w') as f:
                f.write(script_content)
    
    def _create_master_script(self, output_dir):
        """Create a master R script to run all visualizations."""
        
        master_script = f'''
# Master Visualization Script for Crystal Packing Analysis
# Generated automatically by the crystal packing analysis pipeline
# Date: {pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S")}

cat("===========================================\\n")
cat("Crystal Packing Analysis - Visualization\\n")
cat("===========================================\\n\\n")

# Set working directory to script location
script_dir <- dirname(sys.frame(1)$ofile)
if (!is.null(script_dir)) {{
  setwd(script_dir)
}}

# Check for required R scripts and run them
scripts_to_run <- c(
  "summary_statistics.R",
  "correlations.R", 
  "pca_analysis.R",
  "clustering_analysis.R"
)

for (script in scripts_to_run) {{
  if (file.exists(script)) {{
    cat(paste("Running", script, "...\\n"))
    tryCatch({{
      source(script)
      cat(paste("✓", script, "completed successfully\\n\\n"))
    }}, error = function(e) {{
      cat(paste("✗ Error in", script, ":", e$message, "\\n\\n"))
    }})
  }} else {{
    cat(paste("⚠ Script", script, "not found, skipping...\\n\\n"))
  }}
}}

cat("===========================================\\n")
cat("Visualization complete!\\n")
cat("Check the following files for results:\\n")
cat("- summary_statistics.pdf\\n")
cat("- correlations.pdf\\n") 
cat("- pca_analysis.pdf\\n")
cat("- clustering_results.pdf\\n")
cat("- cluster_assignments.txt\\n")
cat("===========================================\\n")
'''
        
        with open(output_dir / 'master_visualization.R', 'w') as f:
            f.write(master_script)
    
    def plot_single_structure_summary(self, results, output_file):
        """Create R script for a single structure summary plot."""
        
        output_path = Path(output_file)
        output_dir = output_path.parent
        script_name = output_path.stem + "_plot.R"
        
        # Extract data for plotting
        structure_id = results.get('structure_id', 'Unknown')
        
        # Prepare data
        plot_data = {}
        
        # Packing metrics
        packing = results.get('packing_metrics', {})
        if isinstance(packing, dict):
            plot_data['packing'] = {
                'Matthews Coefficient': packing.get('matthews_coefficient', 0),
                'Solvent Content (%)': packing.get('solvent_content_percent', 0),
                'Packing Efficiency (%)': packing.get('packing_efficiency_percent', 0)
            }
        
        # Interface metrics
        interface = results.get('interface_analysis', {})
        if isinstance(interface, dict):
            summary = interface.get('summary', {})
            if isinstance(summary, dict):
                plot_data['interface'] = {
                    'Interfaces': summary.get('total_interfaces', 0),
                    'Buried Area (×100 Å²)': summary.get('total_buried_surface_area', 0) / 100
                }
        
        # Save data as JSON for R script
        data_file = output_dir / (output_path.stem + "_data.json")
        with open(data_file, 'w') as f:
            json.dump({'structure_id': structure_id, 'data': plot_data}, f, indent=2)
        
        # Generate R script
        script_content = f'''
{self.r_template_header.format(date=pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"))}

# Single Structure Summary Visualization
cat("Generating single structure summary plot...\\n")

# Load data
library(jsonlite)
data <- fromJSON("{data_file.name}")
structure_id <- data$structure_id

# Create summary plot
pdf("{output_file}", width = 15, height = 12)

par(mfrow = c(2, 2), mar = c(8, 4, 3, 2))

# Packing metrics
if (!is.null(data$data$packing)) {{
  packing_data <- unlist(data$data$packing)
  barplot(packing_data, las = 2, col = c("skyblue", "lightgreen", "lightcoral"),
          main = "Basic Packing Metrics", ylab = "Value")
}}

# Interface metrics  
if (!is.null(data$data$interface)) {{
  interface_data <- unlist(data$data$interface)
  barplot(interface_data, las = 2, col = c("orange", "purple"),
          main = "Interface Analysis", ylab = "Count / Area")
}}

# Add overall title
mtext(paste("Crystal Packing Analysis:", structure_id), 
      outer = TRUE, cex = 1.5, font = 2, line = -2)

dev.off()

cat("Single structure plot saved to {output_file}\\n")
'''
        
        with open(output_dir / script_name, 'w') as f:
            f.write(script_content)
        
        print(f"R script for single structure plot saved: {output_dir / script_name}")
        print(f"To generate plot, run: Rscript {script_name}")

def main():
    """Example usage of PackingVisualizer."""
    print("PackingVisualizer now generates R scripts for publication-quality plots.")
    print("The generated R scripts integrate with your existing R-based workflow.")
    print("Run the scripts using: Rscript master_visualization.R")

if __name__ == "__main__":
    main() 