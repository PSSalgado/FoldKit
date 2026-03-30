#!/usr/bin/env python3
"""
Quick Start Guide for Crystal Packing Analysis
==============================================

This script demonstrates basic usage of the crystal packing analysis pipeline.
"""

import sys
from pathlib import Path

def main():
    print("Crystal Packing Analysis Pipeline - Quick Start")
    print("=" * 50)
    
    print("\nThis pipeline analyzes crystal lattice packing differences.")
    print("Key features:")
    print("• Packing density metrics (Matthews coefficient, solvent content)")
    print("• Interface analysis (buried surface area, contacts)")
    print("• Crystal contact analysis")
    print("• Solvent channel analysis")
    print("• Graph-theoretical analysis")
    print("• Comparative statistics and visualization")
    
    print("\n" + "="*50)
    print("USAGE EXAMPLES:")
    print("="*50)
    
    print("\n1. Analyze a single structure:")
    print("   python crystal_packing_analyzer.py --input model_01.pdb")
    
    print("\n2. Compare multiple structures:")
    print("   python crystal_packing_analyzer.py --input *.pdb --compare")
    
    print("\n3. Specify output directory:")
    print("   python crystal_packing_analyzer.py --input model_01.pdb --output analysis_output")
    
    print("\n4. Use individual modules:")
    print("   python packing_metrics.py model_01.pdb")
    print("   python interface_analyzer.py model_01.pdb")
    
    print("\n5. Generate R visualizations:")
    print("   cd comparison_plots/")
    print("   Rscript master_visualization.R")
    
    print("\n" + "="*50)
    print("INSTALLATION:")
    print("="*50)
    print("\n1. Install dependencies:")
    print("   pip install -r requirements.txt")
    
    print("\n2. Required packages:")
    print("   • numpy, pandas, scipy")
    print("   • biopython")
    print("   • scikit-learn")
    print("   • networkx")
    print("   • R (for visualization)")
    print("   • R packages: ggplot2, viridis, pheatmap, etc.")
    
    print("\n" + "="*50)
    print("FILES CREATED:")
    print("="*50)
    
    files = [
        "crystal_packing_analyzer.py - Main analysis pipeline",
        "packing_metrics.py - Basic packing calculations",
        "interface_analyzer.py - Interface analysis",
        "contact_analyzer.py - Crystal contact analysis",
        "channel_analyzer.py - Solvent channel analysis",
        "graph_analyzer.py - Graph-theoretical analysis",
        "comparative_analyzer.py - Multi-structure comparison",
        "visualization.py - Plotting and visualization",
        "requirements.txt - Required Python packages",
        "README.md - Comprehensive documentation"
    ]
    
    for file_desc in files:
        filename = file_desc.split(" - ")[0]
        if Path(filename).exists():
            print(f"✓ {file_desc}")
        else:
            print(f"✗ {file_desc}")
    
    print("\n" + "="*50)
    print("NEXT STEPS:")
    print("="*50)
    print("\n1. Install dependencies: pip install -r requirements.txt")
    print("2. Install R from https://www.r-project.org/")
    print("3. Test with your PDB files")
    print("4. Generate R visualizations for publication-quality plots")
    print("5. Check README.md for detailed documentation")
    print("6. Customize analysis parameters as needed")
    
    print("\nFor help: python crystal_packing_analyzer.py --help")

if __name__ == "__main__":
    main() 