#!/usr/bin/env python3
"""
Quick Start Guide for Crystal Packing Analysis
==============================================

This script demonstrates basic usage of the crystal packing analysis pipeline.
"""

import os
import sys
from pathlib import Path

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from utils.cli_log import setup_log_from_argv

_argv_no_log, _ = setup_log_from_argv(
    script_path=__file__,
    argv=sys.argv[1:],
    inputs=[],
    pattern=None,
)
sys.argv = [sys.argv[0]] + _argv_no_log

def main():
    print("Crystal Packing Analysis Pipeline - Quick Start")
    print("=" * 50)
    
    print("\nThis pipeline analyses crystal lattice packing differences.")
    print("Key features:")
    print("• Packing density metrics (Matthews coefficient, solvent content)")
    print("• Interface analysis (buried surface area, contacts)")
    print("• Crystal contact analysis")
    print("• Optional batch JSON (--compare)")
    
    print("\n" + "="*50)
    print("USAGE EXAMPLES:")
    print("="*50)
    
    print("\n1. Analyse a single structure:")
    print("   python metrics/crystal_packing_analyser.py --input model_01.pdb")
    
    print("\n2. Combine outputs from multiple structures:")
    print("   python metrics/crystal_packing_analyser.py --input *.pdb --compare")
    
    print("\n3. Specify output directory:")
    print("   python metrics/crystal_packing_analyser.py --input model_01.pdb --output analysis_output")
    
    print("\n4. Use individual modules:")
    print("   python metrics/packing_metrics.py model_01.pdb")
    print("   python metrics/interface_analyser_asu_charge.py model_01.pdb")
    
    print("\n" + "="*50)
    print("INSTALLATION:")
    print("="*50)
    print("\n1. Install dependencies:")
    print("   pip install -r requirements.txt")
    
    print("\n2. Required packages:")
    print("   • numpy, pandas, scipy")
    print("   • biopython")
    print("   • R (optional; only for documented ranking/*.R scripts)")
    
    print("\n" + "="*50)
    print("FILES CREATED:")
    print("="*50)
    
    files = [
        "crystal_packing_analyser.py - Main analysis pipeline",
        "packing_metrics.py - Basic packing calculations",
        "interface_analyser_asu_charge.py - Interface analysis (ASU, charge-tag metrics)",
        "interface_analyser_asu_ec.py - Interface analysis (ASU, electrostatic complementarity, McCoy)",
        "interface_analyser_lattice_charge.py - Interface analysis (lattice, charge-tag metrics)",
        "interface_analyser_lattice_ec.py - Interface analysis (lattice, electrostatic complementarity, McCoy)",
        "contact_analyser.py - Crystal contact analysis",
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
    print("2. Test with representative structure files")
    print("3. Check README.md for detailed documentation")
    print("4. Customise analysis parameters as needed")
    
    print("\nFor help: python metrics/crystal_packing_analyser.py --help")

if __name__ == "__main__":
    main() 