#!/usr/bin/env python3
"""
Comparative Analyzer
===================

Compare crystal packing metrics between multiple structures.
"""

import glob
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

try:
    import sklearn  # noqa: F401
    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

class ComparativeAnalyzer:
    """Analyzer for comparing multiple crystal structures."""
    
    def __init__(self):
        """Initialize the comparative analyzer."""
        pass
    
    def compare_structures(self, structure_results):
        """
        Compare multiple crystal structures.
        
        Parameters:
        -----------
        structure_results : list
            List of analysis results from individual structures
            
        Returns:
        --------
        dict : Comparative analysis results
        """
        if len(structure_results) < 2:
            return {'error': 'Need at least 2 structures for comparison'}
        
        try:
            results = {
                'structure_count': len(structure_results),
                'summary_stats': {},
                'correlations': {},
                'clustering': {},
                'pca_results': {},
                'outliers': []
            }
            
            # Extract metrics into DataFrame
            metrics_df = self._extract_metrics_dataframe(structure_results)
            
            if metrics_df.empty:
                return {'error': 'No comparable metrics found'}
            
            # Calculate summary statistics
            summary_stats = self._calculate_summary_statistics(metrics_df)
            results['summary_stats'] = summary_stats
            
            # Calculate correlations
            correlations = self._calculate_correlations(metrics_df)
            results['correlations'] = correlations
            
            # Perform clustering if sklearn available
            if SKLEARN_AVAILABLE:
                clustering_results = self._perform_clustering(metrics_df)
                results['clustering'] = clustering_results
                
                # PCA analysis
                pca_results = self._perform_pca(metrics_df)
                results['pca_results'] = pca_results
                
                # Outlier detection
                outliers = self._detect_outliers(metrics_df, structure_results)
                results['outliers'] = outliers
            
            return results
            
        except Exception as e:
            return {'error': f"Failed to compare structures: {str(e)}"}
    
    def _extract_metrics_dataframe(self, structure_results):
        """Extract numerical metrics into a pandas DataFrame."""
        metrics_list = []
        
        for result in structure_results:
            if 'error' in result:
                continue
                
            metrics = {}
            metrics['structure_id'] = result.get('structure_id', 'unknown')
            
            # Extract packing metrics
            packing = result.get('packing_metrics', {})
            if isinstance(packing, dict):
                metrics.update({
                    'matthews_coefficient': packing.get('matthews_coefficient', 0),
                    'solvent_content_percent': packing.get('solvent_content_percent', 0),
                    'packing_efficiency_percent': packing.get('packing_efficiency_percent', 0),
                    'molecular_weight': packing.get('molecular_weight', 0),
                    'unit_cell_volume': packing.get('unit_cell_volume', 0),
                    'total_atoms': packing.get('total_atoms', 0)
                })
            
            # Extract interface metrics
            interface = result.get('interface_analysis', {})
            if isinstance(interface, dict):
                summary = interface.get('summary', {})
                if isinstance(summary, dict):
                    metrics.update({
                        'total_interfaces': summary.get('total_interfaces', 0),
                        'total_buried_surface_area': summary.get('total_buried_surface_area', 0),
                        'average_buried_area_per_interface': summary.get('average_buried_area_per_interface', 0),
                        'average_interface_rmsd_ca': summary.get('average_interface_rmsd_ca', 0),
                    })
            
            # Extract contact metrics
            contact = result.get('contact_analysis', {})
            if isinstance(contact, dict):
                contact_summary = contact.get('contact_summary', {})
                if isinstance(contact_summary, dict):
                    asu_stats = contact_summary.get('asu_stats', {})
                    crystal_stats = contact_summary.get('crystal_stats', {})
                    
                    if isinstance(asu_stats, dict):
                        metrics.update({
                            'asu_contact_count': asu_stats.get('total_contacts', 0),
                            'average_contact_distance': asu_stats.get('average_distance', 0)
                        })
                    
                    if isinstance(crystal_stats, dict):
                        metrics.update({
                            'surface_atoms': crystal_stats.get('surface_atoms', 0),
                            'crystal_packing_efficiency': crystal_stats.get('packing_efficiency', 0)
                        })
            
            # Extract graph metrics
            graph = result.get('graph_analysis', {})
            if isinstance(graph, dict):
                residue_graph = graph.get('residue_graph', {})
                if isinstance(residue_graph, dict) and 'error' not in residue_graph:
                    metrics.update({
                        'residue_graph_density': residue_graph.get('density', 0),
                        'residue_graph_clustering': residue_graph.get('clustering_coefficient', 0),
                        'residue_graph_nodes': residue_graph.get('node_count', 0)
                    })
            
            if len(metrics) > 1:  # More than just structure_id
                metrics_list.append(metrics)
        
        if not metrics_list:
            return pd.DataFrame()
        
        df = pd.DataFrame(metrics_list)
        df = df.set_index('structure_id')
        
        # Remove non-numeric columns and handle missing values
        numeric_df = df.select_dtypes(include=[np.number])
        numeric_df = numeric_df.fillna(0)
        
        return numeric_df
    
    def _calculate_summary_statistics(self, df):
        """Calculate summary statistics for all metrics."""
        summary = {}
        
        for column in df.columns:
            values = df[column].values
            summary[column] = {
                'mean': float(np.mean(values)),
                'std': float(np.std(values)),
                'min': float(np.min(values)),
                'max': float(np.max(values)),
                'median': float(np.median(values))
            }
        
        return summary
    
    def _calculate_correlations(self, df):
        """Calculate pairwise correlations between metrics."""
        correlations = {}
        
        columns = df.columns.tolist()
        
        for i, col1 in enumerate(columns):
            for j, col2 in enumerate(columns):
                if i < j:  # Avoid duplicates
                    try:
                        corr = np.corrcoef(df[col1].values, df[col2].values)[0, 1]
                        if not np.isnan(corr):
                            correlations[f"{col1}_vs_{col2}"] = float(corr)
                    except Exception:
                        pass
        
        return correlations
    
    def _perform_clustering(self, df):
        """Perform clustering analysis."""
        if not SKLEARN_AVAILABLE or df.empty:
            return {'error': 'Sklearn not available or empty data'}
        
        try:
            from sklearn.preprocessing import StandardScaler
            from sklearn.cluster import KMeans
            # Standardize the data
            scaler = StandardScaler()
            scaled_data = scaler.fit_transform(df.values)
            
            # Determine optimal number of clusters
            n_structures = len(df)
            max_clusters = min(5, n_structures - 1)
            
            if max_clusters < 2:
                return {'error': 'Not enough structures for clustering'}
            
            # Perform K-means clustering
            best_k = 2
            best_inertia = float('inf')
            
            for k in range(2, max_clusters + 1):
                kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)  # type: ignore[arg-type]
                kmeans.fit_predict(scaled_data)
                inertia = getattr(kmeans, 'inertia_', None)
                if inertia is not None and inertia < best_inertia:
                    best_inertia = inertia
                    best_k = k
            
            # Final clustering with best k
            kmeans = KMeans(n_clusters=best_k, random_state=42, n_init=10)  # type: ignore[arg-type]
            cluster_labels = kmeans.fit_predict(scaled_data)
            
            # Assign clusters to structures
            clusters = {}
            for i, (structure_id, _) in enumerate(df.iterrows()):
                cluster_id = int(cluster_labels[i])
                if cluster_id not in clusters:
                    clusters[cluster_id] = []
                clusters[cluster_id].append(structure_id)
            
            return {
                'optimal_clusters': best_k,
                'cluster_assignments': clusters,
                'inertia': best_inertia
            }
            
        except Exception as e:
            return {'error': f"Clustering failed: {str(e)}"}
    
    def _perform_pca(self, df):
        """Perform Principal Component Analysis."""
        if not SKLEARN_AVAILABLE or df.empty:
            return {'error': 'Sklearn not available or empty data'}
        
        try:
            from sklearn.preprocessing import StandardScaler
            from sklearn.decomposition import PCA
            # Standardize the data
            scaler = StandardScaler()
            scaled_data = scaler.fit_transform(df.values)
            
            # Perform PCA
            n_components = min(len(df.columns), len(df))
            pca = PCA(n_components=n_components)
            pca_result = pca.fit_transform(scaled_data)
            
            # Calculate cumulative variance
            explained_variance = pca.explained_variance_ratio_
            cumulative_variance = np.cumsum(explained_variance)
            
            return {
                'explained_variance_ratio': explained_variance.tolist(),
                'cumulative_variance': cumulative_variance.tolist(),
                'n_components': n_components,
                'pc1_loadings': pca.components_[0].tolist() if len(pca.components_) > 0 else [],
                'pc2_loadings': pca.components_[1].tolist() if len(pca.components_) > 1 else []
            }
            
        except Exception as e:
            return {'error': f"PCA failed: {str(e)}"}
    
    def _detect_outliers(self, df, structure_results):
        """Detect outlier structures."""
        outliers = []
        
        try:
            # Simple outlier detection using z-score
            for column in df.columns:
                values = df[column].values
                if np.std(values) > 0:
                    z_scores = np.abs((values - np.mean(values)) / np.std(values))
                    
                    # Structures with z-score > 2.5 are considered outliers
                    outlier_indices = np.where(z_scores > 2.5)[0]
                    
                    for idx in outlier_indices:
                        structure_id = df.index[idx]
                        outliers.append({
                            'structure_id': structure_id,
                            'metric': column,
                            'value': float(values[idx]),
                            'z_score': float(z_scores[idx])
                        })
        
        except Exception:
            pass

        return outliers

def collect_structure_paths(inputs):
    """Expand inputs (files, directories, glob patterns) to a list of structure file paths."""
    paths = []
    for arg in inputs:
        arg = os.path.abspath(os.path.expanduser(arg))
        if os.path.isdir(arg):
            for ext in ('*.pdb', '*.cif', '*.ent'):
                paths.extend(glob.glob(os.path.join(arg, ext)))
        elif '*' in arg or '?' in arg:
            paths.extend(glob.glob(arg))
        elif os.path.isfile(arg):
            paths.append(arg)
    return sorted(set(paths))


def filter_paths_by_patterns(paths, patterns):
    """Return paths whose basename contains every pattern in `patterns`."""
    if not patterns:
        return list(paths)
    return [p for p in paths if all(pat in os.path.basename(p) for pat in patterns)]


def _run_comparison(calculator, analyzer, paths, out_stream):
    """Compute packing metrics for paths, compare, and write results to out_stream. Requires >= 2 paths."""
    all_results = []
    for pdb_file in paths:
        sid = Path(pdb_file).stem
        res = calculator.calculate_metrics(pdb_file)
        if 'error' in res:
            print(f"Warning: Skipping {pdb_file}: {res['error']}", file=sys.stderr)
            continue
        all_results.append({'structure_id': sid, 'pdb_file': pdb_file, 'packing_metrics': res})

    if len(all_results) < 2:
        print("Need at least 2 successfully analyzed structures to compare.", file=out_stream)
        return

    comparison = analyzer.compare_structures(all_results)
    if 'error' in comparison:
        print(f"Comparison error: {comparison['error']}", file=out_stream)
        return

    print("COMPARATIVE ANALYSIS RESULTS", file=out_stream)
    print("=" * 50, file=out_stream)
    print(f"Structures: {len(all_results)}", file=out_stream)
    summary_stats = comparison.get('summary_stats')
    if isinstance(summary_stats, dict):
        print("\nSummary statistics:", file=out_stream)
        for metric, stats in summary_stats.items():
            print(f"  {metric}: mean={stats.get('mean', 0):.3f} std={stats.get('std', 0):.3f}", file=out_stream)
    correlations = comparison.get('correlations')
    if isinstance(correlations, dict):
        print("\nCorrelations:", file=out_stream)
        for pair, val in correlations.items():
            print(f"  {pair}: {val:.3f}", file=out_stream)
    clustering = comparison.get('clustering', {})
    if isinstance(clustering, dict) and 'error' not in clustering:
        print("\nClustering:", file=out_stream)
        print(f"  Optimal clusters: {clustering.get('optimal_clusters', 0)}", file=out_stream)
        for cid, members in (clustering.get('cluster_assignments') or {}).items():
            print(f"  Cluster {cid}: {members}", file=out_stream)
    outliers = comparison.get('outliers', [])
    if outliers:
        print("\nOutliers:", file=out_stream)
        for o in outliers:
            if isinstance(o, dict):
                print(f"  {o.get('structure_id')} {o.get('metric')}: value={o.get('value')} z={o.get('z_score')}", file=out_stream)
    print("\nDone.", file=out_stream)


def main():
    """CLI: run packing metrics on each input, then compare. Optional multi-set filtering and dry-run."""
    import argparse

    try:
        from packing_metrics import PackingMetricsCalculator
    except ImportError:
        PackingMetricsCalculator = None

    parser = argparse.ArgumentParser(
        description="Compare crystal packing metrics between multiple structures.",
        epilog="""Examples:
  python comparative_analyzer.py model_01.pdb model_02.pdb
  python comparative_analyzer.py dir/
  python comparative_analyzer.py *.pdb -o comparison_results.txt
  python comparative_analyzer.py *.pdb --set set_a,set_b -o by_set.txt
  python comparative_analyzer.py *.pdb --sets set_a set_b set_c -o "output_{}.txt"
  python comparative_analyzer.py *.pdb --sets set_a set_b -o results.txt --dry-run
""",
    )
    parser.add_argument('input', nargs='+', help='PDB/CIF file(s), directory, or glob. Need at least 2 structures per set.')
    parser.add_argument(
        '--set', '-s', action='append', dest='sets', metavar='PATTERNS',
        help='Comma-separated patterns; file included only if basename contains ALL. Repeat for multiple sets.',
    )
    parser.add_argument(
        '--sets', dest='sets_multi', nargs='+', metavar='SET',
        help='Multiple set names in one go (one pattern per set). Each set needs >= 2 structures.',
    )
    parser.add_argument(
        '--output', '-o', metavar='FILE',
        help='Output file. Single file for all sets, or use "{}" for one file per set (e.g. output_{}.txt).',
    )
    parser.add_argument('--dry-run', action='store_true', help='Print which files would be processed per set and exit.')
    args = parser.parse_args()

    if not PackingMetricsCalculator:
        print("packing_metrics module required for standalone comparison.", file=sys.stderr)
        sys.exit(1)

    paths = collect_structure_paths(args.input)
    if not paths:
        print("No structure files found.", file=sys.stderr)
        sys.exit(1)

    set_list = []
    if getattr(args, 'sets_multi', None):
        for name in args.sets_multi:
            patterns = [name.strip()] if name.strip() else []
            if patterns:
                label = patterns[0]
                filtered = filter_paths_by_patterns(paths, patterns)
                set_list.append((label, patterns, filtered))
    if args.sets:
        for s in args.sets:
            patterns = [p.strip() for p in s.split(',') if p.strip()]
            if patterns:
                label = '_'.join(patterns)
                filtered = filter_paths_by_patterns(paths, patterns)
                set_list.append((label, patterns, filtered))
    if not set_list:
        set_list = [("all", [], list(paths))]

    if args.dry_run:
        print("Dry run: no analysis will be performed.\n", file=sys.stderr)
        single_file = args.output and '{}' not in args.output
        for label, patterns, filtered in set_list:
            dest = args.output if single_file else (args.output.replace('{}', label) if args.output and '{}' in args.output else None)
            dest = dest or "stdout"
            note = " (need >= 2 for comparison)" if len(filtered) < 2 else ""
            print(f"Set {label!r} (patterns: {patterns}): {len(filtered)} file(s) -> {dest}{note}", file=sys.stderr)
            for p in filtered:
                print(f"  {os.path.basename(p)}", file=sys.stderr)
            print(file=sys.stderr)
        return

    single_output_file = args.output and '{}' not in args.output
    out = sys.stdout
    if single_output_file:
        out = open(args.output, 'w')
        print(f"Writing all set comparisons to {args.output}", file=sys.stderr)

    calculator = PackingMetricsCalculator()
    analyzer = ComparativeAnalyzer()
    for label, patterns, filtered in set_list:
        if len(filtered) < 2:
            print(f"Set {label!r} has {len(filtered)} file(s); need at least 2 to compare, skipping.", file=sys.stderr)
            continue

        current_out = sys.stdout
        if single_output_file:
            current_out = out
        elif args.output:
            out_path = args.output.replace('{}', label)
            current_out = open(out_path, 'w')
            print(f"Writing comparison for set {label!r} ({len(filtered)} file(s)) to {out_path}", file=sys.stderr)
        elif len(set_list) > 1:
            print(f"\n=== Set {label!r} ({len(filtered)} file(s)) ===\n", file=sys.stderr)

        try:
            if len(set_list) > 1 and (single_output_file or not args.output):
                print(f"\n{'='*50}\nSet {label!r} (patterns: {patterns})\n{'='*50}", file=current_out)
            _run_comparison(calculator, analyzer, filtered, current_out)
        finally:
            if args.output and not single_output_file:
                current_out.close()

    if single_output_file:
        out.close()


if __name__ == "__main__":
    main() 