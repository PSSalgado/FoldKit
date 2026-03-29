#!/usr/bin/env python3
"""
Graph Analyzer
=============

Graph-theoretical analysis of crystal packing using network representations.
"""

import glob
import os
import numpy as np
from pathlib import Path

try:
    from Bio.PDB.PDBParser import PDBParser
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False

class GraphAnalyzer:
    """Graph-theoretical analyzer for crystal packing."""
    
    def __init__(self, contact_threshold=5.0):
        """
        Initialize the graph analyzer.
        
        Parameters:
        -----------
        contact_threshold : float
            Distance threshold for considering residues connected (Å)
        """
        if BIOPYTHON_AVAILABLE:
            self.parser = PDBParser(QUIET=True)
        else:
            self.parser = None
            
        self.contact_threshold = contact_threshold
    
    def analyze_packing_graph(self, pdb_file):
        """
        Analyze crystal packing using graph theory.
        
        Parameters:
        -----------
        pdb_file : str
            Path to PDB file
            
        Returns:
        --------
        dict : Graph analysis results
        """
        if not BIOPYTHON_AVAILABLE or self.parser is None:
            return {'error': 'BioPython not available for graph analysis'}
        
        try:
            structure = self.parser.get_structure('crystal', pdb_file)
            
            results = {
                'residue_graph': {},
                'chain_graph': {},
                'network_metrics': {}
            }
            
            # Build residue contact graph
            residue_graph = self._build_residue_graph(structure)
            results['residue_graph'] = self._analyze_graph(residue_graph, 'residue')
            
            # Build chain interaction graph
            chain_graph = self._build_chain_graph(structure)
            results['chain_graph'] = self._analyze_graph(chain_graph, 'chain')
            
            # Calculate network metrics
            network_metrics = self._calculate_network_metrics(residue_graph, chain_graph)
            results['network_metrics'] = network_metrics
            
            return results
            
        except Exception as e:
            return {'error': f"Failed to analyze graph: {str(e)}"}
    
    def _build_residue_graph(self, structure):
        """Build a graph of residue contacts."""
        if not NETWORKX_AVAILABLE:
            return None
        
        G = nx.Graph()
        
        # Get all residues
        residues = list(structure.get_residues())
        
        # Add nodes (residues)
        for residue in residues:
            res_id = f"{residue.parent.id}_{residue.resname}{residue.id[1]}"
            G.add_node(res_id, 
                      chain=residue.parent.id,
                      resname=residue.resname,
                      resnum=residue.id[1])
        
        # Add edges (contacts)
        for i, res1 in enumerate(residues):
            for j, res2 in enumerate(residues):
                if i >= j:
                    continue
                
                # Calculate minimum distance between residues
                min_dist = self._min_residue_distance(res1, res2)
                
                if min_dist <= self.contact_threshold:
                    res1_id = f"{res1.parent.id}_{res1.resname}{res1.id[1]}"
                    res2_id = f"{res2.parent.id}_{res2.resname}{res2.id[1]}"
                    
                    G.add_edge(res1_id, res2_id, 
                              distance=min_dist,
                              inter_chain=(res1.parent.id != res2.parent.id))
        
        return G
    
    def _build_chain_graph(self, structure):
        """Build a graph of chain interactions."""
        if not NETWORKX_AVAILABLE:
            return None
        
        G = nx.Graph()
        
        chains = list(structure.get_chains())
        
        # Add nodes (chains)
        for chain in chains:
            G.add_node(chain.id, residue_count=len(list(chain.get_residues())))
        
        # Add edges (chain interactions)
        for i, chain1 in enumerate(chains):
            for j, chain2 in enumerate(chains):
                if i >= j:
                    continue
                
                # Count contacts between chains
                contact_count = self._count_chain_contacts(chain1, chain2)
                
                if contact_count > 0:
                    G.add_edge(chain1.id, chain2.id, 
                              contact_count=contact_count)
        
        return G
    
    def _min_residue_distance(self, res1, res2):
        """Calculate minimum distance between two residues."""
        min_dist = float('inf')
        
        for atom1 in res1.get_atoms():
            for atom2 in res2.get_atoms():
                dist = float(np.linalg.norm(atom1.coord - atom2.coord))
                min_dist = min(min_dist, dist)
        
        return min_dist
    
    def _count_chain_contacts(self, chain1, chain2):
        """Count contacts between two chains."""
        contact_count = 0
        
        for res1 in chain1.get_residues():
            for res2 in chain2.get_residues():
                min_dist = self._min_residue_distance(res1, res2)
                if min_dist <= self.contact_threshold:
                    contact_count += 1
        
        return contact_count
    
    def _analyze_graph(self, graph, graph_type):
        """Analyze graph properties."""
        if graph is None or not NETWORKX_AVAILABLE:
            return {'error': 'NetworkX not available or invalid graph'}
        
        if len(graph.nodes()) == 0:
            return {'error': 'Empty graph'}
        
        n_nodes = len(graph.nodes())
        analysis = {
            'node_count': n_nodes,
            'edge_count': len(graph.edges()),
            'density': nx.density(graph),
            'connected': nx.is_connected(graph),
            'average_degree': sum(dict(graph.degree()).values()) / n_nodes if n_nodes else 0
        }
        # Serializable edge list for output (node1, node2, and edge attributes)
        edge_list = []
        for u, v in graph.edges():
            ed = dict(graph.edges[u, v])
            edge_list.append({'node1': u, 'node2': v, **ed})
        analysis['edges'] = edge_list

        if nx.is_connected(graph):
            analysis.update({
                'diameter': nx.diameter(graph),
                'average_path_length': nx.average_shortest_path_length(graph),
                'clustering_coefficient': nx.average_clustering(graph)
            })
        else:
            # Analyze largest connected component
            largest_cc = max(nx.connected_components(graph), key=len)
            subgraph = graph.subgraph(largest_cc)
            
            analysis.update({
                'largest_component_size': len(largest_cc),
                'component_count': nx.number_connected_components(graph),
                'largest_component_diameter': nx.diameter(subgraph),
                'largest_component_avg_path': nx.average_shortest_path_length(subgraph),
                'clustering_coefficient': nx.average_clustering(graph)
            })
        
        # Centrality measures
        try:
            betweenness = nx.betweenness_centrality(graph)
            closeness = nx.closeness_centrality(graph)
            degree_centrality = nx.degree_centrality(graph)
            
            analysis.update({
                'max_betweenness': max(betweenness.values()) if betweenness else 0,
                'avg_betweenness': np.mean(list(betweenness.values())) if betweenness else 0,
                'max_closeness': max(closeness.values()) if closeness else 0,
                'avg_closeness': np.mean(list(closeness.values())) if closeness else 0,
                'max_degree_centrality': max(degree_centrality.values()) if degree_centrality else 0
            })
        except Exception:
            pass
        
        return analysis
    
    def _calculate_network_metrics(self, residue_graph, chain_graph):
        """Calculate comprehensive network metrics."""
        metrics = {}
        
        if NETWORKX_AVAILABLE and residue_graph is not None:
            # Modularity (if possible)
            try:
                # Simple community detection
                communities = nx.algorithms.community.greedy_modularity_communities(residue_graph)
                modularity = nx.algorithms.community.modularity(residue_graph, communities)
                metrics['modularity'] = modularity
                metrics['community_count'] = len(communities)
            except Exception:
                metrics['modularity'] = 0
                metrics['community_count'] = 0
            
            # Small world properties
            try:
                if nx.is_connected(residue_graph):
                    random_graph = nx.erdos_renyi_graph(
                        residue_graph.number_of_nodes(),
                        nx.density(residue_graph)
                    )
                    
                    if nx.is_connected(random_graph):
                        clustering_ratio = (nx.average_clustering(residue_graph) / 
                                          nx.average_clustering(random_graph))
                        path_ratio = (nx.average_shortest_path_length(residue_graph) /
                                    nx.average_shortest_path_length(random_graph))
                        
                        small_world = clustering_ratio / path_ratio if path_ratio > 0 else 0
                        metrics['small_world_coefficient'] = small_world
                        metrics['clustering_ratio'] = clustering_ratio
                        metrics['path_length_ratio'] = path_ratio
            except Exception:
                pass
        
        return metrics


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


LIMIT_INLINE = 80


def _run_analysis(analyzer, paths, out_stream, output_file_path=None):
    """Run graph analysis on paths and write results to out_stream."""
    for i, pdb_file in enumerate(paths):
        stem = os.path.splitext(os.path.basename(pdb_file))[0]
        base = (os.path.splitext(output_file_path)[0] + "_" + stem) if output_file_path else stem
        if len(paths) > 1:
            print(f"\n{'='*50}\n[{i+1}/{len(paths)}] {os.path.basename(pdb_file)}\n{'='*50}", file=out_stream)
        else:
            print(f"Analyzing graph properties of {pdb_file}...", file=out_stream)

        results = analyzer.analyze_packing_graph(pdb_file)

        if 'error' in results:
            print(f"Error: {results['error']}", file=out_stream)
            continue

        print("\nGRAPH ANALYSIS RESULTS:", file=out_stream)
        print("=" * 40, file=out_stream)
        res_graph = results.get('residue_graph', {})
        if 'error' not in res_graph and isinstance(res_graph, dict):
            print("\nResidue contact graph:", file=out_stream)
            print(f"  Nodes: {res_graph.get('node_count', 0)}  Edges: {res_graph.get('edge_count', 0)}", file=out_stream)
            print(f"  Density: {res_graph.get('density', 0):.3f}  Average degree: {res_graph.get('average_degree', 0):.2f}", file=out_stream)
            print(f"  Connected: {res_graph.get('connected', False)}", file=out_stream)
            print(f"  Clustering coefficient: {res_graph.get('clustering_coefficient', 0):.3f}", file=out_stream)
            if res_graph.get('connected'):
                print(f"  Diameter: {res_graph.get('diameter', 0)}  Avg path length: {res_graph.get('average_path_length', 0):.3f}", file=out_stream)
            else:
                print(f"  Components: {res_graph.get('component_count', 0)}  Largest component size: {res_graph.get('largest_component_size', 0)}", file=out_stream)
                print(f"  Largest component diameter: {res_graph.get('largest_component_diameter', 0)}  Avg path: {res_graph.get('largest_component_avg_path', 0):.3f}", file=out_stream)
            for k in ('max_betweenness', 'avg_betweenness', 'max_closeness', 'avg_closeness', 'max_degree_centrality'):
                if k in res_graph:
                    print(f"  {k}: {res_graph[k]:.3f}", file=out_stream)
            edges = res_graph.get('edges', [])
            if isinstance(edges, list) and edges:
                print(f"  Edges (node1, node2, distance Å, inter_chain):", file=out_stream)
                if len(edges) <= LIMIT_INLINE:
                    for e in edges:
                        if isinstance(e, dict):
                            dist = e.get('distance', '')
                            inter = e.get('inter_chain', '')
                            print(f"    {e.get('node1','')} -- {e.get('node2','')}  dist={dist}  inter_chain={inter}", file=out_stream)
                else:
                    for e in edges[:LIMIT_INLINE]:
                        if isinstance(e, dict):
                            dist = e.get('distance', '')
                            inter = e.get('inter_chain', '')
                            print(f"    {e.get('node1','')} -- {e.get('node2','')}  dist={dist}  inter_chain={inter}", file=out_stream)
                    extra_path = f"{base}_residue_edges.txt"
                    with open(extra_path, 'w') as ef:
                        ef.write(f"# All {len(edges)} residue graph edges: node1 node2 distance inter_chain\n")
                        for e in edges:
                            if isinstance(e, dict):
                                ef.write(f"{e.get('node1','')} {e.get('node2','')} {e.get('distance','')} {e.get('inter_chain','')}\n")
                    print(f"  Full list ({len(edges)} edges) written to {extra_path}", file=out_stream)
        chain_graph = results.get('chain_graph', {})
        if 'error' not in chain_graph and isinstance(chain_graph, dict):
            print("\nChain interaction graph:", file=out_stream)
            print(f"  Chains (nodes): {chain_graph.get('node_count', 0)}  Interactions (edges): {chain_graph.get('edge_count', 0)}", file=out_stream)
            print(f"  Density: {chain_graph.get('density', 0):.3f}", file=out_stream)
            edges = chain_graph.get('edges', [])
            if isinstance(edges, list) and edges:
                print(f"  Chain pairs (chain1, chain2, contact_count):", file=out_stream)
                for e in edges:
                    if isinstance(e, dict):
                        print(f"    {e.get('node1','')} -- {e.get('node2','')}  contacts={e.get('contact_count', 0)}", file=out_stream)
        metrics = results.get('network_metrics', {})
        if metrics and isinstance(metrics, dict):
            print("\nNetwork metrics:", file=out_stream)
            print(f"  Modularity: {metrics.get('modularity', 0):.3f}  Communities: {metrics.get('community_count', 0)}", file=out_stream)
            if 'small_world_coefficient' in metrics:
                print(f"  Small-world coefficient: {metrics.get('small_world_coefficient', 0):.3f}", file=out_stream)
            if 'clustering_ratio' in metrics:
                print(f"  Clustering ratio: {metrics.get('clustering_ratio', 0):.3f}  Path length ratio: {metrics.get('path_length_ratio', 0):.3f}", file=out_stream)

    if len(paths) > 1:
        print(f"\nDone. Processed {len(paths)} file(s).", file=out_stream)


def main():
    """CLI: single file, multiple files, directory, or glob patterns; optional multi-set filtering and dry-run."""
    import argparse
    import sys

    parser = argparse.ArgumentParser(
        description="Graph-theoretical analysis of crystal packing (residue/chain networks).",
        epilog="""Examples:
  python graph_analyzer.py structure.pdb
  python graph_analyzer.py dir/
  python graph_analyzer.py *.pdb -o results.txt
  python graph_analyzer.py *.pdb --set tag1,tag2 -o output_set.txt
  python graph_analyzer.py *.pdb --sets set1 set2 set3 -o "output_{}.txt"
  python graph_analyzer.py *.pdb --per-structure -o "{}_graph.txt"
  python graph_analyzer.py *.pdb --sets set1 set2 -o results.txt --dry-run
""",
    )
    parser.add_argument('input', nargs='+', help='PDB/CIF file(s), directory, or glob pattern (e.g. *.pdb)')
    parser.add_argument(
        '--set', '-s', action='append', dest='sets', metavar='PATTERNS',
        help='Comma-separated patterns; file included only if basename contains ALL. Repeat for multiple sets.',
    )
    parser.add_argument(
        '--sets', dest='sets_multi', nargs='+', metavar='SET',
        help='Multiple set names in one go (one pattern per set). E.g. --sets set1 set2 set3.',
    )
    parser.add_argument(
        '--output', '-o', metavar='FILE',
        help='Output file. Single file for all sets, or use "{}" for one file per set. Per-structure: use --per-structure with "{}" for file stem.',
    )
    parser.add_argument(
        '--per-structure', '-p', action='store_true',
        help='Write one output file per structure; -o must contain "{}" (replaced by file stem). Can combine with --set/--sets to filter which files are processed.',
    )
    parser.add_argument('--dry-run', action='store_true', help='Print which files would be processed per set and exit.')
    args = parser.parse_args()

    paths = collect_structure_paths(args.input)
    if not paths:
        print("No structure files found.", file=sys.stderr)
        sys.exit(1)

    # Per-structure mode: one output file per input file (optionally filtered by --set/--sets)
    if args.per_structure:
        if not args.output or '{}' not in args.output:
            print("Error: --per-structure requires -o with '{}' in the path (e.g. -o '{}_graph.txt').", file=sys.stderr)
            sys.exit(1)
        paths_to_process = list(paths)
        if getattr(args, 'sets_multi', None) or args.sets:
            filtered_set = set()
            if getattr(args, 'sets_multi', None):
                for name in args.sets_multi:
                    patterns = [name.strip()] if name.strip() else []
                    if patterns:
                        filtered_set.update(filter_paths_by_patterns(paths, patterns))
            if args.sets:
                for s in args.sets:
                    patterns = [p.strip() for p in s.split(',') if p.strip()]
                    if patterns:
                        filtered_set.update(filter_paths_by_patterns(paths, patterns))
            paths_to_process = sorted(filtered_set)
        if not paths_to_process:
            print("No structure files match the filter.", file=sys.stderr)
            sys.exit(1)
        if args.dry_run:
            print("Dry run (per-structure mode): no analysis will be performed.\n", file=sys.stderr)
            for p in paths_to_process:
                stem = os.path.splitext(os.path.basename(p))[0]
                dest = args.output.replace('{}', stem)
                print(f"  {os.path.basename(p)} -> {dest}", file=sys.stderr)
            return
        analyzer = GraphAnalyzer()
        for p in paths_to_process:
            stem = os.path.splitext(os.path.basename(p))[0]
            out_path = args.output.replace('{}', stem)
            print(f"Writing {os.path.basename(p)} -> {out_path}", file=sys.stderr)
            with open(out_path, 'w') as out:
                _run_analysis(analyzer, [p], out, output_file_path=None)
        return

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
            print(f"Set {label!r} (patterns: {patterns}): {len(filtered)} file(s) -> {dest}", file=sys.stderr)
            for p in filtered:
                print(f"  {os.path.basename(p)}", file=sys.stderr)
            print(file=sys.stderr)
        return

    single_output_file = args.output and '{}' not in args.output
    out = sys.stdout
    if single_output_file:
        out = open(args.output, 'w')
        print(f"Writing all sets to {args.output}", file=sys.stderr)

    analyzer = GraphAnalyzer()
    for idx, (label, patterns, filtered) in enumerate(set_list):
        if not filtered:
            print(f"No files match set {label!r} (patterns: {patterns}), skipping.", file=sys.stderr)
            continue

        out_path = None
        current_out = sys.stdout
        if single_output_file:
            current_out = out
        elif args.output:
            out_path = args.output.replace('{}', label)
            current_out = open(out_path, 'w')
            print(f"Writing set {label!r} ({len(filtered)} file(s)) to {out_path}", file=sys.stderr)
        elif len(set_list) > 1:
            print(f"\n=== Set {label!r} ({len(filtered)} file(s)) ===\n", file=sys.stderr)

        try:
            if len(set_list) > 1 and (single_output_file or not args.output):
                print(f"\n{'='*50}\nSet {label!r} (patterns: {patterns})\n{'='*50}", file=current_out)
            _run_analysis(analyzer, filtered, current_out, out_path if out_path is not None else (args.output if single_output_file else None))
        finally:
            if args.output and not single_output_file:
                current_out.close()

    if single_output_file:
        out.close()


if __name__ == "__main__":
    main() 