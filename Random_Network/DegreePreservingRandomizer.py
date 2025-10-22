# -*- coding: utf-8 -*-

import networkx as nx
import pandas as pd
import random
import numpy as np
from collections import defaultdict
import sys

def load_network_data(file_path):
    """Load entire network data and count node information"""
    print("Loading network data...")
    
    # Read data
    df = pd.read_csv(file_path, sep='\t', header=None, names=['RE', 'Gene'])
    
    # Create list of all edges
    all_edges = list(zip(df['RE'], df['Gene']))
    
    print(f"Total network edges: {len(all_edges)}")
    
    # Count unique RE and gene nodes
    RE_nodes = set(df['RE'].unique())
    gene_nodes = set(df['Gene'].unique())
    
    print(f"RE nodes: {len(RE_nodes)}")
    print(f"Gene nodes: {len(gene_nodes)}")
    print(f"Total nodes: {len(RE_nodes) + len(gene_nodes)}")
    
    # Create directed graph
    full_network = nx.DiGraph()
    full_network.add_edges_from(all_edges)
    
    return full_network, RE_nodes, gene_nodes, len(all_edges)

def create_random_RE_gene_network(RE_nodes, gene_nodes, total_edges):
    """
    Generate random RE→gene network with precise edge count control
    Use two-phase approach: ensure connectivity first, then fill remaining edges
    """
    
    print(f"Generating random RE→gene network: {len(RE_nodes)} REs, {len(gene_nodes)} genes, {total_edges} edges")
    
    # Convert to lists for random selection
    RE_list = list(RE_nodes)
    gene_list = list(gene_nodes)
    
    random_network = nx.DiGraph()
    
    # Phase 1: Create basic connections to ensure each RE and gene has at least one edge
    print("Phase 1: Creating basic connections...")
    
    # Assign one outgoing edge for each RE node
    for re_node in RE_list:
        target_gene = random.choice(gene_list)
        random_network.add_edge(re_node, target_gene)
    
    # Assign one incoming edge for each gene node (if not already present)
    for gene_node in gene_list:
        if random_network.in_degree(gene_node) == 0:
            source_RE = random.choice(RE_list)
            random_network.add_edge(source_RE, gene_node)
    
    base_edges = random_network.number_of_edges()
    print(f"Base connection edges: {base_edges}")
    
    # Phase 2: Randomly add remaining edges, precisely controlling total edge count
    print("Phase 2: Randomly adding remaining edges...")
    
    remaining_edges = total_edges - base_edges
    
    if remaining_edges > 0:
        # Use set to track existing edges and avoid duplicates
        existing_edges = set(random_network.edges())
        edges_added = 0
        
        while edges_added < remaining_edges:
            source = random.choice(RE_list)
            target = random.choice(gene_list)
            
            if (source, target) not in existing_edges:
                random_network.add_edge(source, target)
                existing_edges.add((source, target))
                edges_added += 1
                
                if edges_added % 100000 == 0:
                    print(f"  Added {edges_added}/{remaining_edges} edges")
    
    # Verify final edge count
    final_edges = random_network.number_of_edges()
    print(f"Final edges: {final_edges} (Target: {total_edges})")
    
    if final_edges != total_edges:
        print(f"Warning: Edge count mismatch! Actual: {final_edges}, Target: {total_edges}")
    
    # Verify network connectivity
    isolated_RE = [node for node in RE_list if random_network.out_degree(node) == 0]
    isolated_genes = [node for node in gene_list if random_network.in_degree(node) == 0]
    
    print(f"Isolated RE nodes: {len(isolated_RE)}")
    print(f"Isolated gene nodes: {len(isolated_genes)}")
    
    return random_network

def calculate_basic_metrics(G, RE_nodes, gene_nodes, name=""):
    """Calculate basic network metrics"""
    print(f"\n{name} network metrics:")
    print(f"  Nodes: {G.number_of_nodes()}")
    print(f"  Edges: {G.number_of_edges()}")
    
    if G.is_directed() and G.number_of_edges() > 0:
        # Calculate degree information for REs and genes
        RE_out_degrees = []
        gene_in_degrees = []
        
        for node in G.nodes():
            if node in RE_nodes:
                RE_out_degrees.append(G.out_degree(node))
            elif node in gene_nodes:
                gene_in_degrees.append(G.in_degree(node))
        
        if RE_out_degrees:
            print(f"  RE average out-degree: {np.mean(RE_out_degrees):.2f} ± {np.std(RE_out_degrees):.2f}")
            print(f"  RE maximum out-degree: {max(RE_out_degrees)}")
            print(f"  RE minimum out-degree: {min(RE_out_degrees)}")
            print(f"  REs with zero out-degree: {sum(1 for d in RE_out_degrees if d == 0)}")
        
        if gene_in_degrees:
            print(f"  Gene average in-degree: {np.mean(gene_in_degrees):.2f} ± {np.std(gene_in_degrees):.2f}")
            print(f"  Gene maximum in-degree: {max(gene_in_degrees)}")
            print(f"  Gene minimum in-degree: {min(gene_in_degrees)}")
            print(f"  Genes with zero in-degree: {sum(1 for d in gene_in_degrees if d == 0)}")
    
    print(f"  Network density: {nx.density(G):.6f}")

def save_network_to_file(G, filename):
    """Save network as edge list file"""
    with open(filename, 'w') as f:
        for edge in G.edges():
            f.write(f"{edge[0]}\t{edge[1]}\n")
    print(f"Network saved to: {filename}")
    print(f"Actual saved edges: {G.number_of_edges()}")

def compare_degree_distributions(original_G, random_G, RE_nodes, gene_nodes):
    """Compare degree distributions between original and random networks"""
    print("\n" + "="*50)
    print("Degree distribution comparison:")
    
    # Original network degree distribution
    original_RE_out_degrees = [original_G.out_degree(node) for node in RE_nodes]
    original_gene_in_degrees = [original_G.in_degree(node) for node in gene_nodes]
    
    # Random network degree distribution
    random_RE_out_degrees = [random_G.out_degree(node) for node in RE_nodes]
    random_gene_in_degrees = [random_G.in_degree(node) for node in gene_nodes]
    
    print(f"Original network - RE average out-degree: {np.mean(original_RE_out_degrees):.2f} ± {np.std(original_RE_out_degrees):.2f}")
    print(f"Random network - RE average out-degree: {np.mean(random_RE_out_degrees):.2f} ± {np.std(random_RE_out_degrees):.2f}")
    
    print(f"Original network - Gene average in-degree: {np.mean(original_gene_in_degrees):.2f} ± {np.std(original_gene_in_degrees):.2f}")
    print(f"Random network - Gene average in-degree: {np.mean(random_gene_in_degrees):.2f} ± {np.std(random_gene_in_degrees):.2f}")

def main():
    """Main function"""
    # 1. Load data and count node information
    file_path = sys.argv[1]
    full_network, RE_nodes, gene_nodes, total_edges = load_network_data(file_path)
    
    # 2. Calculate original network metrics
    print("\n" + "="*50)
    calculate_basic_metrics(full_network, RE_nodes, gene_nodes, "Original")
    
    # 3. Generate RE→gene random network
    print("\n" + "="*50)
    print("Generating RE→gene random network...")
    random_network = create_random_RE_gene_network(RE_nodes, gene_nodes, total_edges)
    
    # 4. Calculate random network metrics
    print("\n" + "="*50)
    calculate_basic_metrics(random_network, RE_nodes, gene_nodes, "Random")
    
    # 5. Compare degree distributions
    compare_degree_distributions(full_network, random_network, RE_nodes, gene_nodes)
    
    # 6. Save random network
    print("\n" + "="*50)
    save_network_to_file(random_network, sys.argv[2])
    
    print("\nRandom network generation completed!")
    print(f"Generated network contains:")
    print(f"  - RE nodes: {len(RE_nodes)}")
    print(f"  - Gene nodes: {len(gene_nodes)}") 
    print(f"  - Edges: {random_network.number_of_edges()} (Target: {total_edges})")
    print(f"  - All edges are directed RE→gene edges")

if __name__ == "__main__":
    # Run main program
    main()