### Author: Edward Huang

import file_operations
import random
import math
import sys

### Creates a file for the network with the GO terms and the network without.
### Follows the following format:
### SPECIES_INDEX
### NUM_NODES
### gene_a gene_b weight
### gene_b gene_a weight
### Each edge twice. SPECIES_INDEX should just be 0 for single species.
### We can sample a subgraph to speed up clustering.

MIN_GO_SIZE = 10 # Minimum number of genes in a GO term to consider that term.
MAX_GO_SIZE = 1000

# Return a set of sampled genes from the sampled edges.
def get_sampled_genes(sampled_edges):
    sampled_genes = set([])
    for gene_a, gene_b in sampled_edges:
        sampled_genes.add(gene_a)
        sampled_genes.add(gene_b)
    return sampled_genes

# Find the size of the largest GO term.
def find_max_go_size(go_dct):
    max_go_size = 0
    for go in go_dct:
        num_go_genes = len(go_dct[go])
        if num_go_genes < MIN_GO_SIZE or num_go_genes > MAX_GO_SIZE:
            continue
        max_go_size = max(max_go_size, num_go_genes)
    return max_go_size

# Writes the networks without GO.
def write_no_go_files(run_num, num_genes, sampled_edges, edge_dct):
    no_go_out = open('./data/network_no_go_%s.txt' % (run_num), 'w')
    no_go_out.write('0\n%d\n' % num_genes)
    # Real network file for cluster evaluation.
    ng_real = open('./data/real_network_no_go_%s.txt' % (run_num), 'w')
    ng_real.write('Real network\n')
    for gene_a, gene_b in sampled_edges:
        weight = edge_dct[(gene_a, gene_b)]
        # Write in each edge twice to make it undirected.
        no_go_out.write('%s\t%s\t%s\n' % (gene_a, gene_b, weight))
        no_go_out.write('%s\t%s\t%s\n' % (gene_b, gene_a, weight))
        ng_real.write('0\t%s\t%s\t%s\n' % (gene_a, gene_b, weight))
        ng_real.write('0\t%s\t%s\t%s\n' % (gene_b, gene_a, weight))
    no_go_out.close()
    ng_real.close()

# Compute the GO weight.
def compute_go_weight(lamb, max_go_size, go_genes):
    return lamb * math.log(max_go_size / float(len(go_genes)))

# Writes the networks with GO.
def write_go_files(run_num, lamb, max_go_size, num_genes, go_dct, edge_dct,
    sampled_edges):
    go_out = open('./data/network_go_%s.txt' % (run_num), 'w')
    g_real = open('./data/real_network_go_%s.txt' % (run_num), 'w')
    g_real.write('Real network\n')
    go_out.write('0\n')
    # Number of nodes starts with number of genes. Increment each time we see
    # a valid GO term.
    num_nodes = num_genes
    for go in go_dct:
        go_genes = go_dct[go]
        num_go_genes = len(go_genes)
        if num_go_genes < MIN_GO_SIZE or num_go_genes > MAX_GO_SIZE:
            continue
        num_nodes += 1
    go_out.write('%d\n' % num_nodes)
    # Now write all of the gene-GO edges.
    for go in go_dct:
        go_genes = go_dct[go]
        num_go_genes = len(go_genes)
        if num_go_genes < MIN_GO_SIZE or num_go_genes > MAX_GO_SIZE:
            continue
        # Here we penalize GO terms that have many genes.
        go_weight = compute_go_weight(lamb, max_go_size, go_genes)

        for gene in go_genes:
            go_out.write('%s\t%s\t%f\n' % (gene, go, go_weight))
            go_out.write('%s\t%s\t%f\n' % (go, gene, go_weight))
            g_real.write('0\t%s\t%s\t%f\n' % (gene, go, go_weight))
            g_real.write('0\t%s\t%s\t%f\n' % (go, gene, go_weight))
    # Now write all of the gene-gene edges.
    for gene_a, gene_b in sampled_edges:
        weight = edge_dct[(gene_a, gene_b)]
        # Write in each edge twice to make it undirected.
        go_out.write('%s\t%s\t%s\n' % (gene_a, gene_b, weight))
        go_out.write('%s\t%s\t%s\n' % (gene_b, gene_a, weight))
        g_real.write('0\t%s\t%s\t%s\n' % (gene_a, gene_b, weight))
        g_real.write('0\t%s\t%s\t%s\n' % (gene_b, gene_a, weight))
    go_out.close()
    g_real.close()


if __name__ == '__main__':
    if len(sys.argv) != 5:
        print 'Usage:python %s run_num lambda subgraph_frac pearson/embedding' % sys.argv[0]
        exit()
    run_num = sys.argv[1]
    lamb = float(sys.argv[2])
    subgraph_frac = float(sys.argv[3]) # Fraction of graph to randomly sample.
    edge_method = sys.argv[4]
    assert edge_method in ['pearson', 'embedding']

    # Keys are pairs of genes, values are the edge weights.
    edge_dct = file_operations.get_raw_edge_dct()

    # Sample the fraction subgraph.
    num_samp_edges = int(math.ceil(subgraph_frac * len(edge_dct)))
    # We seed so that networks of the same percentage of raw network will have
    # the same non-GO edges.
    random.seed(5191993)
    sampled_edges = random.sample(edge_dct.keys(), num_samp_edges)

    # Keeps track of all the unique genes in the network.
    sampled_genes = get_sampled_genes(sampled_edges)

    # First count how many sampled nodes there are.
    num_genes = len(sampled_genes)

    if edge_method == 'embedding':
        edge_dct = file_operations.get_embedding_edge_dct()

    # Write networks without GO.
    write_no_go_files(run_num, num_genes, sampled_edges, edge_dct)

    # Get the GO labels.
    go_dct = file_operations.get_go_labels(sampled_genes)

    # Find the size of the largest GO term.
    max_go_size = find_max_go_size(go_dct)

    write_go_files(run_num, lamb, max_go_size, num_genes, go_dct, edge_dct,
        sampled_edges)