### Author: Edward Huang

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

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print 'Usage:python %s run_num lambda subgraph_frac' % sys.argv[0]
        exit()
    run_num = sys.argv[1]
    lamb = float(sys.argv[2])
    subgraph_frac = float(sys.argv[3]) # Fraction of graph to randomly sample.

    # Keys are pairs of genes, values are the edge weights.
    edge_dct = {}
    print 'Extracting network data...'
    f = open('./data/raw_network.txt', 'r')
    for line in f:
        gene_a, gene_b, pcc = line.split()
        edge_dct[(gene_a, gene_b)] = pcc
    f.close()

    # Sample the fraction subgraph.
    num_samp_edges = int(math.ceil(subgraph_frac * len(edge_dct)))
    # We seed so that networks of the same percentage of raw network will have
    # the same non-GO edges.
    random.seed(5191993)
    sampled_edges = random.sample(edge_dct.keys(), num_samp_edges)

    # Keeps track of all the unique genes in the network.
    genes = set([])
    for gene_a, gene_b in sampled_edges:
        genes.add(gene_a)
        genes.add(gene_b)

    print 'Writing to network without GO labels...'
    no_go_out = open('./data/network_no_go_%s.txt' % (run_num), 'w')
    no_go_out.write('0\n%d\n' % len(genes))
    # Real network file... for cluster evaluation.
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

    print 'Extracting GO labels...'
    go_f = open('./data/go_edges.txt', 'r')
    go_dct = {}
    for line in go_f:
        gene, go = line.split()
        if gene not in genes:
            continue
        if go not in go_dct:
            go_dct[go] = [gene]
        else:
            go_dct[go] += [gene]
    go_f.close()

    # Find the size of the largest GO term.
    max_go_size = 0
    for go in go_dct:
        num_go_genes = len(go_dct[go])
        # if num_go_genes < MIN_GO_SIZE or num_go_genes > 0.25 * len(genes):
        if num_go_genes < MIN_GO_SIZE or num_go_genes > MAX_GO_SIZE:
            continue
        max_go_size = max(max_go_size, num_go_genes)

    print 'Writing to network with GO labels...'
    go_out = open('./data/network_go_%s.txt' % (run_num), 'w')
    g_real = open('./data/real_network_go_%s.txt' % (run_num), 'w')
    g_real.write('Real network\n')
    go_out.write('0\n')
    # First count how many nodes there are after adding in GO labels.
    num_nodes = len(genes)
    for go in go_dct:
        go_genes = go_dct[go]
        # if len(go_genes) < MIN_GO_SIZE or len(go_genes) > 0.25 * len(genes):
        if len(go_genes) < MIN_GO_SIZE or len(go_genes) > MAX_GO_SIZE:
            continue
        num_nodes += 1
    go_out.write('%d\n' % num_nodes)
    # Now write all of the gene-GO edges.
    for go in go_dct:
        go_genes = go_dct[go]
        # if len(go_genes) < MIN_GO_SIZE or len(go_genes) > 0.25 * len(genes):
        if len(go_genes) < MIN_GO_SIZE or len(go_genes) > MAX_GO_SIZE:
            continue
        # Here we penalize GO terms that have many genes.
        go_weight = max(math.log(lamb * max_go_size / float(len(go_genes))), 1.0)
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