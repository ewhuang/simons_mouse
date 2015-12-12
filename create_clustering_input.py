### Author: Edward Huang

import random
import math

### Creates a file for the network with the GO terms and the network without.
### Follows the following format:
### SPECIES_INDEX
### NUM_NODES
### gene_a gene_b weight
### gene_b gene_a weight
### Each edge twice. SPECIES_INDEX should just be 0 for single species.
### We can sample a subgraph to speed up clustering.

SUBGRAPH_FRAC = 0.01 # Fraction of graph to randomly sample.
lamb = 1.5 # Tunable weight for all GO edge weights.
# Maybe punish big GO nodes by inversely weighting the lambda.

if __name__ == '__main__':
    # Keys are pairs of genes, values are the edge weights.
    edge_dct = {}
    print 'Extracting network data...'
    f = open('./data/raw_network.txt', 'r')
    for line in f:
        gene_a, gene_b, pcc = line.split()
        edge_dct[(gene_a, gene_b)] = pcc
    f.close()

    # Sample the fraction subgraph.
    num_samp_edges = int(math.ceil(SUBGRAPH_FRAC * len(edge_dct)))
    sampled_edges = random.sample(edge_dct.keys(), num_samp_edges)

    # Keeps track of all the unique genes in the network.
    genes = set([])
    for gene_a, gene_b in sampled_edges:
        genes.add(gene_a)
        genes.add(gene_b)

    print 'Writing to network without GO labels...'
    no_go_out = open('./data/network_no_go.txt', 'w')
    no_go_out.write('0\n%d\n' % len(genes))
    for gene_a, gene_b in sampled_edges:
        weight = edge_dct[(gene_a, gene_b)]
        # Write in each edge twice to make it undirected.
        no_go_out.write('%s\t%s\t%s\n' % (gene_a, gene_b, weight))
        no_go_out.write('%s\t%s\t%s\n' % (gene_b, gene_a, weight))
    no_go_out.close()

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

    print 'Writing to network with GO labels...'
    go_out = open('./data/network_go.txt', 'w')
    go_out.write('0\n')
    # First count how many nodes there are after adding in GO labels.
    num_nodes = len(genes)
    for go in go_dct:
        go_genes = go_dct[go]
        if len(go_genes) < 5 or len(go_genes) > 0.25 * len(genes):
            continue
        num_nodes += 1
    go_out.write('%d\n' % num_nodes)
    # Now write all of the gene-GO edges.
    for go in go_dct:
        go_genes = go_dct[go]
        if len(go_genes) < 5 or len(go_genes) > 0.25 * len(genes):
            continue
        for gene in go_genes:
            go_out.write('%s\t%s\t%f\n' % (gene, go, lamb))
            go_out.write('%s\t%s\t%f\n' % (go, gene, lamb))
    # Now write all of the gene-gene edges.
    for gene_a, gene_b in sampled_edges:
        weight = edge_dct[(gene_a, gene_b)]
        # Write in each edge twice to make it undirected.
        go_out.write('%s\t%s\t%s\n' % (gene_a, gene_b, weight))
        go_out.write('%s\t%s\t%s\n' % (gene_b, gene_a, weight))
    go_out.close()