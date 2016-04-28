### Author: Edward Huang

import file_operations
import math
import operator
import sys
import time

### Creates a file for the network with the GO terms and the network without.
### Follows the following format:
### SPECIES_INDEX
### NUM_NODES
### gene_a gene_b weight
### gene_b gene_a weight
### Each edge twice. SPECIES_INDEX should just be 0 for single species.
### We can sample a subgraph to speed up clustering.
### Run time: 18 minutes.

MIN_GO_SIZE = 10 # Minimum number of genes in a GO term to consider that term.
MAX_GO_SIZE = 1000

# Return the set of genes from a set of edges.
def get_genes_from_edges(edges):
    genes = set([])
    for gene_a, gene_b in edges:
        genes.add(gene_a)
        genes.add(gene_b)
    return genes

# Find the size of the largest GO term. The second argument is the set of terms
# to actually iterate over.
def find_max_go_size(go_dct, go_terms):
    max_go_size = 0
    for go in go_dct:
        if go not in go_terms:
            continue
        num_go_genes = len(go_dct[go])
        if num_go_genes < MIN_GO_SIZE or num_go_genes > MAX_GO_SIZE:
            continue
        max_go_size = max(max_go_size, num_go_genes)
    return max_go_size

# Writes the networks without GO.
def write_no_go_files(run_num, num_genes, edges):
    no_go_out = open('./data/network_no_go_%s.txt' % (run_num), 'w')
    no_go_out.write('0\n%d\n' % num_genes)
    # Real network file for cluster evaluation.
    ng_real = open('./data/real_network_no_go_%s.txt' % (run_num), 'w')
    ng_real.write('Real network\n')
    for gene_a, gene_b in edges:
        # Write in each edge twice to make it undirected.
        no_go_out.write('%s\t%s\t1\n' % (gene_a, gene_b))
        no_go_out.write('%s\t%s\t1\n' % (gene_b, gene_a))
        ng_real.write('0\t%s\t%s\t1\n' % (gene_a, gene_b))
        ng_real.write('0\t%s\t%s\t1\n' % (gene_b, gene_a))
    no_go_out.close()
    ng_real.close()

# Compute the GO weight.
def compute_go_weight(lamb, max_go_size, num_go_genes):
    return lamb * math.log(max_go_size / float(num_go_genes))

# Writes the networks with GO.
def write_go_files(run_num, num_genes, edges, go_dct, lamb):
    # 3-fold clustering. Each sub-list contains GO terms.
    go_term_chunks = file_operations.chunkify_go_terms()

    for chunk_index, current_go_terms in enumerate(go_term_chunks):
        # Find the size of the largest GO term.
        max_go_size = find_max_go_size(go_dct, current_go_terms)

        go_out = open('./data/network_go_%s_%d.txt' % (run_num, chunk_index),
            'w')
        go_out.write('0\n')
        g_real = open('./data/real_network_go_%s_%d.txt' % (run_num,
            chunk_index), 'w')
        g_real.write('Real network\n')
        # Number of nodes starts with number of genes. Increment each time we
        # see a valid GO term.
        num_nodes = num_genes
        for go in current_go_terms:
            go_genes = go_dct[go]
            num_go_genes = len(go_genes)
            if num_go_genes < MIN_GO_SIZE or num_go_genes > MAX_GO_SIZE:
                continue
            num_nodes += 1
        go_out.write('%d\n' % num_nodes)
        # Now write all of the gene-GO edges.
        for go in current_go_terms:
            go_genes = go_dct[go]
            num_go_genes = len(go_genes)
            if num_go_genes < MIN_GO_SIZE or num_go_genes > MAX_GO_SIZE:
                continue
            # Here we penalize GO terms that have many genes.
            go_weight = compute_go_weight(lamb, max_go_size, num_go_genes)
            if go_weight == 0.0:
                continue
            for gene in go_genes:
                go_out.write('%s\t%s\t%.3f\n' % (gene, go, go_weight))
                go_out.write('%s\t%s\t%.3f\n' % (go, gene, go_weight))
                g_real.write('0\t%s\t%s\t%.3f\n' % (gene, go, go_weight))
                g_real.write('0\t%s\t%s\t%.3f\n' % (go, gene, go_weight))
        # Now write all of the gene-gene edges.
        for gene_a, gene_b in edges:
            # Write in each edge twice to make it undirected.
            go_out.write('%s\t%s\t1\n' % (gene_a, gene_b))
            go_out.write('%s\t%s\t1\n' % (gene_b, gene_a))
            g_real.write('0\t%s\t%s\t1\n' % (gene_a, gene_b))
            g_real.write('0\t%s\t%s\t1\n' % (gene_b, gene_a))
        go_out.close()
        g_real.close()

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s run_num' % sys.argv[0]
        exit()
    run_num = sys.argv[1]

    config_dct = file_operations.read_config_file()[run_num]
    lamb = float(config_dct['lamb'])
    subgraph_decimal = float(config_dct['subgraph_decimal'])
    edge_method = config_dct['edge_method']
    assert edge_method in ['pearson', 'embedding']

    # Keys are pairs of genes, values are the edge weights.
    # Right now, only looking at genes with high standard deviation.
    edge_dct = file_operations.get_raw_edge_dct()

    # Calculate how many edges to keep.
    num_edges = int(math.ceil(subgraph_decimal * len(edge_dct)))

    ### THIS IS WHERE YOU LEFT OFF.

    # Keep the top num_edges weights from the edge_dictionary, first sort.
    top_edges = sorted(edge_dct.items(), key=operator.itemgetter(1),
        reverse=True)[:num_edges]
    edges = [edge for (edge, weight) in top_edges]

    # Keeps track of all the unique genes in the network.
    genes = get_genes_from_edges(edges)

    # First count how many sampled nodes there are.
    num_genes = len(genes)

    if edge_method == 'embedding':
        edge_dct = file_operations.get_embedding_edge_dct()

    # Write networks without GO.
    write_no_go_files(run_num, num_genes, edges)

    # Get the GO labels.
    go_dct = file_operations.get_go_labels(genes)

    # write_go_files(run_num, num_genes, sampled_edges, edge_dct, go_dct, lamb)
    write_go_files(run_num, num_genes, edges, go_dct, lamb)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))