### Author: Edward Huang

import file_operations
import json
import math
import operator
import sys
import time

### Creates files for the network with the GO terms and the network without.
### Follows the following format:
### SPECIES_INDEX
### NUM_NODES
### gene_a gene_b weight
### gene_b gene_a weight
### Each edge twice. SPECIES_INDEX should just be 0 for single species.
### We can sample a subgraph to speed up clustering.
### Run time: 10 minutes.

MIN_GO_SIZE = 10 # Minimum number of genes in a GO term to consider that term.
MAX_GO_SIZE = 1000

# Return the set of genes from a set of edges.
def get_genes_from_edges(edges):
    genes = set([])
    for gene_a, gene_b in edges:
        genes.add(gene_a)
        genes.add(gene_b)
    return genes

# Writes the networks without GO.
def write_no_go_files(run_num, num_genes, edges):
    no_go_folder = './data/networks_no_go/'

    # Regular network file for clustering.
    no_go_out = open('%snetwork_no_go_%s.txt' % (no_go_folder, run_num), 'w')
    no_go_out.write('0\n%d\n' % num_genes)
    # Real network file for cluster evaluation.
    ng_real = open('%sreal_network_no_go_%s.txt' % (no_go_folder, run_num), 'w')
    ng_real.write('Real network\n')
    for gene_a, gene_b in edges:
        # Write in each edge twice to make it undirected. Edge weights are 1.
        no_go_out.write('%s\t%s\t1\n' % (gene_a, gene_b))
        no_go_out.write('%s\t%s\t1\n' % (gene_b, gene_a))
        ng_real.write('0\t%s\t%s\t1\n' % (gene_a, gene_b))
        ng_real.write('0\t%s\t%s\t1\n' % (gene_b, gene_a))
    no_go_out.close()
    ng_real.close()

# Find the size of the largest GO term. The second argument is the set of terms
# to actually iterate over.
def find_max_go_size(go_dct):
    max_go_size = 0
    for go in go_dct:
        num_go_genes = len(go_dct[go])
        if num_go_genes < MIN_GO_SIZE or num_go_genes > MAX_GO_SIZE:
            continue
        max_go_size = max(max_go_size, num_go_genes)
    return max_go_size

# Compute the GO weight. Weights are at least 0.5.
def compute_go_weight(lamb, max_go_size, num_go_genes):
    return lamb * math.log(max_go_size / float(num_go_genes))

def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z

# Writes the networks with GO.
def write_go_files(run_num, num_genes, edges, lamb):

    # First, load all of the GO dictionaries.
    with open('./data/biological_process.json', 'r') as fp:
        bp_go_gene_dct = json.load(fp)
    fp.close()

    with open('./data/cellular_component_go.json', 'r') as fp:
        cc_go_gene_dct = json.load(fp)
    fp.close()

    with open('./data/molecular_function.json', 'r') as fp:
        mf_go_gene_dct = json.load(fp)
    fp.close()

    domain_dictionary_list = [bp_go_gene_dct, cc_go_gene_dct, mf_go_gene_dct]

    # 3-fold clustering. Each domain contains GO terms.
    num_domains = len(domain_dictionary_list)
    for domain_index in range(num_domains):
        # domain_index is the index of the domain that we do not add into the 
        # network, and that we evaluate on.

        # Cluster on the two domains that we do not evaluate on.
        go_dct_list = []
        for i in [num for num in range(num_domains) if num != domain_index]:
            go_dct_list += [domain_dictionary_list[i]]
        go_dct = merge_two_dicts(go_dct_list[0], go_dct_list[1])

        # Find the size of the largest GO term.
        max_go_size = find_max_go_size(go_dct)

        # Number of nodes is number of genes + number of GO terms.
        num_nodes = num_genes

        # Calculate the weights for each gene-GO edge.
        gene_go_weights = {}
        for go in go_dct:
            go_genes = go_dct[go]
            num_go_genes = len(go_genes)
            if num_go_genes < MIN_GO_SIZE or num_go_genes > MAX_GO_SIZE:
                continue
            
            # Calculate the weights between a GO term and all of its genes.
            go_weight = compute_go_weight(lamb, max_go_size, num_go_genes)

            if go_weight < 0.5:
                continue

            for gene in go_genes:
                gene_go_weights[(gene, go)] = go_weight
            num_nodes += 1

        # Open up and initialize the network files for the current domain.
        go_folder = './data/networks_go/'
        go_out = open('%snetwork_go_%s_%d.txt' % (go_folder, run_num,
            domain_index), 'w')
        go_out.write('0\n')
        g_real = open('%sreal_network_go_%s_%d.txt' % (go_folder, run_num,
            domain_index), 'w')
        g_real.write('Real network\n')
        go_out.write('%d\n' % num_nodes)

        # Write gene-GO edges.
        for gene, go in gene_go_weights:
            go_weight = gene_go_weights[(gene, go)]
            go_out.write('%s\t%s\t%.3f\n' % (gene, go, go_weight))
            go_out.write('%s\t%s\t%.3f\n' % (go, gene, go_weight))
            g_real.write('0\t%s\t%s\t%.3f\n' % (gene, go, go_weight))
            g_real.write('0\t%s\t%s\t%.3f\n' % (go, gene, go_weight))

        # Write all of the gene-gene edges.
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

    # Extracting configuration options.
    config_dct = file_operations.read_config_file()[run_num]
    lamb = float(config_dct['lamb'])
    edge_method = config_dct['edge_method']
    assert edge_method in ['pearson', 'embedding']

    high_std_edge_list = file_operations.get_high_std_edge_list()
    
    # Count the number of unique genes in the network.
    num_genes = len(get_genes_from_edges(high_std_edge_list))

    # if edge_method == 'embedding':
    #     edge_dct = file_operations.get_embedding_edge_dct()

    # Write networks without GO.
    write_no_go_files(run_num, num_genes, high_std_edge_list)
    write_go_files(run_num, num_genes, high_std_edge_list, lamb)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))