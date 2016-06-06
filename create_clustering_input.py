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
### Run time: Depends on the Pearson coefficient threshold set in gene_gene_edge
### weights.py. 1-15 minutes.

MIN_GO_SIZE = 10
MAX_GO_SIZE = 1000

def get_genes_from_edges(edges):
    '''
    Return the set of genes from a set of edges.
    '''
    genes = set([])
    for gene_a, gene_b in edges:
        genes.add(gene_a)
        genes.add(gene_b)
    return genes

def write_no_go_files(run_num, edge_genes, edge_dct):
    '''
    Writes the networks without GO.
    '''
    no_go_folder = './data/networks_no_go/'

    # Regular network file for clustering.
    no_go_out = open('%snetwork_no_go_%s.txt' % (no_go_folder, run_num), 'w')
    no_go_out.write('0\n%d\n' % len(edge_genes))
    # Real network file for cluster evaluation.
    ng_real = open('%sreal_network_no_go_%s.txt' % (no_go_folder, run_num), 'w')
    ng_real.write('Real network\n')
    for gene_a, gene_b in edge_dct:
        # Write in each edge twice to make it undirected. Edge weights are 1.
        edge_weight = edge_dct[(gene_a, gene_b)]
        no_go_out.write('%s\t%s\t%s\n' % (gene_a, gene_b, edge_weight))
        no_go_out.write('%s\t%s\t%s\n' % (gene_b, gene_a, edge_weight))
        ng_real.write('0\t%s\t%s\t%s\n' % (gene_a, gene_b, edge_weight))
        ng_real.write('0\t%s\t%s\t%s\n' % (gene_b, gene_a, edge_weight))
    no_go_out.close()
    ng_real.close()

def find_max_go_size(go_dct):
    '''
    Find the size of the largest GO term. The second argument is the set of
    terms to actually iterate over.
    '''
    max_go_size = 0
    for go in go_dct:
        num_go_genes = len(go_dct[go])
        if num_go_genes < MIN_GO_SIZE or num_go_genes > MAX_GO_SIZE:
            continue
        max_go_size = max(max_go_size, num_go_genes)
    return max_go_size

def compute_go_weight(lamb, max_go_size, num_go_genes):
    '''
    Compute the GO weight. Weights are at least 0.5.
    '''
    return lamb * math.log(max_go_size / float(num_go_genes))

def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = {}
    z.update(x)
    z.update(y)
    return z

def get_go_dictionaries():
    '''
    Fetches the GO dictionaries from the three domains. Returns a list of the
    dictionaries, by alphabetical order.
    '''
    # First, load all of the GO dictionaries.
    with open('./data/bp_ensmusg.json', 'r') as fp:
        bp_go_gene_dct = json.load(fp)
    fp.close()

    with open('./data/cc_ensmusg.json', 'r') as fp:
        cc_go_gene_dct = json.load(fp)
    fp.close()

    with open('./data/mf_ensmusg.json', 'r') as fp:
        mf_go_gene_dct = json.load(fp)
    fp.close()

    return [bp_go_gene_dct, cc_go_gene_dct, mf_go_gene_dct]

def write_go_files(run_num, edge_genes, edge_dct, lamb):
    '''
    Writes the networks with GO.
    '''
    # Extract the three GO dictionaries.
    domain_dictionary_list = get_go_dictionaries()

    # 3-fold clustering. Each domain contains GO terms.
    for domain_index in range(3):
        # domain_index is the index of the domain that we do not add into the 
        # network, and that we evaluate on.

        # HOLD ONE OUT CODE.
        # Cluster on the two domains that we do not evaluate on.
        go_dct = {}
        for i in [num for num in range(3) if num != domain_index]:
            go_dct.update(domain_dictionary_list[i])
        # go_dct = domain_dictionary_list[domain_index]        

        go_size_dct = {} # Keys are GO terms, values are the sizes of the terms.
        max_go_size = 0 # Find the size of the largest GO term.
        for go in go_dct:
            num_go_genes = len(go_dct[go])
            if num_go_genes < MIN_GO_SIZE or num_go_genes > MAX_GO_SIZE:
                continue

            go_size_dct[go] = num_go_genes
            # Update the largest GO size.
            max_go_size = max(max_go_size, num_go_genes)

        # Open up and initialize the network files for the current domain.
        go_folder = './data/networks_go/'
        go_out = open('%snetwork_go_%s_%d.txt' % (go_folder, run_num,
            domain_index), 'w')
        go_out.write('0\n%d\n' % (len(edge_genes) + len(go_size_dct)))
        g_real = open('%sreal_network_go_%s_%d.txt' % (go_folder, run_num,
            domain_index), 'w')
        g_real.write('Real network\n')

        # Write gene-GO edges.
        for go in go_size_dct:
            num_go_genes = go_size_dct[go]
            for gene in go_dct[go]:
                if gene not in edge_genes:
                    continue
                go_weight = compute_go_weight(lamb, max_go_size, num_go_genes)
                go_out.write('%s\t%s\t%f\n' % (gene, go, go_weight))
                go_out.write('%s\t%s\t%f\n' % (go, gene, go_weight))
                g_real.write('0\t%s\t%s\t%f\n' % (gene, go, go_weight))
                g_real.write('0\t%s\t%s\t%f\n' % (go, gene, go_weight))

        # Write all of the gene-gene edges.
        for gene_a, gene_b in edge_dct:
            # Write in each edge twice to make it undirected.
            edge_weight = edge_dct[(gene_a, gene_b)]
            go_out.write('%s\t%s\t%s\n' % (gene_a, gene_b, edge_weight))
            go_out.write('%s\t%s\t%s\n' % (gene_b, gene_a, edge_weight))
            g_real.write('0\t%s\t%s\t%s\n' % (gene_a, gene_b, edge_weight))
            g_real.write('0\t%s\t%s\t%s\n' % (gene_b, gene_a, edge_weight))
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

    high_std_edge_dct = file_operations.get_high_std_edge_dct()
    
    # Count the number of unique genes in the network.
    edge_genes = get_genes_from_edges(high_std_edge_dct.keys())

    # if edge_method == 'embedding':
    #     edge_dct = file_operations.get_embedding_edge_dct()

    # Write networks without GO.
    write_no_go_files(run_num, edge_genes, high_std_edge_dct)
    write_go_files(run_num, edge_genes, high_std_edge_dct, lamb)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))