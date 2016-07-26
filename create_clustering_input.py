### Author: Edward Huang

import file_operations
import json
import math
import operator
import random
import sys
import time

### Creates files for the network with the GO terms and the network without.
### Follows the following format:
### SPECIES_INDEX
### NUM_NODES
### gene_a gene_b weight
### gene_b gene_a weight
### Each edge twice. SPECIES_INDEX should just be 0 for single species.

### real_network_go format:
### Real network
### 0   gene_a  gene_b  edge_weight
### 0   gene_b  gene_a  edge_weight
### Run time < 1 minute.

def get_genes_from_edges(edges):
    '''
    Return the set of genes from a set of edges.
    '''
    genes = set([])
    for gene_a, gene_b in edges:
        genes.add(gene_a)
        genes.add(gene_b)
    return genes

def write_no_go_files(edge_genes, edge_dct):
    '''
    Writes the networks without GO.
    '''
    no_go_folder = './data/%s_networks_no_go/' % data_type

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

def compute_go_weight(largest_go_size, num_go_genes):
    '''
    Compute the GO weight. Weights are at least 1.0.
    '''
    return max(lamb * math.log(largest_go_size / float(num_go_genes)), 1.0)

def get_go_dictionaries():
    '''
    Fetches the GO dictionaries from the three domains. Returns a list of the
    dictionaries, by alphabetical order.
    '''
    # First, load all of the GO dictionaries.
    with open('./data/bp_%s.json' % data_type, 'r') as fp:
        bp_go_gene_dct = json.load(fp)
    fp.close()

    with open('./data/mf_%s.json' % data_type, 'r') as fp:
        mf_go_gene_dct = json.load(fp)
    fp.close()

    return [bp_go_gene_dct, mf_go_gene_dct]

def bootstrap_reduce_go_dct(go_dct, edge_genes):
    new_go_dct = {}
    for go_term in go_dct:
        original_gene_set = go_dct[go_term]
        new_go_dct[go_term] = set(edge_genes).intersection(original_gene_set)
    return new_go_dct

def read_go_overlap():
    '''
    Gets the BP and MF GO terms that are too similar to each other. We exclude
    them from training, but still use them to evaluate.
    '''
    overlap_list = []
    f = open('./data/overlapping_bp_mf_go_labels_%s.txt' % data_type, 'r')
    for line in f:
        bp_label, mf_label, p_value = line.strip().split('\t')
        overlap_list += [(bp_label, mf_label)]
    f.close()
    return overlap_list

def write_go_files(edge_genes, edge_dct, bootstrap_idx=0):
    '''
    Writes the networks with GO.
    '''
    # Extract the three GO dictionaries.
    domain_dictionary_list = get_go_dictionaries()
    overlap_list = read_go_overlap()

    # Each domain contains GO terms. domain_index indicates the index of the GO
    # domain we are evaluating on.
    # for domain_index in range(len(domain_dictionary_list)):
    for domain_index in [0]:
        # This line assumes we only use two domains. Train on the GO domain that
        # we aren't evaluating on.
        go_dct = domain_dictionary_list[1 - domain_index]

        # Remove the terms from MF that overlap too much with terms from BP.
        overlapping_go_terms = set([tup[1 - domain_index] for tup in
            overlap_list])
        for overlapping_go in overlapping_go_terms:
            del go_dct[overlapping_go]

        if bootstrap:
            go_dct = bootstrap_reduce_go_dct(go_dct, edge_genes)

        go_size_dct = {} # Keys are GO terms, values are the sizes of the terms.
        largest_go_size = 0 # Find the size of the largest GO term.
        for go in go_dct:
            num_go_genes = len(go_dct[go])
            if num_go_genes < min_go_size or num_go_genes > max_go_size:
                continue

            go_size_dct[go] = num_go_genes
            # Update the largest GO size.
            largest_go_size = max(largest_go_size, num_go_genes)

        # Open up and initialize the network files for the current domain.
        if bootstrap:
            go_folder = './data/bootstrapped_%s_networks_go/' % data_type
            go_out = open('%snetwork_go_%s_%d_%d.txt' % (go_folder, run_num,
                domain_index, bootstrap_idx), 'w')
            g_real = open('%sreal_network_go_%s_%d_%d.txt' % (go_folder,
                run_num, domain_index, bootstrap_idx), 'w')
        else:
            go_folder = './data/%s_networks_go/' % data_type
            go_out = open('%snetwork_go_%s_%d.txt' % (go_folder, run_num,
                domain_index), 'w')
            g_real = open('%sreal_network_go_%s_%d.txt' % (go_folder, run_num,
                domain_index), 'w')
        go_out.write('0\n%d\n' % (len(edge_genes) + len(go_size_dct)))
        g_real.write('Real network\n')

        # Write gene-GO edges.
        for go in go_size_dct:
            num_go_genes = go_size_dct[go]
            for gene in go_dct[go]:
                if gene not in edge_genes:
                    continue
                go_weight = compute_go_weight(largest_go_size, num_go_genes)
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
    if len(sys.argv) not in [3, 4]:
        print 'Usage:python %s data_type run_num -b <bootstrap>' % sys.argv[0]
        exit()
    global data_type, run_num, bootstrap, lamb, max_go_size, min_go_size
    data_type = sys.argv[1]
    assert data_type in ['mouse', 'tcga']
    run_num = sys.argv[2]
    assert run_num.isdigit()
    bootstrap = '-b' in sys.argv

    # Extracting configuration options.
    config_dct = file_operations.read_config_file(data_type)[run_num]
    lamb, max_go_size = config_dct['lamb'], config_dct['max_go_size']
    min_go_size = config_dct['min_go_size']

    high_std_edge_dct = file_operations.get_high_std_edge_dct(data_type)
    
    # Count the number of unique genes in the network.
    edge_genes = list(get_genes_from_edges(high_std_edge_dct.keys()))

    if bootstrap:
        for bootstrap_idx in range(1000):
            # If in bootstrap mode, sample 90% of the edges.
            new_high_std_edge_dct = dict(random.sample(
                high_std_edge_dct.items(), int(0.9 * len(high_std_edge_dct))))
            write_go_files(edge_genes, new_high_std_edge_dct, bootstrap_idx)
    else:
        # Write networks without GO.
        write_no_go_files(edge_genes, high_std_edge_dct)
        write_go_files(edge_genes, high_std_edge_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))