### Author: Edward Huang

import file_operations
import json
import math
import operator
import os
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
    # First, write out the genes in its own separate file for DCA.
    # Make the folder if it doesn't exist.
    dca_subfolder = './data/%s_data/dca_networks_no_go' % data_type
    if not os.path.exists(dca_subfolder):
        os.makedirs(dca_subfolder)

    dca_genes_out = open('%s/dca_genes_%s.txt' % (dca_subfolder, run_num), 'w')
    # Write out the genes in its own file.
    for gene in edge_genes:
        dca_genes_out.write('%s\tA\n' % gene)
    dca_genes_out.close()

    no_go_folder = './data/%s_data/networks_no_go' % data_type
    # Regular network file for clustering.
    no_go_out = open('%s/network_no_go_%s.txt' % (no_go_folder, run_num), 'w')
    no_go_out.write('0\n%d\n' % len(edge_genes))
    # Real network file for cluster evaluation.
    ng_real = open('%s/real_network_no_go_%s.txt' % (no_go_folder, run_num),
        'w')
    ng_real.write('Real network\n')
    # Edge network for DCA.
    dca_edges_out = open('%s/dca_edges_%s.txt' % (dca_subfolder, run_num), 'w')

    for gene_a, gene_b in edge_dct:
        # Write in each edge twice to make it undirected. Edge weights are 1.
        edge_weight = edge_dct[(gene_a, gene_b)]
        no_go_out.write('%s\t%s\t%s\n' % (gene_a, gene_b, edge_weight))
        no_go_out.write('%s\t%s\t%s\n' % (gene_b, gene_a, edge_weight))
        ng_real.write('0\t%s\t%s\t%s\n' % (gene_a, gene_b, edge_weight))
        ng_real.write('0\t%s\t%s\t%s\n' % (gene_b, gene_a, edge_weight))
        dca_edges_out.write('%s\t%s\t%s\tn\n' % (gene_a, gene_b, edge_weight))
    no_go_out.close()
    ng_real.close()
    dca_edges_out.close()

    # Write out the orth file.
    orth_out = open('./data/%s_data/orth.txt' % data_type, 'w')
    orth_out.write('0\t0\t%s\t%s' % (gene_a, gene_a))
    orth_out.close()

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
    with open('./data/%s_data/bp_dct.json' % data_type, 'r') as fp:
        bp_go_gene_dct = json.load(fp)
    fp.close()

    with open('./data/%s_data/mf_dct.json' % data_type, 'r') as fp:
        mf_go_gene_dct = json.load(fp)
    fp.close()

    return [bp_go_gene_dct, mf_go_gene_dct]

def get_go_go_edges():
    '''
    Returns a dictionary.
    Key: GO:XXXXX -> str
    Value: list of neighboring GO terms -> list(str)
    '''
    with open('./data/%s_data/mf_go_go_dct.json' % data_type, 'r') as fp:
        mf_go_go_dct = json.load(fp)
    fp.close()

    return mf_go_go_dct

def bootstrap_reduce_go_dct(go_dct, edge_genes):
    new_go_dct = {}
    for go_term in go_dct:
        original_gene_set = go_dct[go_term]
        new_go_dct[go_term] = set(edge_genes).intersection(original_gene_set)
    return new_go_dct

def write_go_files(edge_genes, edge_dct, bootstrap_idx=0):
    '''
    Writes the networks with GO.
    '''
    # First, write out the genes in its own separate file for DCA.
    # Make the folder if it doesn't exist.
    dca_subfolder = './data/%s_data/dca_networks_go' % data_type
    if not os.path.exists(dca_subfolder):
        os.makedirs(dca_subfolder)

    # We specify MF go because the domain_index variable only loops through [0].
    dca_genes_out = open('%s/dca_genes_mf_go_%s.txt' % (dca_subfolder, run_num),
        'w')
    # Write out the genes in its own file.
    for gene in edge_genes:
        dca_genes_out.write('%s\tA\n' % gene)
    dca_edges_out = open('%s/dca_edges_mf_go_%s.txt' % (dca_subfolder, run_num),
        'w')

    # Extract the three GO dictionaries.
    domain_dictionary_list = get_go_dictionaries()
    # TODO. Read overlap list in order to add in only partial MF terms.
    # overlap_list = file_operations.read_go_overlap()

    # Get the GO-GO edges.
    mf_go_go_dct = get_go_go_edges()

    # Each domain contains GO terms. domain_index indicates the index of the GO
    # domain we are evaluating on.
    for domain_index in [0]:
        # This line assumes we only use two domains. Train on the GO domain that
        # we aren't evaluating on.
        go_dct = domain_dictionary_list[1 - domain_index]

        # TODO Right now, adding in all MF terms, and then evaluating on
        # partial BP terms.
        # Remove the terms from MF that overlap too much with terms from BP.
        # overlapping_go_terms = set([tup[1 - domain_index] for tup in
        #     overlap_list])
        # for overlapping_go in overlapping_go_terms:
        #     del go_dct[overlapping_go]

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
            go_folder = './data/%s_data/bootstrapped_networks_go/' % data_type
            go_out = open('%snetwork_go_%s_%d_%d.txt' % (go_folder, run_num,
                domain_index, bootstrap_idx), 'w')
            g_real = open('%sreal_network_go_%s_%d_%d.txt' % (go_folder,
                run_num, domain_index, bootstrap_idx), 'w')
        else:
            go_folder = './data/%s_data/networks_go/' % data_type
            go_out = open('%snetwork_go_%s_%d.txt' % (go_folder, run_num,
                domain_index), 'w')
            g_real = open('%sreal_network_go_%s_%d.txt' % (go_folder, run_num,
                domain_index), 'w')
        go_out.write('0\n%d\n' % (len(edge_genes) + len(go_size_dct)))
        g_real.write('Real network\n')

        # Write GO-GO edges. Don't write them for simulated annealing.
        for go in mf_go_go_dct:
            if go not in go_size_dct:
                continue
            go_neighbor_list = mf_go_go_dct[go]
            for go_neighbor in go_neighbor_list:
                if go_neighbor not in go_size_dct:
                    continue
                    # Can we add weights into DCA?
                dca_edges_out.write('%s\t%s\t1\tg\n' % (go, go_neighbor))

        # Write gene-GO edges.
        for go in go_size_dct:
            dca_genes_out.write('%s\tB\n' % go)

            num_go_genes = go_size_dct[go]
            for gene in go_dct[go]:
                if gene not in edge_genes:
                    continue
                go_weight = compute_go_weight(largest_go_size, num_go_genes)
                go_out.write('%s\t%s\t%f\n' % (gene, go, go_weight))
                go_out.write('%s\t%s\t%f\n' % (go, gene, go_weight))
                g_real.write('0\t%s\t%s\t%f\n' % (gene, go, go_weight))
                g_real.write('0\t%s\t%s\t%f\n' % (go, gene, go_weight))
                dca_edges_out.write('%s\t%s\t%f\to\n' % (gene, go,
                    go_weight))

        # Write gene-gene edges.
        for gene_a, gene_b in edge_dct:
            # Write in each edge twice to make it undirected.
            edge_weight = edge_dct[(gene_a, gene_b)]
            go_out.write('%s\t%s\t%s\n' % (gene_a, gene_b, edge_weight))
            go_out.write('%s\t%s\t%s\n' % (gene_b, gene_a, edge_weight))
            g_real.write('0\t%s\t%s\t%s\n' % (gene_a, gene_b, edge_weight))
            g_real.write('0\t%s\t%s\t%s\n' % (gene_b, gene_a, edge_weight))
            dca_edges_out.write('%s\t%s\t%s\tn\n' % (gene_a, gene_b,
                edge_weight))
        go_out.close()
        g_real.close()
    dca_edges_out.close()
    dca_genes_out.close()

def main():
    if len(sys.argv) not in [3, 4]:
        print 'Usage:python %s mouse/tcga_cancer_index run_num -b <bootstrap>' % sys.argv[0]
        exit()
    global data_type, run_num, bootstrap, lamb, max_go_size, min_go_size
    data_type = sys.argv[1]
    assert data_type == 'mouse' or data_type.isdigit()
    run_num = sys.argv[2]
    assert run_num.isdigit()
    bootstrap = '-b' in sys.argv

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    # Extracting configuration options.
    config_dct = file_operations.read_config_file(data_type)[run_num]
    lamb, max_go_size = config_dct['lamb'], config_dct['max_go_size']
    min_go_size = config_dct['min_go_size']

    high_std_edge_dct = file_operations.get_high_std_edge_dct(data_type)
    
    # Get the unique genes in the network.
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