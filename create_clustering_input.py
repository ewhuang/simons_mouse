### Author: Edward Huang

import file_operations
import json
import math
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

def write_dca_gene_file(edge_genes):
    # First, write out the genes in its own separate file for DCA.
    # Make the folder if it doesn't exist.
    if not os.path.exists(dca_no_go_subfolder):
        os.makedirs(dca_no_go_subfolder)

    dca_genes_out = open('%s/dca_genes_no_go_%s.txt' % (dca_no_go_subfolder,
        run_num), 'w')
    # Write out the genes in its own file.
    for gene in edge_genes:
        dca_genes_out.write('%s\tA\n' % gene)
    dca_genes_out.close()

def write_orth_file(gene_example):
    '''
    Writes out the orthology file.
    '''
    # Write out the orth file.
    orth_out = open('./data/%s_data/orth.txt' % data_type, 'w')
    orth_out.write('0\t0\t%s\t%s' % (gene_example, gene_example))
    orth_out.close()

def write_no_go_files(edge_genes, edge_dct):
    '''
    Writes the networks without GO.
    '''
    no_go_folder = './data/%s_data/networks_no_go' % data_type
    # Regular network file for clustering.
    no_go_out = open('%s/network_no_go_%s.txt' % (no_go_folder, run_num), 'w')
    no_go_out.write('0\n%d\n' % len(edge_genes))
    # Real network file for cluster evaluation.
    ng_real = open('%s/real_network_no_go_%s.txt' % (no_go_folder, run_num),
        'w')
    ng_real.write('Real network\n')
    # Edge network for DCA.
    dca_edges_out = open('%s/dca_edges_no_go_%s.txt' % (dca_no_go_subfolder,
        run_num), 'w')

    for gene_a, gene_b in edge_dct:
        # Write in each edge twice to make it undirected. Edge weights are 1.
        edge_weight = edge_dct[(gene_a, gene_b)]
        forward, reverse = (gene_a, gene_b, edge_weight), (gene_b, gene_a,
            edge_weight)
        no_go_out.write('%s\t%s\t%s\n%s\t%s\t%s\n' % (gene_a, gene_b,
            edge_weight, gene_b, gene_a, edge_weight))
        ng_real.write('0\t%s\t%s\t%s\n0\t%s\t%s\t%s\n' % (gene_a, gene_b,
            edge_weight, gene_b, gene_a, edge_weight))
        dca_edges_out.write('%s\t%s\t%s\tn\n' % (gene_a, gene_b, edge_weight))
    no_go_out.close()
    ng_real.close()
    dca_edges_out.close()

    write_orth_file(gene_a)

def get_mf_dct():
    '''
    Returns the MF annotation dictionary.
    '''
    with open('./data/%s_data/mf_dct.json' % data_type, 'r') as fp:
        mf_go_gene_dct = json.load(fp)
    fp.close()

    return mf_go_gene_dct

def get_go_go_edges():
    '''
    Returns a dictionary.
    Key: GO:XXXXX -> str
    Value: list of neighboring GO terms -> list(str)
    '''
    with open('./data/mf_go_go_dct.json', 'r') as fp:
        mf_go_go_dct = json.load(fp)
    fp.close()

    return mf_go_go_dct

def bootstrap_reduce_go_dct(go_dct, edge_genes):
    new_go_dct = {}
    for go_term in go_dct:
        original_gene_set = go_dct[go_term]
        new_go_dct[go_term] = set(edge_genes).intersection(original_gene_set)
    return new_go_dct

def get_go_size_dct(mf_dct):
    '''
    Returns a (dictionary, int) pair.
    Dictionary is a GO to (size of GO) dct. Int is the size of the largest GO.
    Key: GO term -> str
    Value: Number of genes annotated by the key -> int
    '''
    go_size_dct = {}
    for go in mf_dct:
        num_go_genes = len(mf_dct[go])
        # Skip GOs that have sizes outside of our specifications in the config.
        if num_go_genes < min_go_size or num_go_genes > max_go_size:
            continue
        go_size_dct[go] = num_go_genes
    # Float the maximum size for later computations.
    return go_size_dct, float(max(go_size_dct.values()))

def write_go_files(edge_genes, edge_dct, bootstrap_idx=0):
    '''
    Writes the networks with GO.
    '''
    # First, write out the genes in its own separate file for DCA.
    # Make the folder if it doesn't exist.
    if not os.path.exists(dca_go_subfolder):
        os.makedirs(dca_go_subfolder)
    dca_genes_out = open('%s/dca_genes_go_%s.txt' % (dca_go_subfolder, run_num),
        'w')
    # Write out the genes in its own file.
    for gene in edge_genes:
        dca_genes_out.write('%s\tA\n' % gene)
    dca_edges_out = open('%s/dca_edges_go_%s.txt' % (dca_go_subfolder, run_num),
        'w')

    # Extract the MF GO dictionary and GO-GO edges.
    mf_dct, mf_go_go_dct = get_mf_dct(), get_go_go_edges()

    # TODO. Fix this to work with the mf_dct (bootstrapping stuff).
    if bootstrap:
        mf_dct = bootstrap_reduce_go_dct(mf_dct, edge_genes)

    go_size_dct, largest_go_size = get_go_size_dct(mf_dct)

    # Open up and initialize the network files for the current domain.
    if bootstrap:
        go_folder = './data/%s_data/bootstrapped_networks_go' % data_type
        go_out = open('%s/network_go_%s_%d.txt' % (go_folder, run_num,
            bootstrap_idx), 'w')
        g_real = open('%s/real_network_go_%s_%d.txt' % (go_folder, run_num,
            bootstrap_idx), 'w')
    else:
        go_folder = './data/%s_data/networks_go' % data_type
        go_out = open('%s/network_go_%s.txt' % (go_folder, run_num), 'w')
        g_real = open('%s/real_network_go_%s.txt' % (go_folder, run_num), 'w')
    # Write out the total number of nodes (genes + GO terms).
    go_out.write('0\n%d\n' % (len(edge_genes) + len(go_size_dct)))
    g_real.write('Real network\n')

    # Write GO-GO edges. Only write them for the DCA network.
    for go in mf_go_go_dct:
        if go not in go_size_dct:
            continue
        go_neighbor_list = mf_go_go_dct[go]
        for go_neighbor in go_neighbor_list:
            if go_neighbor not in go_size_dct:
                continue
            dca_edges_out.write('%s\t%s\t1\tg\n' % (go, go_neighbor))

    # Write gene-GO edges.
    for go in go_size_dct:
        dca_genes_out.write('%s\tB\n' % go)

        num_go_genes = go_size_dct[go]
        for gene in mf_dct[go]:
            if gene not in edge_genes:
                continue
            # Write in each edge twice to make it undirected.
            go_weight = max(lamb * math.log(largest_go_size / num_go_genes), 1)
            go_out.write('%s\t%s\t%f\n%s\t%s\t%f\n' % (gene, go, go_weight, go,
                gene, go_weight))
            g_real.write('0\t%s\t%s\t%f\n0\t%s\t%s\t%f\n' % (gene, go,
                go_weight, go, gene, go_weight))
            dca_edges_out.write('%s\t%s\t%f\to\n' % (gene, go, go_weight))

    # Write gene-gene edges.
    for gene_a, gene_b in edge_dct:
        edge_weight = float(edge_dct[(gene_a, gene_b)])
        # Skip 0 weighted edges.
        if edge_weight == 0:
            continue
        # Write in each edge twice to make it undirected.
        go_out.write('%s\t%s\t%f\n%s\t%s\t%f\n' % (gene_a, gene_b, edge_weight,
            gene_b, gene_a, edge_weight))
        g_real.write('0\t%s\t%s\t%f\n0\t%s\t%s\t%f\n' % (gene_a, gene_b,
            edge_weight, gene_b, gene_a, edge_weight))
        dca_edges_out.write('%s\t%s\t%f\tn\n' % (gene_a, gene_b, edge_weight))
    
    go_out.close()
    g_real.close()
    dca_edges_out.close()
    dca_genes_out.close()

def main():
    if len(sys.argv) not in [3, 4]:
        print ('Usage:python %s mouse/tcga_cancer_index run_num -b<bootstrap>'
            % sys.argv[0])
        exit()
    global data_type, run_num, bootstrap, lamb, max_go_size, min_go_size
    data_type = sys.argv[1]
    assert data_type == 'mouse' or data_type.isdigit()
    run_num = sys.argv[2]
    assert run_num.isdigit()
    bootstrap = '-b' in sys.argv

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    global dca_no_go_subfolder, dca_go_subfolder
    dca_no_go_subfolder = './data/%s_data/dca_networks_no_go' % data_type
    dca_go_subfolder = './data/%s_data/dca_networks_go' % data_type

    # Extracting configuration options.
    config_dct = file_operations.read_config_file(data_type)[run_num]
    lamb, max_go_size = config_dct['lamb'], config_dct['max_go_size']
    min_go_size = config_dct['min_go_size']

    high_std_edge_dct = file_operations.get_high_std_edge_dct(data_type)
    
    # Get the unique genes in the network.
    edge_genes = list(get_genes_from_edges(high_std_edge_dct.keys()))
    # Write out the DCA gene file for use with PROSNET.
    write_dca_gene_file(edge_genes)

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