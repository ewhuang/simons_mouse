### Author: Edward Huang

import file_operations
import json
import operator
import os
from scipy.stats import fisher_exact
import sys
import time

### This script calculates the GO enrichment score for each clustering
### result. For each method, for every cluster, for every possible GO
### label's related genes, compute Fisher's test. Thus, for every cluster,
### obtain the lowest possible p-value corresponding to that cluster's
### most related GO label. These p-values (GO enrichments) should improve
### as a whole going from coexpression network to GO network without
### sacrificing in-density.
### Run time: 2 minutes.

def get_bp_dct(base_data_type):
    '''
   Gets the BP dictionary.
    '''
    with open('./data/%s_data/bp_dct.json' % base_data_type, 'r') as fp:
        bp_dct = json.load(fp)
    fp.close()
    return bp_dct

def read_dbgap_file():
    '''
    Gets the DBGAP dictionary. Maps a dbgap ID to a list of genes.
    Key: DBGAP ID -> str
    Value: list of ENSMUSG IDs -> list(str)
    '''
    def get_ensg_to_ensmusg_dct():
        '''
        Gets a dictionary mapping ENSG ID's to their mouse homologs.
        Key: ENSG ID -> str
        Value: list of ENSMUSG IDs -> list(str)
        '''
        ensg_to_ensmusg_dct = {}
        f = open('./data/mouse_data/mart_export.txt', 'r')
        for i, line in enumerate(f):
            # Skip header.
            if i == 0:
                continue
            ensg_id, ensmusg_id = line.split()
            if ensg_id not in ensg_to_ensmusg_dct:
                ensg_to_ensmusg_dct[ensg_id] = []
            ensg_to_ensmusg_dct[ensg_id] += [ensmusg_id]
        f.close()
        return ensg_to_ensmusg_dct

    if data_type == 'mouse':
        ensg_to_ensmusg_dct = get_ensg_to_ensmusg_dct()

    dbgap_to_gene_dct = {}
    f = open('./data/dbgap.txt', 'r')
    for line in f:
        dbgap_id, ensg_id, bloat_1, bloat_2 = line.split()

        # ENSG values are single genes.
        value = [ensg_id]
        if data_type == 'mouse':
            # Convert human to mouse homolog list.
            value = 'null'
            if ensg_id in ensg_to_ensmusg_dct:
                value = ensg_to_ensmusg_dct[ensg_id]
        if value == 'null':
            continue

        if dbgap_id not in dbgap_to_gene_dct:
            dbgap_to_gene_dct[dbgap_id] = []
        dbgap_to_gene_dct[dbgap_id] += value
    f.close()
    return dbgap_to_gene_dct

def get_sorted_fisher_dct(clus_genes, label_dct):
    '''
    Returns a dictionary of labels and their p-values of Fisher's tests to the
    input cluster genes.
    Key: Label term -> str
    Value: p-values of Fisher's test -> float
    '''
    # TODO: Change the cluster genes.
    # clus_genes = gene_universe.intersection(clus_genes)
    fisher_dct = {}
    for label in label_dct:
        # labeled_genes = set(label_dct[label])
        labeled_genes = gene_universe.intersection(label_dct[label])

        # Compute the four sets for Fisher's test.
        clus_and_label = len(clus_genes.intersection(labeled_genes))
        clus_not_label = len(clus_genes.difference(labeled_genes))
        label_not_clus = len(labeled_genes.difference(clus_genes))
        neither = len(gene_universe) - len(labeled_genes.union(clus_genes))

        # Run Fisher's test.
        f_table = ([[clus_and_label, clus_not_label], [label_not_clus,
            neither]])
        o_r, p_value = fisher_exact(f_table, alternative='greater')
        # Handle overflow issues.
        fisher_dct[label] = max(p_value, 1e-300)

    return sorted(fisher_dct.items(), key=operator.itemgetter(1))

def compute_label_enrichments(in_fname, out_fname, label_dct):
    '''
    Takes in a cluster name, reads it, and writes to the out_fname.
    '''
    cluster_dct = file_operations.get_cluster_dictionary(in_fname)
    out = open(out_fname, 'w')
    # Loop through the clusters.
    for clus_id in cluster_dct:
        clus_genes = set(cluster_dct[clus_id])

        sorted_fisher_dct = get_sorted_fisher_dct(clus_genes, label_dct)
        
        # Get the log of the top 5 enrichment p-values.
        top_label_terms, top_p_values = [], []
        for (label_term, p_value) in sorted_fisher_dct[:5]:
            top_label_terms += [label_term]
            top_p_values += [str(p_value)]
        out.write('Cluster %s\n%s\n%s\n' % (clus_id, '\t'.join(top_label_terms),
            '\t'.join(top_p_values)))
    out.close()

def main():
    if len(sys.argv) != 5:
        print ('Usage:python %s data_type objective_function run_num '
            'label_type' % sys.argv[0])
        exit()
    global data_type, objective_function, run_num, gene_universe, label_type
    data_type, objective_function, run_num, label_type = sys.argv[1:]
    assert objective_function in ['oclode', 'schaeffer', 'wlogv', 'prosnet']
    assert run_num.isdigit()
    assert label_type in ['go', 'dbgap']

    if 'prosnet_' in data_type:
        base_data_type = data_type.split('_')[1]
        if base_data_type.isdigit():
            base_data_type = file_operations.get_tcga_disease_list()[int(
                base_data_type)]
        data_type = 'prosnet_' + base_data_type
    else:
        if data_type.isdigit():
            data_type = file_operations.get_tcga_disease_list()[int(data_type)]
        base_data_type = data_type

    # Grab the BP dictionary.
    if label_type == 'go':
        label_dct = get_bp_dct(base_data_type)
    else:
        label_dct = read_dbgap_file()

    # Gene universe is intersection of network genes and GO-labeled genes.
    network_genes = file_operations.get_high_std_genes(base_data_type)
    # bp_labeled_genes = [gene for sublist in bp_dct.values() for gene in sublist]
    # TODO: change the gene universe.
    # gene_universe = set(network_genes).intersection(bp_labeled_genes)
    gene_universe = set(network_genes)

    results_folder = './results/%s_results/%s' % (data_type, objective_function)    
    for network in ['go', 'no_go']:
        if network == 'go':
            cluster_fname = '%s/clusters_%s/clusters_%s_clean_%s.txt' % (
                results_folder, network, network, run_num)
        else:
            cluster_fname = '%s/clusters_%s/clusters_%s_%s.txt' % (
                results_folder, network, network, run_num)
                    
        subfolder = '%s/%s_enrichment_terms_%s' % (results_folder,
            label_type, network)
        if not os.path.exists(subfolder):
            os.makedirs(subfolder)

        out_fname = ('%s/%s_enrichment_terms_%s_%s.txt' % (
            subfolder, label_type, network, run_num))
        compute_label_enrichments(cluster_fname, out_fname, label_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))