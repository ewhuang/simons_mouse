### Author: Edward Huang

import file_operations
import itertools
import json
from multiprocessing import Pool
import operator
import os
from scipy.stats import fisher_exact
import sys
import time

### For each method, for every cluster, for every possible GO label's related
### genes, compute Fisher's test. Thus, for every cluster, obtain the lowest
### possible p-value corresponding to that cluster's most related GO label.
### Run time: 30 seconds.

def get_label_dct(base_data_type):
    '''
    Gets the annotation dictionary.
    '''
    if label_type == 'go':
        fname = 'bp_dct'
    elif label_type == 'dbgap':
        fname = 'dbgap_dct'
    elif label_type == 'gwas':
        fname = 'gwas_dct'
    elif label_type == 'kegg':
        fname = 'kegg_dct'
    elif label_type == 'ctd':
        fname = 'ctd_dct'
    with open('./data/%s_data/%s.json' % (base_data_type, fname), 'r') as fp:
        label_dct = json.load(fp)
    fp.close()

    # TODO: This block makes us evaluate on the partial BP terms.

    if label_type == 'go':
        overlap_list = file_operations.read_go_overlap(base_data_type)
        # Remove the overlapping BP terms.
        overlapping_go_terms = set([tup[0] for tup in overlap_list])
        for overlapping_go in overlapping_go_terms:
            # GO might not be in BP dictionary because we computed overlap using
            # full set of genes, but bp_dct only contains genes in the high_std.
            if overlapping_go in label_dct:
                del label_dct[overlapping_go]

    return label_dct

def compute_fisher_table(*args):
    '''
    Computes Fisher's test between cluster genes and annotations. Returns a
    p-value.
    '''
    clus_genes, label = args[0]
    labeled_genes = gene_universe.intersection(label_dct[label])

    # Compute the four sets for Fisher's test.
    clus_and_label = len(clus_genes.intersection(labeled_genes))
    clus_not_label = len(clus_genes.difference(labeled_genes))
    label_not_clus = len(labeled_genes.difference(clus_genes))
    neither = len(gene_universe) - len(labeled_genes.union(clus_genes))

    # Run Fisher's test.
    f_table = ([[clus_and_label, clus_not_label], [label_not_clus, neither]])
    o_r, p_value = fisher_exact(f_table, alternative='greater')
    # Handle overflow issues.
    return max(p_value, 1e-300)

def get_sorted_fisher_dct(clus_genes):
    '''
    Returns a dictionary of labels and their p-values of Fisher's tests to the
    input cluster genes.
    Key: Label term -> str
    Value: p-values of Fisher's test -> float
    '''
    fisher_dct = {}
    pool = Pool(processes=20)
    label_list = label_dct.keys()
    # p_values = pool.map(compute_fisher_table_star, itertools.izip(
    #     itertools.repeat(clus_genes), label_list))
    p_values = pool.map(compute_fisher_table, itertools.izip(itertools.repeat(
        clus_genes), label_list))
    for label_idx, label in enumerate(label_list):
        fisher_dct[label] = p_values[label_idx]
    pool.close()
    pool.join()
    return sorted(fisher_dct.items(), key=operator.itemgetter(1))

def compute_label_enrichments(in_fname, out_fname):
    '''
    Takes in a cluster name, reads it, and writes to the out_fname.
    '''
    cluster_dct = file_operations.get_cluster_dictionary(in_fname)
    out = open(out_fname, 'w')
    # Loop through the clusters.
    for clus_id in cluster_dct:
        clus_genes = set(cluster_dct[clus_id])
        sorted_fisher_dct = get_sorted_fisher_dct(clus_genes)

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
            'go/dbgap/gwas/kegg/ctd' % sys.argv[0])
        exit()
    global data_type, objective_function, run_num, label_type
    data_type, objective_function, run_num, label_type = sys.argv[1:]
    assert (objective_function in ['oclode', 'schaeffer', 'wlogv', 'prosnet',
        'wgcna'])
    assert run_num.isdigit()
    assert label_type in ['go', 'dbgap', 'gwas', 'kegg', 'ctd']

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

    global label_dct, gene_universe
    label_dct = get_label_dct(base_data_type)

    # Gene universe is intersection of network genes and GO-labeled genes.
    gene_universe = set(file_operations.get_high_std_genes(base_data_type))

    results_folder = './results/%s_results/%s' % (data_type, objective_function)    

    for network in ['go', 'no_go']:
        # No 'no_go' situation for WGCNA.
        if objective_function == 'wgcna' and network == 'no_go':
            continue
        if network == 'go':
            cluster_fname = '%s/clusters_%s/clusters_%s_clean_%s.txt' % (
                results_folder, network, network, run_num)
        else:
            cluster_fname = '%s/clusters_%s/clusters_%s_%s.txt' % (
                results_folder, network, network, run_num)
                    
        subfolder = '%s/%s_enrichment_terms_%s' % (results_folder, label_type,
            network)
        if not os.path.exists(subfolder):
            os.makedirs(subfolder)

        out_fname = ('%s/%s_enrichment_terms_%s_%s.txt' % (subfolder,
            label_type, network, run_num))
        compute_label_enrichments(cluster_fname, out_fname)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))