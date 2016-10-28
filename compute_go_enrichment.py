### Author: Edward Huang

import file_operations
import json
import operator
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

def get_sorted_fisher_dct(clus_genes, go_dct):
    '''
    Returns a dictionary, where keys=GO terms, values=p-values of Fisher's test
    GO enrichment for the set of genes in the cluster.
    '''
    # Keys are GO terms, and values are the enrichment p-values.
    fisher_dct = {}
    for go_label in go_dct:
        go_genes = set(go_dct[go_label])

        # Compute the four sets for Fisher's test.
        clus_and_go = len(clus_genes.intersection(go_genes))
        clus_not_go = len(clus_genes.difference(go_genes))
        go_not_clus = len(go_genes.difference(clus_genes))
        neither = len(gene_universe) - len(go_genes.union(clus_genes))

        # Run Fisher's test.
        f_table = ([[clus_and_go, clus_not_go], [go_not_clus, neither]])
        o_r, p_value = fisher_exact(f_table, alternative='greater')

        # Handle overflow issues.
        p_value = max(p_value, 1e-300)

        fisher_dct[go_label] = p_value

    return sorted(fisher_dct.items(), key=operator.itemgetter(1))

def compute_go_enrichments(fname, cluster_dct, go_dct):
    '''
    Takes in a filename, a cluster dictionary, and a GO dictionary. Writes out
    to the filename the top 5 GO enrichments for each cluster.
    '''
    out = open(fname, 'w')

    # Loop through the clusters.
    for clus_id in cluster_dct:
        clus_genes = set(cluster_dct[clus_id])

        sorted_fisher_dct = get_sorted_fisher_dct(clus_genes, go_dct)
        
        # Get the log of the top 5 enrichment p-values.
        top_go_terms, top_p_values = [], []
        for (go_term, p_value) in sorted_fisher_dct[:5]:
            top_go_terms += [go_term]
            top_p_values += [str(p_value)]
        out.write('Cluster %s\n' % clus_id)
        out.write('\t'.join(top_go_terms) + '\n')
        out.write('\t'.join(top_p_values) + '\n')
    out.close()

def write_enrichment_files(in_fname, out_fname, go_dct):
    '''
    Takes in a cluster name, reads it, and writes to the out_fname.
    '''
    # Compute GO enrichment for networks without GO.
    cluster_dct = file_operations.get_cluster_dictionary(in_fname)
    compute_go_enrichments(out_fname, cluster_dct, go_dct)

def get_bp_dct(base_data_type):
    '''
    Fetch the BP dictionary, and remove GO terms that overlap too much with MF.
    '''
    with open('./data/%s_data/bp_dct.json' % base_data_type, 'r') as fp:
        bp_dct = json.load(fp)
    fp.close()

    # Get the overlapping MF and BP terms.
    overlap_list = file_operations.read_go_overlap(base_data_type)

    # Remove the overlapping BP terms.
    overlapping_go_terms = set([tup[0] for tup in overlap_list])
    for overlapping_go in overlapping_go_terms:
        del bp_dct[overlapping_go]

    return bp_dct

def main():
    if len(sys.argv) != 4:
        print 'Usage:python %s data_type objective_function run_num' % sys.argv[0]
        exit()
    global data_type, objective_function, run_num
    data_type, objective_function, run_num = sys.argv[1:]
    assert data_type in ['mouse', 'prosnet_mouse'] or data_type.isdigit()
    assert objective_function in ['oclode', 'schaeffer', 'wlogv', 'prosnet']
    assert run_num.isdigit()

    if 'prosnet_' in data_type:
        base_data_type = data_type[len('prosnet_'):]
    else:
        base_data_type = data_type

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    global gene_universe
    gene_universe = file_operations.get_high_std_genes(base_data_type)

    bp_dct = get_bp_dct(base_data_type)

    results_folder = './results/%s_results/%s' % (data_type, objective_function)
    
    for network in ['go', 'no_go']:
        if network == 'go':
            cluster_fname = '%s/clusters_%s/clusters_%s_clean_%s.txt' % (
                results_folder, network, network, run_num)
        else:
            cluster_fname = '%s/clusters_%s/clusters_%s_%s.txt' % (
                results_folder, network, network, run_num)
        out_fname = ('%s/cluster_enrichment_terms_%s/cluster_enrichment_terms_'
                        '%s_%s.txt' % (results_folder, network, network,
                            run_num))
        write_enrichment_files(cluster_fname, out_fname, bp_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))