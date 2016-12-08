### Author: Edward Huang

import file_operations
import json
import operator
import os
from scipy.stats import fisher_exact
import sys
import time

### This script finds the enrichment with DBGAP gene annotations.
### Obtain mart_export.txt by going to ensembl.org/biomart
### Dataset -> Homo sapiens genes (GRCh38.p7)
### Filters -> Orthologous Mouse Genes: Only
### Attributes -> Ensembl Gene ID, Mouse Ensembl Gene ID
### Export as TSV, tick "unique results only".

def get_sorted_fisher_dct(clus_genes):
    '''
    Returns a dictionary, where keys=GO terms, values=p-values of Fisher's test
    GO enrichment for the set of genes in the cluster.
    '''
    # Keys are GO terms, and values are the enrichment p-values.
    fisher_dct = {}
    for dbgap_label in dbgap_dct:
        dbgap_genes = set(dbgap_dct[dbgap_label])

        # Compute the four sets for Fisher's test.
        clus_and_dbgap = len(clus_genes.intersection(dbgap_genes))
        clus_not_dbgap = len(clus_genes.difference(dbgap_genes))
        dbgap_not_clus = len(dbgap_genes.difference(clus_genes))
        neither = len(gene_universe) - len(dbgap_genes.union(clus_genes))

        # Run Fisher's test.
        f_table = ([[clus_and_dbgap, clus_not_dbgap], [dbgap_not_clus,
            neither]])
        o_r, p_value = fisher_exact(f_table, alternative='greater')

        # Handle overflow issues.
        p_value = max(p_value, 1e-300)

        fisher_dct[(dbgap_label, clus_and_dbgap, clus_not_dbgap, dbgap_not_clus,
            neither)] = p_value

    return sorted(fisher_dct.items(), key=operator.itemgetter(1))

def compute_dbgap_enrichments(fname, cluster_dct, side_fname):
    '''
    Takes in a filename, a cluster dictionary, and a GO dictionary. Writes out
    to the filename the top 5 GO enrichments for each cluster.
    '''
    out = open(fname, 'w')
    side =open(side_fname, 'w')
    # Loop through the clusters.
    for clus_id in cluster_dct:
        clus_genes = set(cluster_dct[clus_id])

        sorted_fisher_dct = get_sorted_fisher_dct(clus_genes)
        
        # Get the log of the top 5 enrichment p-values.
        top_dbgap_terms, top_p_values = [], []
        for ((dbgap_term, clus_and_dbgap, clus_not_dbgap, dbgap_not_clus,
            neither), p_value) in sorted_fisher_dct[:5]:
            top_dbgap_terms += [dbgap_term]
            top_p_values += [str(p_value)]
            side.write('%s\t%s\t%d\t%d\t%d\t%d\t%g\n' % (clus_id, dbgap_term,
                clus_and_dbgap, clus_not_dbgap, dbgap_not_clus, neither,
                p_value))
        out.write('Cluster %s\n' % clus_id)
        out.write('\t'.join(top_dbgap_terms) + '\n')
        out.write('\t'.join(top_p_values) + '\n')
    side.close()
    out.close()

def write_enrichment_files(in_fname, out_fname, side_fname):
    '''
    Takes in a cluster name, reads it, and writes to the out_fname.
    '''
    # Compute GO enrichment for networks without GO.
    cluster_dct = file_operations.get_cluster_dictionary(in_fname)
    compute_dbgap_enrichments(out_fname, cluster_dct, side_fname)

def main():
    if len(sys.argv) != 4:
        print 'Usage:python %s data_type objective_function run_num' % (
            sys.argv[0])
        exit()
    global data_type, objective_function, run_num
    data_type, objective_function, run_num = sys.argv[1:]
    assert objective_function in ['oclode', 'schaeffer', 'wlogv', 'prosnet']
    assert run_num.isdigit()

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

    global gene_universe, dbgap_dct
    gene_universe = file_operations.get_high_std_genes(base_data_type)

    if 'mouse' in data_type:
        dbgap_dct = file_operations.read_dbgap_file()
    else:
        dbgap_dct = file_operations.read_ensg_dbgap_file(base_data_type)

    # No GO network.

    results_folder = './results/%s_results/%s' % (data_type, objective_function)
    for network in ['go', 'no_go']:
        if network == 'go':
            cluster_fname = '%s/clusters_%s/clusters_%s_clean_%s.txt' % (
                results_folder, network, network, run_num)
        else:
            cluster_fname = '%s/clusters_%s/clusters_%s_%s.txt' % (
                results_folder, network, network, run_num)
        subfolder = '%s/dbgap_enrichment_terms_%s' % (results_folder, network)
        if not os.path.exists(subfolder):
            os.makedirs(subfolder)
        out_fname = '%s/dbgap_enrichment_terms_%s_%s.txt' % (subfolder, network,
            run_num)
        side_fname = '%s/dbgap_f_table_%s_%s.txt' % (subfolder, network, run_num)
        write_enrichment_files(cluster_fname, out_fname, side_fname)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))