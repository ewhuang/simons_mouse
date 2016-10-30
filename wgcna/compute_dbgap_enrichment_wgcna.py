### Author: Edward Huang

import file_operations
import json
import math
import operator
from scipy.stats import fisher_exact
import sys
import time

### This script calculates the DBGAP enrichment score for each cluster using
### Fisher's exact test. Writes out a file, where each cluster shows the five
### DBGAP terms that are the most highly enriched for each cluster, followed by
### their p-values.

def get_cluster_dictionary():
    '''
    Returns a dictionary, keys=cluster ID's, values=lists of genes in the
    corresponding clusters.
    '''
    cluster_wgcna_dct = {}
    f = open('./results/%s_results/%s/clusters_%s.txt' % (data_type,
        network_type, network_type), 'r')
    # Read in the cluster file to create the cluster dictionary.
    for i, line in enumerate(f):
        if i == 0:
            continue
        newline = line.strip().split('\t')
        cluster_id = newline[2][len('Cluster '):]
        # Skip trashcan clusters.
        if cluster_id == '0':
            continue
        gene = newline[1][len('Gene '):]
        if cluster_id not in cluster_wgcna_dct:
            cluster_wgcna_dct[cluster_id] = [gene]
        else:
            cluster_wgcna_dct[cluster_id] += [gene]
    f.close()
    return cluster_wgcna_dct

def get_high_std_genes():
    '''
    Retrieves the list of genes that have high standard deviations across their
    gene expression vectors.
    '''
    high_std_genes = []
    f = open('../data/%s_data/high_std_genes.txt' % data_type, 'r')
    for line in f:
        gene = line.strip()
        assert 'ENSMUSG' in gene or 'ENSG' in gene
        high_std_genes += [gene]
    f.close()
    return high_std_genes

def get_sorted_fisher_dct(clus_genes, gene_universe):
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

        fisher_dct[dbgap_label] = p_value

    return sorted(fisher_dct.items(), key=operator.itemgetter(1))

def compute_dbgap_enrichments(cluster_wgcna_dct):
    gene_universe = get_high_std_genes()

    out = open('./results/%s_results/%s/dbgap_enrichment_terms_%s.txt' % (
        data_type, network_type, network_type), 'w')

    # Loop through the clusters.
    for clus_id in cluster_wgcna_dct:
        clus_genes = set(cluster_wgcna_dct[clus_id])

        sorted_fisher_dct = get_sorted_fisher_dct(clus_genes, gene_universe)

        top_go_terms, top_p_values = [], []
        for (go_term, p_value) in sorted_fisher_dct[:5]:
            top_go_terms += [go_term]
            top_p_values += [str(p_value)]
        out.write('Cluster %s\n' % clus_id)
        out.write('\t'.join(top_go_terms) + '\n')
        out.write('\t'.join(top_p_values) + '\n')
    out.close()

def main():
    if len(sys.argv) != 3:
        print 'Usage: %s data_type genes_only/pca/mean/median' % sys.argv[0]
        exit()
    global data_type, network_type
    data_type = sys.argv[1]
    assert data_type == 'mouse' or data_type.isdigit()
    network_type = sys.argv[2]
    assert network_type in ['genes_only', 'pca', 'mean', 'median']

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    # Fetch the DBGAP dictionary.
    global dbgap_dct
    dbgap_dct = file_operations.read_dbgap_file()

    if network_type == 'genes_only':
        cluster_wgcna_dct = get_cluster_dictionary()
    else:
        cluster_wgcna_dct = get_cluster_dictionary()
    compute_dbgap_enrichments(cluster_wgcna_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))