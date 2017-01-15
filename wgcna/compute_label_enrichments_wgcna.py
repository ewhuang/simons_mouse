### Author: Edward Huang

import file_operations
import json
import operator
from scipy.stats import fisher_exact
import sys
import time

### This script calculates the enrichment score for each cluster using Fisher's
### exact test. Writes out a file, where each cluster shows the five label terms
### that are the most highly enriched for each cluster, along with their
### p-values.
### Run time: 52 seconds.

def get_label_dct():
    '''
    Gets the annotation dictionary.
    '''
    if label_type == 'go':
        fname = 'bp_dct'
    elif label_type == 'dbgap':
        fname = 'dbgap_dct'
    elif label_type == 'gwas':
        fname = 'gwas_dct'
    with open('../data/%s_data/%s.json' % (data_type, fname), 'r') as fp:
        label_dct = json.load(fp)
    fp.close()
    return label_dct

def get_cluster_dictionary():
    '''
    Returns a dictionary of clusters.
    Key: cluster ID -> str
    Value: lists of genes in the cluster-> list(str)
    '''
    cluster_wgcna_dct = {}
    f = open('./results/%s_results/clusters.txt' % data_type, 'r')
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
            cluster_wgcna_dct[cluster_id] = []
        cluster_wgcna_dct[cluster_id] += [gene]
    f.close()
    return cluster_wgcna_dct

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
        fisher_dct[label] = max(p_value, 1e-300)
    return sorted(fisher_dct.items(), key=operator.itemgetter(1))

def compute_label_enrichments(label_dct, cluster_wgcna_dct):
    '''
    Write out the label enrichments for each cluster.
    '''
    out = open('./results/%s_results/%s_enrichment_terms.txt' % (data_type,
        label_type), 'w')
    # Loop through the clusters.
    for clus_id in cluster_wgcna_dct:
        clus_genes = set(cluster_wgcna_dct[clus_id])

        sorted_fisher_dct = get_sorted_fisher_dct(clus_genes, label_dct)

        top_label_terms, top_p_values = [], []
        for (label_term, p_value) in sorted_fisher_dct[:5]:
            top_label_terms += [label_term]
            top_p_values += [str(p_value)]
        out.write('Cluster %s\n%s\n%s\n' % (clus_id, '\t'.join(top_label_terms),
            '\t'.join(top_p_values)))
    out.close()

def main():
    if len(sys.argv) != 3:
        print 'Usage: %s data_type go/dbgap/gwas' % sys.argv[0]
        exit()
    global data_type, gene_universe, label_type
    data_type, label_type = sys.argv[1:]
    assert data_type == 'mouse' or data_type.isdigit()
    assert label_type in ['go', 'dbgap', 'gwas']

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    label_dct = get_label_dct()

    network_genes = file_operations.get_high_std_genes(data_type)
    # labeled_genes = [gene for sublist in label_dct.values() for gene in sublist]
    
    # TODO: Switch up gene universe.
    # gene_universe = set(network_genes).intersection(bp_labeled_genes)
    gene_universe = set(network_genes)

    cluster_wgcna_dct = get_cluster_dictionary()
    compute_label_enrichments(label_dct, cluster_wgcna_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))