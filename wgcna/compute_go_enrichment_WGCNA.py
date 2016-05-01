### Author: Edward Huang

import json
import math
import operator
from scipy.stats import fisher_exact
import sys
import time

### This script calculates the GO enrichment score for each cluster using
### Fisher's exact test. Writes out a file, where each cluster shows the five
### GO terms that are the most highly enriched for each cluster, followed by
### their p-values.
### Run time: 52 seconds.

def merge_GO_enrichments():
    '''
    Merges the three GO dictionaries corresponding to each GO domain.
    '''
    # First, load all of the GO dictionaries.
    with open('../data/biological_process.json', 'r') as fp:
        bp_go_gene_dct = json.load(fp)
    fp.close()

    with open('../data/cellular_component_go.json', 'r') as fp:
        cc_go_gene_dct = json.load(fp)
    fp.close()

    with open('../data/molecular_function.json', 'r') as fp:
        mf_go_gene_dct = json.load(fp)
    fp.close()

    bp_go_gene_dct.update(cc_go_gene_dct)
    bp_go_gene_dct.update(mf_go_gene_dct)
    return bp_go_gene_dct

def get_cluster_dictionary():
    '''
    Returns a dictionary, keys=cluster ID's, values=lists of genes in the
    corresponding clusters.
    '''
    cluster_wgcna_dct = {}
    f = open('./results/WGCNA_clusters_high_std_genes.txt', 'r')
    # Read in the cluster file to create the cluster dictionary.
    for i, line in enumerate(f):
        if i == 0:
            continue
        newline = line.strip().split('\t')
        cluster = newline[2][len('Cluster '):]
        # Skip trashcan clusters.
        if cluster == '0':
            continue
        gene = newline[1][len('Gene '):]
        if cluster not in cluster_wgcna_dct:
            cluster_wgcna_dct[cluster] = [gene]
        else:
            cluster_wgcna_dct[cluster] += [gene]
    f.close()
    return cluster_wgcna_dct

def get_high_std_genes():
    '''
    Retrieves the list of genes that have high standard deviations across their
    gene expression vectors.
    '''
    high_std_genes = []
    f = open('../data/high_std_genes.txt', 'r')
    for line in f:
        gene = line.strip()        
        assert 'ENSMUSG' in gene
        high_std_genes += [gene]
    f.close()
    return high_std_genes

def get_sorted_fisher_dct(clus_genes, go_dct, gene_universe):
    '''
    Returns a dictionary, where keys=GO terms, values=p-values of Fisher's test
    GO enrichment for the set of genes in the cluster.
    '''

    # Keys are GO terms, and values are the enrichment p-values.
    fisher_dct = {}
    for go_label in go_dct:
        go_genes = set(go_dct[go_label])

        # Skip bad GO terms.
        if len(go_genes) > 1000 or len(go_genes) < 10:
            continue

        # Compute the four sets for Fisher's test.
        clus_and_go = len(clus_genes.intersection(go_genes))
        clus_not_go = len(clus_genes.difference(go_genes))
        go_not_clus = len(go_genes.difference(clus_genes))
        neither = len(gene_universe) - len(go_genes.union(clus_genes))

        # Run Fisher's test.
        f_table = ([[clus_and_go, clus_not_go], [go_not_clus, neither]])
        o_r, p_value = fisher_exact(f_table)

        fisher_dct[go_label] = p_value

    return sorted(fisher_dct.items(), key=operator.itemgetter(1))

def compute_go_enrichments(go_dct, cluster_wgcna_dct):
    gene_universe = get_high_std_genes()

    out = open('./results/cluster_enrichment_terms_wgcna.txt', 'w')

    # Loop through the clusters.
    for i in range(len(cluster_wgcna_dct)):
        clus_id = str(i + 1)
        clus_genes = set(cluster_wgcna_dct[clus_id])

        sorted_fisher_dct = get_sorted_fisher_dct(clus_genes, go_dct,
            gene_universe)

        top_go_terms, top_p_values = [], []
        for (go_term, p_value) in sorted_fisher_dct[:5]:
            top_go_terms += [go_term]
            top_p_values += [str(p_value)]
        out.write('Cluster %s\n' % clus_id)
        out.write('\t'.join(top_go_terms) + '\n')
        out.write('\t'.join(top_p_values) + '\n')
    out.close()

def main():
    go_dct = merge_GO_enrichments()
    cluster_wgcna_dct = get_cluster_dictionary()

    compute_go_enrichments(go_dct, cluster_wgcna_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))