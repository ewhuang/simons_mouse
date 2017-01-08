### Author: Edward Huang

import file_operations
import json
import operator
from scipy.stats import fisher_exact
import sys
import time

### This script calculates the enrichment score for each cluster using Fisher's
### exact test. Writes out a file, where each cluster shows the five GO terms
### that are the most highly enriched for each cluster, along with their
### p-values.
### Run time: 52 seconds.

def get_bp_dct():
    '''
   Gets the GO dictionary corresponding to the GO domain.
    '''
    with open('../data/%s_data/bp_dct.json' % (data_type), 'r') as fp:
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
        f = open('../data/mouse_data/mart_export.txt', 'r')
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
    f = open('../data/dbgap.txt', 'r')
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
    Returns a dictionary of Go terms and their p-values of Fisher's tests to the
    input cluster genes.
    Key: GO term -> str
    Value: p-values of Fisher's test -> float
    '''
    # TODO: Change the cluster genes.
    # clus_genes = gene_universe.intersection(clus_genes)
    # Keys are GO terms, and values are the enrichment p-values.
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
        fisher_dct[label] = p_value

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

        top_go_terms, top_p_values = [], []
        for (go_term, p_value) in sorted_fisher_dct[:5]:
            top_go_terms += [go_term]
            top_p_values += [str(p_value)]
        out.write('Cluster %s\n%s\n%s\n' % (clus_id, '\t'.join(top_go_terms),
            '\t'.join(top_p_values)))
    out.close()

def main():
    if len(sys.argv) != 3:
        print 'Usage: %s data_type go/dbgap' % sys.argv[0]
        exit()
    global data_type, gene_universe, label_type
    data_type, label_type = sys.argv[1:]
    assert data_type == 'mouse' or data_type.isdigit()
    assert label_type in ['go', 'dbgap']

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    if label_type == 'go':
        label_dct = get_bp_dct()
    else:
        label_dct = read_dbgap_file()

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