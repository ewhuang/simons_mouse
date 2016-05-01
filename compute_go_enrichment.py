### Author: Edward Huang

import file_operations
import math
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

MAX_GO_SIZE = 1000
MIN_GO_SIZE = 10

gene_universe = file_operations.get_high_std_genes()

def get_sorted_fisher_dct(clus_genes, go_dct):
    '''
    Returns a dictionary, where keys=GO terms, values=p-values of Fisher's test
    GO enrichment for the set of genes in the cluster.
    '''

    # Keys are GO terms, and values are the enrichment p-values.
    fisher_dct = {}
    for go_label in go_dct:
        go_genes = set(go_dct[go_label])

        # Skip bad GO terms.
        if len(go_genes) > MAX_GO_SIZE or len(go_genes) < MIN_GO_SIZE:
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
    for i in range(len(cluster_dct)):
        clus_id = str(i + 1)
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
    compute_go_enrichments(out_fname, cluster_dct)

def read_go_dictionaries():
    with open('./data/biological_process.json', 'r') as fp:
        bp_go_gene_dct = json.load(fp)
    fp.close()

    with open('./data/cellular_component_go.json', 'r') as fp:
        cc_go_gene_dct = json.load(fp)
    fp.close()

    with open('./data/molecular_function.json', 'r') as fp:
        mf_go_gene_dct = json.load(fp)
    fp.close()

    return (bp_go_gene_dct, cc_go_gene_dct, mf_go_gene_dct)

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s run_num' % sys.argv[0]
        exit()
    run_num = sys.argv[1]

    bp_go_gene_dct, cc_go_gene_dct, mf_go_gene_dct = read_go_dictionaries()

    # Compute GO enrichment for networks with GO.
    go_dct_list = [bp_go_gene_dct, cc_go_gene_dct, mf_go_gene_dct]

    for domain_index in range(len(go_dct_list)):
        go_dct = go_dct_list[domain_index]
        cluster_fname = './results/clusters_go/clusters_go_clean_%s_%d.txt' % (
            run_num, domain_index)
        go_fname = './results/cluster_enrichment_terms_go/cluster_enrichment_terms_go_%s_%d.txt' % (
            run_num, domain_index)
        write_enrichment_files(cluster_fname, go_fname, go_dct)

    # Compute GO enrichment for networks without GO.

    # The GO dictionary for the network without GO contains all three domains.
    # Merge the three GO dictionaries.
    bp_go_gene_dct.update(cc_go_gene_dct)
    bp_go_gene_dct.update(mf_go_gene_dct)

    no_go_cluster_fname = './results/clusters_no_go/clusters_no_go_%s.txt' % (
        run_num)
    no_go_fname = './results/cluster_enrichment_terms_no_go/cluster_enrichment_terms_no_go_%s.txt' % run_num
    write_enrichment_files(no_go_cluster_fname, no_go_fname, bp_go_gene_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))