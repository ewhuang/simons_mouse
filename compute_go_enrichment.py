### Author: Edward Huang

# from scipy.stats import fisher_exact
import file_operations
import fisher_test
import math
import operator
import sys
# import matplotlib.pyplot as plt
# import numpy

### This script calculates the GO enrichment score for each clustering
### result. For each method, for every cluster, for every possible GO
### label's related genes, compute Fisher's test. Thus, for every cluster,
### obtain the lowest possible p-value corresponding to that cluster's
### most related GO label. These p-values (GO enrichments) should improve
### as a whole going from coexpression network to GO network without
### sacrificing in-density.

MAX_GO_SIZE = 1000
MIN_GO_SIZE = 10

sampled_genes = file_operations.get_sampled_genes()
go_dct = file_operations.get_go_labels(sampled_genes)

def write_out_go_enrichments(fname, cluster_dct):
    # Find GO enrichment for each cluster.
    go_p_vals = []
    go_top_labels = {}
    out = open(fname, 'w')
    for i in range(len(cluster_dct)):
        clus_id = str(i + 1)
        fisher_dct = {}
        clus_genes = set(cluster_dct[clus_id])
        for go_label in go_dct:
            go_genes = set(go_dct[go_label])
            # Skip giant or tiny GO terms.
            if len(go_genes) > MAX_GO_SIZE or len(go_genes) < MIN_GO_SIZE:
                continue
            # Compute the four disjoint set sizes of Venn diagram for Fisher's.
            clus_and_go = len(clus_genes.intersection(go_genes))
            clus_not_go = len(clus_genes.difference(go_genes))
            go_not_clus = len(go_genes.difference(clus_genes))
            neither = len(sampled_genes) - len(go_genes.union(clus_genes))
            # Compute Fisher's exact test.
            f_table = ([[clus_and_go, clus_not_go], [go_not_clus, neither]])
            p_value = fisher_test.FishersExactTest(f_table).two_tail_p()
            if p_value == 0:
                p_value = 1e-300
            fisher_dct[go_label] = p_value

        top_go = sorted(fisher_dct.items(), key=operator.itemgetter(1))
        # Get the log of the top 5 enrichment p-values.
        go_p_vals += [math.log(x[1], 10) for x in top_go[:5]]
        out_p = []
        out_go = []
        for (label, p_value) in top_go[:5]:
            out_p += [str(p_value)]
            out_go += [label]
            if label not in go_top_labels:
                go_top_labels[label] = 1
            else:
                go_top_labels[label] += 1
        out.write('Cluster %s\n' % clus_id)
        out.write('\t'.join(out_go) + '\n')
        out.write('\t'.join(out_p) + '\n')
    out.close()
    return go_p_vals, go_top_labels

if __name__ == '__main__':
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print 'Usage:python %s run_num predicted?' % sys.argv[0]
        exit()
    run_num = sys.argv[1]

    # Cluster dictionary generation.
    # Compute GO enrichment without GO nodes, so we use cleaned file.
    go_cluster_fname = './results/clusters_go_clean_%s.txt' % run_num
    cluster_go_dct = file_operations.create_cluster_dct(go_cluster_fname)
    assert len(cluster_go_dct) == 20

    no_go_cluster_fname = './results/clusters_no_go_%s.txt' % run_num
    cluster_no_go_dct = file_operations.create_cluster_dct(no_go_cluster_fname)
    assert len(cluster_no_go_dct) == 20

    go_fname = './results/cluster_enrichment_terms_go_%s.txt' % run_num
    go_p_vals, go_top_labels = write_out_go_enrichments(go_fname,
        cluster_go_dct)

    no_go_fname = './results/cluster_enrichment_terms_no_go_%s.txt' % run_num
    no_go_p_vals, no_go_top_labels = write_out_go_enrichments(no_go_fname,
        cluster_no_go_dct)

    # # Plotting GO enrichment histograms.
    # bins = numpy.linspace(-150, 0, 100)
    # plt.hist(go_p_vals, bins, alpha=0.5, label='GO')
    # plt.hist(no_go_p_vals, bins, alpha=0.5, label='no GO')

    # plt.xlabel('log of p values')
    # plt.ylabel('number of enrichments')
    # plt.title('GO Enrichment p-values for top 5 labels of each cluster')
    # plt.legend(loc='upper left')
    # plt.show()