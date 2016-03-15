### Author: Edward Huang

from scipy.stats import fisher_exact
import operator
import sys
import matplotlib.pyplot as plt
import math
import numpy

### This script calculates the GO enrichment score for each clustering
### result. For each method, for every cluster, for every possible GO
### label's related genes, compute Fisher's test. Thus, for every cluster,
### obtain the lowest possible p-value corresponding to that cluster's
### most related GO label. These p-values (GO enrichments) should improve
### as a whole going from coexpression network to GO network without
### sacrificing in-density.

if __name__ == '__main__':
    all_genes = set([])
    print 'Extracting GO labels...'
    go_dct = {}
    f = open('./data/go_edges.txt', 'r')
    for line in f:
        gene, go_label = line.split()
        all_genes.add(gene)
        if go_label not in go_dct:
            go_dct[go_label] = [gene]
        else:
            go_dct[go_label] += [gene]
    f.close()

    print 'Creating cluster dictionary from WGCNA results...'
    cluster_wgcna_dct = {}
    f = open('./results/WGCNA_results/WGCNA_clusters_all_genes.txt', 'r')
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
        all_genes.add(gene)
        if cluster not in cluster_wgcna_dct:
            cluster_wgcna_dct[cluster] = [gene]
        else:
            cluster_wgcna_dct[cluster] += [gene]
    f.close()

    cluster_no_go_dct = {}
    f = open('./results/clusters_no_go_20.txt', 'r')
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
        all_genes.add(gene)
        if cluster not in cluster_no_go_dct:
            cluster_no_go_dct[cluster] = [gene]
        else:
            cluster_no_go_dct[cluster] += [gene]
    f.close()

    wgcna_p_vals = []
    wgcna_top_labels = {}
    out = open('./results/WGCNA_results/cluster_enrichment_terms_wgcna.txt', 'w')
    for i in range(len(cluster_wgcna_dct)):
        clus_id = str(i + 1)
        fisher_dct = {}
        clus_genes = set(cluster_wgcna_dct[clus_id])
        for go_label in go_dct:
            go_genes = set(go_dct[go_label])
            if len(go_genes) > 1000 or len(go_genes) < 10:
                continue
            clus_and_go = len(clus_genes.intersection(go_genes))
            clus_not_go = len(clus_genes.difference(go_genes))
            go_not_clus = len(go_genes.difference(clus_genes))
            neither = len(all_genes) - len(go_genes.union(clus_genes))
            o_r, p_value = fisher_exact([[clus_and_go, clus_not_go],
                [go_not_clus, neither]])
            fisher_dct[go_label] = p_value
        top_go = sorted(fisher_dct.items(), key=operator.itemgetter(1))
        wgcna_p_vals += [math.log(x[1], 10) for x in top_go[:5]]
        out_p = []
        out_go = []
        for (label, p_value) in top_go[:5]:
            out_p += [str(p_value)]
            out_go += [label]
            if label not in wgcna_top_labels:
                wgcna_top_labels[label] = 1
            else:
                wgcna_top_labels[label] += 1
        out.write('Cluster %s\n' % clus_id)
        out.write('\t'.join(out_go) + '\n')
        out.write('\t'.join(out_p) + '\n')
    out.close()

    no_go_p_vals = []
    no_go_top_labels = {}
    for i in range(len(cluster_no_go_dct)):
        clus_id = str(i + 1)
        fisher_dct = {}
        clus_genes = set(cluster_no_go_dct[clus_id])
        for go_label in go_dct:
            go_genes = set(go_dct[go_label])
            if len(go_genes) > 1000 or len(go_genes) < 10:
                continue
            clus_and_go = len(clus_genes.intersection(go_genes))
            clus_not_go = len(clus_genes.difference(go_genes))
            go_not_clus = len(go_genes.difference(clus_genes))
            neither = len(all_genes) - len(go_genes.union(clus_genes))
            o_r, p_value = fisher_exact([[clus_and_go, clus_not_go],
                [go_not_clus, neither]])
            fisher_dct[go_label] = p_value
        top_go = sorted(fisher_dct.items(), key=operator.itemgetter(1))
        no_go_p_vals += [math.log(x[1], 10) for x in top_go[:5]]
        out_p = []
        out_go = []
        for (label, p_value) in top_go[:5]:
            out_p += [str(p_value)]
            out_go += [label]
            if label not in no_go_top_labels:
                no_go_top_labels[label] = 1
            else:
                no_go_top_labels[label] += 1

    bins = numpy.linspace(-150, 0, 100)
    plt.hist(wgcna_p_vals, bins, alpha=0.5, label='WGCNA')
    plt.hist(no_go_p_vals, bins, alpha=0.5, label='no GO')

    plt.xlabel('log of p values')
    plt.ylabel('number of enrichments')
    plt.title('GO Enrichment p-values for top 5 labels of each cluster.')
    plt.legend(loc='upper left')
    plt.show()