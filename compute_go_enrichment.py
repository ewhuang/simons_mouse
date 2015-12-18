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

    print 'Creating cluster dictionary...'
    cluster_go_dct = {}
    f = open('./results/clusters_go_clean_1.txt', 'r')
    # Read in the cluster file to create the cluster dictionary.
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.strip().split()
        clus_id = line[5]
        gene = line[3]
        all_genes.add(gene)
        if clus_id not in cluster_go_dct:
            cluster_go_dct[clus_id] = [gene]
        else:
            cluster_go_dct[clus_id] += [gene]
    f.close()

    cluster_no_go_dct = {}
    f = open('./results/clusters_no_go_1.txt', 'r')
    # Read in the cluster file to create the cluster dictionary.
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.strip().split()
        clus_id = line[5]
        gene = line[3]
        all_genes.add(gene)
        if clus_id not in cluster_no_go_dct:
            cluster_no_go_dct[clus_id] = [gene]
        else:
            cluster_no_go_dct[clus_id] += [gene]
    f.close()

    # Find GO enrichment for each cluster.
    go_p_vals = []
    go_top_labels = {}
    out = open('./results/go_top_go_1.txt', 'w')
    for i in range(len(cluster_go_dct)):
        clus_id = str(i + 1)
        fisher_dct = {}
        clus_genes = set(cluster_go_dct[clus_id])
        for go_label in go_dct:
            go_genes = set(go_dct[go_label])
            if len(go_genes) > 1000 or len(go_genes) < 10:
                continue
            clus_and_go = len(clus_genes.intersection(go_genes))
            clus_not_go = len(clus_genes.difference(go_genes))
            go_not_clus = len(go_genes.difference(clus_genes))
            neither = len(all_genes) - len(go_genes.union(clus_genes))
            o_r, p_value = fisher_exact([[clus_and_go, clus_not_go], [go_not_clus, neither]])
            fisher_dct[go_label] = p_value
        top_go = sorted(fisher_dct.items(), key=operator.itemgetter(1))
        go_p_vals += [math.log(x[1], 10) for x in top_go[:5]]
        for (label, p_value) in top_go[:5]:
            out.write(str(p_value) + '\t')
            if label not in go_top_labels:
                go_top_labels[label] = 1
            else:
                go_top_labels[label] += 1
        out.write('\n')
    out.close()

    no_go_p_vals = []
    no_go_top_labels = {}
    out = open('./results/go_top_no_go_1.txt', 'w')
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
            o_r, p_value = fisher_exact([[clus_and_go, clus_not_go], [go_not_clus, neither]])
            fisher_dct[go_label] = p_value
        top_go = sorted(fisher_dct.items(), key=operator.itemgetter(1))
        no_go_p_vals += [math.log(x[1], 10) for x in top_go[:5]]
        for label in sorted(top_go[:5]):
            label = label[0]
            out.write(label + '\t')
            if label not in no_go_top_labels:
                no_go_top_labels[label] = 1
            else:
                no_go_top_labels[label] += 1
        out.write('\n')
    out.close()
        
    bins = numpy.linspace(-25, 0, 100)
    plt.hist(go_p_vals, bins, alpha=0.5, label='GO')
    plt.hist(no_go_p_vals, bins, alpha=0.5, label='no GO')

    plt.xlabel('log of p values')
    plt.ylabel('number of clusters')
    plt.title('GO Enrichment p-values for top 5 labels of each cluster')
    plt.legend(loc='upper left')
    plt.show()