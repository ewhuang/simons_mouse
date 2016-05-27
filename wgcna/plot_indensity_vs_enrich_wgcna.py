### Author: Edward Huang

import math
import sys
import matplotlib
import numpy as np
import sys

### This script plots a scatterplot for both networks with and without GO
### the in-density of each cluster against the negative log of its best
### enrichment p-value.

matplotlib.use('Agg')
THRESHOLD_MULT = 0.8
import pylab

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: %s genes_only/pca/mean/median mf/cc/bp' % sys.argv[0]
        exit()
    go_method = sys.argv[1]
    assert go_method in ['genes_only', 'pca', 'mean', 'median']
    domain = sys.argv[2]
    assert domain in ['mf', 'cc', 'bp']

    domain_pts = []
    f = open('./%s_results/clus_info_%s_%s.txt' % (go_method, go_method,
        domain), 'r')            
    for i, line in enumerate(f):
        if i < 3:
            continue
        line = line.split()
        if float(line[4]) < 50:
            continue
        in_dens = float(line[1])
        in_out_ratio = math.log(float(line[3]), 10)
        top_enrichment_p = -math.log(float(line[8]), 10)

        domain_pts += [(top_enrichment_p, in_out_ratio)]
    f.close()

    gene_pts = []
    f = open('./genes_only_results/clus_info_genes_only_%s.txt' % domain,
        'r')            
    for i, line in enumerate(f):
        if i < 3:
            continue
        line = line.split()
        if float(line[4]) < 50:
            continue
        in_dens = float(line[1])
        in_out_ratio = math.log(float(line[3]), 10)
        top_enrichment_p = -math.log(float(line[8]), 10)

        gene_pts += [(top_enrichment_p, in_out_ratio)]
    f.close()

    matplotlib.pyplot.scatter(*zip(*domain_pts), color='red', label=domain)
    med_in_dens = np.median([pt[1] for pt in gene_pts]) * THRESHOLD_MULT
    matplotlib.pyplot.axhline(med_in_dens)

    matplotlib.pyplot.scatter(*zip(*gene_pts), color='black', label='no_go')

    matplotlib.pyplot.title('In-Density/Out-Density Ratio vs. Best GO enric'
        'hment, WGCNA')
    matplotlib.pyplot.xlabel('negative log of lowest GO enrichment p-value')
    matplotlib.pyplot.ylabel('in-density/out-density ratio')
    matplotlib.pyplot.legend(loc='upper right')
    matplotlib.pyplot.legend(loc='upper right')
    matplotlib.pyplot.xlim(0, 50)
    matplotlib.pyplot.ylim(0, 5)
    matplotlib.pyplot.show()
    pylab.savefig('./%s_results/wgcna_cluster_plots_%s_%s.png' % (go_method,
        go_method,domain))