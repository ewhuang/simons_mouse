### Author: Edward Huang

import math
import sys
import matplotlib
import numpy as np

### This script plots a scatterplot for four networks. Network without GO,
### network with GO for lambda=0.9, WGCNA clusters, and preliminary clusters
### from Sheng's prediction networks.

matplotlib.use('Agg')
THRESHOLD_MULT = 0.8
import pylab

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: %s go_domain run_num' % sys.argv[0]
        exit()
    go_domain = sys.argv[1]
    go_domain_list = ['bp', 'cc', 'mf']
    assert go_domain in go_domain_list
    go_domain_index = go_domain_list.index(go_domain)
    go_run_num = sys.argv[2]

    for mode in ['wgcna', 'go', 'no_go']:
        pts = []
        if mode == 'wgcna':
            f = open('./wgcna/genes_only_results/clus_info_genes_only_%s.txt' %(
                go_domain), 'r')
        elif mode == 'no_go':
            f = open('./results/clus_info_no_go/clus_info_no_go_46.txt', 'r')
        elif mode == 'go':
            f = open('./results/clus_info_go/clus_info_go_%s_%d.txt' % (
                go_run_num, go_domain_index), 'r')
        # elif mode == 'ppi_go':
        #     f = open('./results/clus_info_go_35.txt', 'r')
        # elif mode == 'ppi_no_go':
        #     f = open('./results/clus_info_no_go_36.txt', 'r')
        # elif mode == 'high_std_no_go':
        #     f = open('./results/clus_info_no_go_42.txt', 'r')
        for i, line in enumerate(f):
            if i < 3:
                continue
            line = line.split()
            # Skip the clusters with fewer than 50 genes.
            if float(line[4]) < 10:
                continue
            weighted_in_out_ratio = math.log(float(line[3]), math.e)
            # weighted_in_out_ratio = float(line[3])
            in_dens = float(line[1])
            out_dens = float(line[2])
            in_out_ratio = in_dens / (in_dens + out_dens)
            top_enrichment_p = -math.log(float(line[8]), 10)
            pts += [(top_enrichment_p, weighted_in_out_ratio)]
        f.close()

        if mode == 'wgcna':
            matplotlib.pyplot.scatter(*zip(*pts), color='blue', label=mode)
            med_in_dens = math.log(np.median([pow(math.e, pt[1]) for pt in pts]
                ) * THRESHOLD_MULT, math.e)
            matplotlib.pyplot.axhline(med_in_dens)
        elif mode == 'no_go':
            matplotlib.pyplot.scatter(*zip(*pts), color='purple', label=mode)
        elif mode == 'go':
            matplotlib.pyplot.scatter(*zip(*pts), color='red', label=mode)
        # elif mode == 'ppi_go':
        #     matplotlib.pyplot.scatter(*zip(*pts), color='black', label=mode)
        # elif mode == 'ppi_no_go':
        #     matplotlib.pyplot.scatter(*zip(*pts), color='orange', label=mode)
        # elif mode == 'high_std_no_go':
        #     matplotlib.pyplot.scatter(*zip(*pts), color='green', label=mode)

    matplotlib.pyplot.axvline(5)


    matplotlib.pyplot.title('Cluster in-densities vs. Best enrichment p-values')
    matplotlib.pyplot.xlabel('negative log of lowest GO enrichment p-value')
    matplotlib.pyplot.ylabel('log of in-density/out-density ratio')
    matplotlib.pyplot.legend(loc='upper right')
    matplotlib.pyplot.show()
    pylab.savefig('./results/best_plots.png')