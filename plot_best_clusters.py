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
        print 'Usage: %s data_type run_num' % sys.argv[0]
        exit()
    data_type = sys.argv[1]
    assert data_type in ['mouse', 'tcga']
    run_num = sys.argv[2]
    assert run_num.isdigit()
    go_domain_list = ['bp', 'mf']

    for go_domain_index in [0]:
        go_domain = go_domain_list[go_domain_index]
        colors = ['blue', 'purple', 'red', 'black', 'orange', 'green', 'yellow']
        # for mode_index, mode in enumerate(['wgcna', 'wlogv_go', 'wlogv_no_go',
        #     'oclode_go', 'oclode_no_go', 'schaeffer_go', 'schaeffer_no_go']):
        for mode_index, mode in enumerate(['wgcna', 'wlogv_go', 'wlogv_no_go']):
            pts = []

            no_go_fname = 'clus_info_no_go/clus_info_no_go_%s_%d.txt' % (
                run_num, go_domain_index)
            go_fname = 'clus_info_go/clus_info_go_%s_%d.txt' % (
                run_num, go_domain_index)

            if mode == 'wgcna':
                f = open('./wgcna/%s_results/genes_only/clus_info_genes_only_%s.txt' %(
                    data_type, go_domain), 'r')
            elif mode == 'wlogv_go':
                f = open('./%s_results/wlogv/' % data_type + go_fname, 'r')
            elif mode == 'wlogv_no_go':
                f = open('./%s_results/wlogv/' % data_type + no_go_fname, 'r')
            elif mode == 'oclode_go':
                f = open('./%s_results/oclode/' % data_type + go_fname, 'r')
            elif mode == 'oclode_no_go':
                f = open('./%s_results/oclode/' % data_type + no_go_fname, 'r')
            elif mode == 'schaeffer_go':
                f = open('./%s_results/schaeffer/' % data_type + go_fname, 'r')
            elif mode == 'schaeffer_no_go':
                f = open('./%s_results/schaeffer/' % data_type + no_go_fname,
                    'r')

            for i, line in enumerate(f):
                if i < 3:
                    continue
                line = line.split()
                # Skip the clusters with fewer than 50 genes.
                if float(line[4]) < 10:
                    continue
                weighted_in_out_ratio = math.log(float(line[3]), math.e)
                in_dens = float(line[1])
                out_dens = float(line[2])
                in_out_ratio = in_dens / (in_dens + out_dens)
                top_enrichment_p = -math.log(float(line[8]), 10)
                # pts += [(top_enrichment_p, weighted_in_out_ratio)]
                pts += [(top_enrichment_p, in_out_ratio)]
            f.close()

            matplotlib.pyplot.scatter(*zip(*pts), color=colors[mode_index],
                label=mode)

            # WGCNA is our baseline, so get the median of their clusters' in-
            # densities and plot it as a horizontal line.
            if mode == 'wgcna':
                # med_in_dens = math.log(np.median([pow(math.e, pt[1]) for pt in pts]
                #     ) * THRESHOLD_MULT, math.e)
                med_in_dens = np.median([pt[1] for pt in pts]) * THRESHOLD_MULT
                matplotlib.pyplot.axhline(med_in_dens)

        matplotlib.pyplot.axvline(10)
        matplotlib.pyplot.title('Cluster in-densities vs. Best enrichment p-values')
        matplotlib.pyplot.xlabel('negative log of lowest GO enrichment p-value')
        matplotlib.pyplot.ylabel('in-density/(in-density + out-density)')
        matplotlib.pyplot.legend(loc='lower right')
        matplotlib.pyplot.ylim(0, 1.2)
        if data_type == 'mouse':
            matplotlib.pyplot.xlim(0, 50)
        elif data_type == 'tcga':
            matplotlib.pyplot.xlim(0, 150)            
        matplotlib.pyplot.show()
        pylab.savefig('./%s_results/comparison_plots/comparison_plot_%s_%s.png' % (
            data_type, go_domain, run_num))
        matplotlib.pyplot.close()