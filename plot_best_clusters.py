### Author: Edward Huang

import file_operations
import math
import sys
import matplotlib
import numpy as np

### This script plots a scatterplot for four networks. Network without GO,
### network with GO for lambda=0.9, WGCNA clusters, and preliminary clusters
### from Sheng's prediction networks.

matplotlib.use('Agg')
THRESHOLD_MULT = 0.8
import matplotlib.pyplot as plt
import pylab

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: %s data_type run_num' % sys.argv[0]
        exit()
    data_type, run_num = sys.argv[1:]
    assert data_type == 'mouse' or data_type.isdigit()
    assert run_num.isdigit()

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    # 8 and 9 are the indices of GO and DBGAP enrichments, respectively.
    for dbgap_or_go in [8, 9]:
        highest_p = 0
        colors = ['blue', 'purple', 'red', 'black', 'orange', 'green', 'yellow']
        for mode_index, mode in enumerate(['wgcna', 'prosnet_go',
            'prosnet_no_go', 'wlogv_no_go', 'wlogv_go']):
            pts = []

            no_go_fname = 'clus_info_no_go/clus_info_no_go_%s.txt' % run_num
            go_fname = 'clus_info_go/clus_info_go_%s.txt' % run_num

            if mode == 'wgcna':      
                f = open('./wgcna/results/%s_results/genes_only/clus_info_'
                            'genes_only_bp.txt' % data_type, 'r')
            elif mode == 'wlogv_go':
                f = open('./results/%s_results/wlogv/' % data_type + go_fname,
                    'r')
                cheating_pts = []
            elif mode == 'wlogv_no_go':
                f = open('./results/%s_results/wlogv/' % data_type +
                    no_go_fname, 'r')                
            elif mode == 'prosnet_go':
                f = open('./results/prosnet_%s_results/wlogv/' % data_type +
                    go_fname, 'r')
            elif mode == 'prosnet_no_go':
                f = open('./results/prosnet_%s_results/wlogv/' % data_type +
                    no_go_fname, 'r')

            for i, line in enumerate(f):
                if i < 3:
                    continue
                line = line.split()
                # Skip the clusters with fewer than 50 genes.
                if float(line[4]) < 10:
                    continue
                ratio = float(line[3]) # For wlogv_go, ratio is a p-value.
                in_dens = float(line[1])
                out_dens = float(line[2])
                in_out_ratio = in_dens / (in_dens + out_dens)
                top_enrichment_p = -math.log(float(line[dbgap_or_go]), 10)
                if mode == 'wlogv_go':
                    if ratio < 1e-5:
                        cheating_pts += [(top_enrichment_p, in_out_ratio)]
                    else:
                        pts += [(top_enrichment_p, in_out_ratio)]
                else:
                    pts += [(top_enrichment_p, in_out_ratio)]
            f.close()
            highest_p = max(highest_p, max(pt[0] for pt in pts))

            plt.scatter(*zip(*pts), color=colors[mode_index], label=mode)
            if mode == 'wlogv_go':
                if len(cheating_pts) == 0:
                    continue
                plt.scatter(*zip(*cheating_pts), color=colors[mode_index],
                    label=mode, marker='>')                

            # WGCNA is our baseline, so get the median of their clusters' in-
            # densities and plot it as a horizontal line.
            if mode == 'wgcna':
                med_in_dens = np.median([pt[1] for pt in pts]) * THRESHOLD_MULT
                plt.axhline(med_in_dens)

        plt.axvline(10)
        plt.title('Cluster in-densities vs. Best enrichment p-values')
        plt.xlabel('negative log of lowest GO enrichment p-value')
        plt.ylabel('in-density/(in-density + out-density)')
        plt.legend(loc='lower right')
        plt.ylim(0, 1.2)
        plt.xlim(0, highest_p * 1.2)
        plt.show()
        if dbgap_or_go == 8:
            pylab.savefig('./results/%s_results/comparison_plots/go_comparison_'
                            'plot_%s.png' % (data_type, run_num))
        else:
            pylab.savefig('./results/%s_results/comparison_plots/dbgap_'
                            'comparison_plot_%s.png' % (data_type, run_num))
        plt.close()