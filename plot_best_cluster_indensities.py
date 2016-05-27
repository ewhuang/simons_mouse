### Author: Edward Huang

import math
import sys
import matplotlib.pyplot as plt
import numpy as np

### This script plots a scatterplot for four networks. Network without GO,
### network with GO for lambda=0.9, WGCNA clusters, and preliminary clusters
### from Sheng's prediction networks.

THRESHOLD_MULT = 0.8

if __name__ == '__main__':
    for mode in ['wgcna', 'no_go', 'go', 'ppi_go', 'ppi_no_go', 'high_std_no_go']:
        pts = []
        if mode == 'wgcna':
            f = open('./results/WGCNA_results/clus_info_WGCNA.txt', 'r')
        elif mode == 'no_go':
            f = open('./results/clus_info_no_go_40.txt', 'r')
        elif mode == 'go':
            f = open('./results/clus_info_go_41.txt', 'r')
        elif mode == 'ppi_go':
            f = open('./results/clus_info_go_35.txt', 'r')
        elif mode == 'ppi_no_go':
            f = open('./results/clus_info_no_go_36.txt', 'r')
        elif mode == 'high_std_no_go':
            f = open('./results/clus_info_no_go_42.txt', 'r')           
        for i, line in enumerate(f):
            if i < 3:
                continue
            line = line.split()
            # Skip the clusters with fewer than 50 genes.
            if float(line[4]) < 50:
                continue
            in_dens = float(line[1])
            top_enrichment_p = -math.log(float(line[8]), 10)
            pts += [(top_enrichment_p, in_dens)]
        f.close()

        if mode == 'wgcna':
            plt.scatter(*zip(*pts), color='blue', label=mode)
        elif mode == 'no_go':
            plt.scatter(*zip(*pts), color='gray', label=mode)
            med_in_dens = np.median([pt[1] for pt in pts]) * THRESHOLD_MULT
            plt.axhline(med_in_dens)
        elif mode == 'go':
            plt.scatter(*zip(*pts), color='red', label=mode)
        elif mode == 'ppi_go':
            plt.scatter(*zip(*pts), color='black', label=mode)
        elif mode == 'ppi_no_go':
            plt.scatter(*zip(*pts), color='orange', label=mode)
        elif mode == 'high_std_no_go':
            plt.scatter(*zip(*pts), color='green', label=mode)



    plt.title('Cluster in-densities vs. Best enrichment p-values')
    plt.xlabel('negative log of lowest GO enrichment p-value')
    plt.ylabel('in-density')
    plt.legend(loc='upper right')
    plt.show()