### Author: Edward Huang

import math
import sys
import matplotlib.pyplot as plt
import numpy as np

### This script plots a scatterplot for both networks with and without GO
### the in-density of each cluster against the negative log of its best
### enrichment p-value.

THRESHOLD_MULT = 0.8

if __name__ == '__main__':
    for mode in ['wgcna', 'no_go']:
        pts = []
        if mode == 'wgcna':
            f = open('./results/WGCNA_results/clus_info_WGCNA.txt', 'r')
        else:
            f = open('./results/clus_info_no_go_20.txt', 'r')            
        for i, line in enumerate(f):
            if i < 3:
                continue
            line = line.split()
            in_dens = float(line[1])
            top_enrichment_p = -math.log(float(line[8]), 10)
            pts += [(top_enrichment_p, in_dens)]
        f.close()


        if mode == 'wgcna':
            plt.scatter(*zip(*pts), color='r', label=mode)
        else:
            plt.scatter(*zip(*pts), color='k', label=mode)
            med_in_dens = np.median([pt[1] for pt in pts]) * THRESHOLD_MULT
            plt.axhline(med_in_dens)


    plt.title('Cluster in-densities vs. Best enrichment p-values, lambda=2.5')
    plt.xlabel('negative log of lowest GO enrichment p-value')
    plt.ylabel('in-density')
    plt.legend(loc='upper right')
    plt.show()