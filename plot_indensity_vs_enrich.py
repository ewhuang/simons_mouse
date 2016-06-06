### Author: Edward Huang

import math
import sys
import matplotlib
import numpy as np

### This script plots a scatterplot for both networks with and without GO
### the in-density of each cluster against the negative log of its best
### enrichment p-value.

matplotlib.use('Agg')
THRESHOLD_MULT = 0.8
import pylab

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage:python %s run_num' % sys.argv[0]
        exit()
    run_num = sys.argv[1]

    for mode in ['go', 'no_go']:
        pts = []
        f = open('./results/clus_info_%s/clus_info_%s_%s.txt' % (mode,
            mode, run_num), 'r')
        for i, line in enumerate(f):
            if i < 3:
                continue
            line = line.split()
            in_out_ratio = math.log(float(line[3]), math.e)
            top_enrichment_p = -math.log(float(line[8]), 10)
            pts += [(top_enrichment_p, in_out_ratio)]
        f.close()
        if mode == 'go':
            fig = matplotlib.pyplot.figure()
            matplotlib.pyplot.scatter(*zip(*pts), color='r', label=mode)
            ax = fig.add_subplot(111)
            # for i, xy in enumerate(pts):
            #     if xy[0] < 10:
            #         continue
            #     # Annotate clusters that have high GO enrichment p values.
            #     ax.annotate('%d' % (i+1), xy=xy)
        else:
            med_in_dens = np.median([pt[1] for pt in pts]) * THRESHOLD_MULT
            matplotlib.pyplot.axhline(med_in_dens)
            matplotlib.pyplot.scatter(*zip(*pts), color='k', label=mode)
    print med_in_dens
    matplotlib.pyplot.title('Simulated Annealing 20 Clusters, lambda=1.0')
    matplotlib.pyplot.xlabel('negative log of lowest GO enrichment p-value')
    matplotlib.pyplot.ylabel('log of weighted in-density/out-density ratio')
    matplotlib.pyplot.ylim(0, 10)
    matplotlib.pyplot.legend(loc='upper right')
    matplotlib.pyplot.show()
    pylab.savefig('./results/wgcna_cluster_plots_%s.png' % (run_num))