# find the median of in-density for clusterings without GO terms. Take something
# like 0.8 of that number. Look at the clusters WITH GO terms. Clusters that meet
# our threshold will be considered. We want to have a high number of GO terms
# for those clusters that meet that threshold. This is our "objective function".

# To do this, change the lambda values for edges between gene-go edges, and redo
# clustering.

### Author: Edward Huang

import sys
import numpy as np

### This script looks at the output files from analyze_clusters.py. First we
### look at clus_inf_no_go_RUNNUM.txt, and then find the median of the
### in-density for all of the clusters. 

THRESHOLD_MULT = 0.8

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage:python %s run_num' % sys.argv[0]
        exit()
    run_num = sys.argv[1]

    # The list of in-densities for the clusters without GO terms.
    in_dens = []
    # First open cluster information for networks without GO terms.
    f = open('./results/clus_info_no_go_%s.txt' % run_num, 'r')
    for i, line in enumerate(f):
        if i < 3:
            continue
        line = line.split()
        in_dens += [float(line[0])]
    f.close()
    # The in-density threshold that clusters with GO terms must meet. 
    in_dens_threshold = THRESHOLD_MULT * np.median(in_dens)

    # Now, check the clusters with GO terms.
    f = open('./results/clus_info_go_%s.txt' % run_num, 'r')
    for i, line in enumerate(f):
        if i < 3:
            continue
        line = line.split()
        # Only consider the clusters that meet the threshold
        if float(line[0]) < in_dens_threshold:
            continue
        print line[4]
    f.close()