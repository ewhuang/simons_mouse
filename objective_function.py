### Author: Edward Huang

import sys
import numpy as np

### This script looks at the output files from analyze_clusters.py. First we
### look at clus_inf_no_go_RUNNUM.txt, and then find the median of the
### in-density for all of the clusters. We take a threshold of that median,
### and then look at the lusters with GO that meet that threshold. We then count
### the number of GO terms in these clusters, and print that information out.

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
        line = line.split('\t')
        in_dens += [float(line[0])]
    f.close()
    # The in-density threshold that clusters with GO terms must meet. 
    in_dens_threshold = THRESHOLD_MULT * np.median(in_dens)


    print 'Clusters with GO for run number %s...' % run_num
    print 'Median in-density for clusters without GO: %f' % in_dens_threshold

    # Now, check the clusters with GO terms.
    num_go_terms_by_cluster = []
    go_terms_by_cluster = []
    f = open('./results/clus_info_go_%s.txt' % run_num, 'r')
    for i, line in enumerate(f):
        if i < 3:
            continue
        line = line.strip().split('\t')
        # Only consider the clusters that meet the threshold
        if float(line[0]) < in_dens_threshold:
            continue
        num_go_terms_by_cluster += [line[4]]
        go_terms_by_cluster += [line[7:]]

    print 'Number of GO terms in clusters with in-density above threshold: '
    print ', '.join(num_go_terms_by_cluster)
    print (go_terms_by_cluster)
    f.close()