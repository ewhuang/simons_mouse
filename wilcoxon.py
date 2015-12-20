### Author: Edward Huang

import sys
from scipy.stats import ranksums

### This script takes the results of analyze_clusters.py and then analyzes the
### in-density and out-density columns for cluster with and without GO to
### determine the change.

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage:python %s run_num' % sys.argv[0]
        exit()
    run_num = sys.argv[1]

    go_f = open('./results/clus_info_go_%s.txt' % run_num, 'r')
    in_go = []
    out_go = []
    ratio_go = []
    for i, line in enumerate(go_f):
        if i < 3:
            continue
        line = line.split()
        in_dens, out_dens, ratio = line[0], line[1], line[2]
        in_go += [float(in_dens)]
        out_go += [float(out_dens)]
        ratio_go += [float(ratio)]

    go_f.close()

    no_go_f = open('./results/clus_info_no_go_%s.txt' % run_num, 'r')
    in_no_go = []
    out_no_go = []
    ratio_no_go = []
    for i, line in enumerate(no_go_f):
        if i < 3:
            continue
        line = line.split()
        in_dens, out_dens, ratio = line[0], line[1], line[2]
        in_no_go += [float(in_dens)]
        out_no_go += [float(out_dens)]
        ratio_no_go += [float(ratio)]
    no_go_f.close()

    print 'In-Density z-stat, p-value: ', ranksums(in_go, in_no_go)
    print 'Out-Density z-stat, p-value: ', ranksums(out_go, out_no_go)
    print 'In/(In + Out) z-stat, p-value: ', ranksums(ratio_go, ratio_no_go)