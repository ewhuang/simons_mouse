### Author: Edward Huang

import sys

### This file removes the GO nodes from the network_go.txt files
### so that we can evaluate them without the GO nodes inside the
### the clusters.

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage:python %s run_num' % sys.argv[0]
        exit()
    run_num = sys.argv[1]

    f = open('./results/clusters_go_%s.txt' % run_num, 'r')
    out = open('./results/clusters_go_clean_%s.txt' % run_num, 'w')
    for i, line in enumerate(f):
        if i == 0:
            out.write(line)
            continue
        newline = line.split()
        if 'ENSMUSG' not in newline[3]:
            continue
        out.write(line)
    out.close()
    f.close()