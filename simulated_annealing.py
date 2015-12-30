### Author: Edward Huang

import subprocess
import sys

### This script calls the C++ code to cluster on the network files.
### Command line arguments: temperature, number of clusters, and run number.

if __name__ == '__main__':
    if len(sys.argv) != 5:
        print 'Usage:python %s go/no_go temp num_clusters run_num' % sys.argv[0]
        exit()

    network, temp, num_clusters, run_num = sys.argv[1:]
    command = './sim_anneal/bin/cs-grn %s 1 0 ' % num_clusters
    command += './data/orth.txt 1 '
    assert (network in ['go', 'no_go'])
    command += './data/network_%s_%s.txt ' % (network, run_num)
    command += '-t 1 2> log > '
    command += './results/clusters_%s_%s.txt' % (network, run_num)

    print command
    subprocess.call(command, shell=True)