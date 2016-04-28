### Author: Edward Huang

import file_operations
import subprocess
import sys

### This script calls the C++ code to cluster on the network files.
### Command line arguments: temperature, number of clusters, and run number.

if __name__ == '__main__':
    # if len(sys.argv) != 5:
    #     print 'Usage:python %s go/no_go temp num_clusters run_num' % sys.argv[0]
    #     exit()

    # network, temp, num_clusters, run_num = sys.argv[1:]
    # # Old non-trashcan clustering.
    # # command = './sim_anneal/bin/cs-grn %s 1 0 ' % num_clusters
    # # New trashcan clustering.
    # grn_folder = '/co-expression_network_clustering/InOutRatioModel_withNoiseCl'
    # grn_folder += 'us_noiseOrthoEdge_fast/src/bin/'
    # command = '.%sclustering %s 1 0 ' % (grn_folder, num_clusters)

    # command += './data/orth.txt 1 '
    # assert (network in ['go', 'no_go'])
    # command += './data/network_%s_%s.txt -n ' % (network, run_num)
    # command += '-t 1 2> log > '
    # command += './results/clusters_%s_%s.txt' % (network, run_num)

    # print command
    # subprocess.call(command, shell=True)
    # subprocess.call('rm log', shell=True)

    if len(sys.argv) != 3:
        print 'Usage:python %s run_num go/no_go' % sys.argv[0]
        exit()

    run_num, network = sys.argv[1:]
    config_dct = file_operations.read_config_file()[run_num]
    temp = config_dct['temp']
    num_clusters = config_dct['num_clusters']

    command = './sim_anneal/bin/cs-grn %s 1 0 ' % num_clusters
    command += './data/orth.txt 1 '
    assert (network in ['go', 'no_go'])
    if network == 'no_go':
        command += './data/network_%s_%s.txt ' % (network, run_num)
        command += '-t %s 2> log > ' % temp
        command += './results/clusters_%s_%s.txt' % (network, run_num)

        print command
        subprocess.call(command, shell=True)
    else:
        for i in range(3):
            command += './data/network_%s_%s_%d.txt ' % (network, run_num, i)
            command += '-t %s 2> log > ' % temp
            command += './results/clusters_%s_%s_%d.txt' % (network, run_num, i)

            print command
            subprocess.call(command, shell=True)