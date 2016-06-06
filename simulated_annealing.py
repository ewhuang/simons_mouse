### Author: Edward Huang

import file_operations
import subprocess
import sys
import time

### This script calls the C++ code to cluster on the network files.
### Command line arguments: temperature, number of clusters, and run number.
### Run time: 40 minutes.

binary = '''./OCLODE/makedir/OCLODE_efficient'''
# binary = './sim_anneal/bi/cs-grn'

def main():
    if len(sys.argv) != 3 and len(sys.argv) != 4:
        print 'Usage:python %s run_num go/no_go optional:go_num' % sys.argv[0]
        exit()

    run_num, network = sys.argv[1], sys.argv[2]
    config_dct = file_operations.read_config_file()[run_num]
    temp = config_dct['temp']
    num_clusters = config_dct['num_clusters']

    assert (network in ['go', 'no_go'])
    if network == 'no_go':
        assert len(sys.argv) == 3
        command = '%s %s 1 0 ' % (binary, num_clusters)
        command += './data/orth.txt 1 '
        command += './data/networks_%s/network_%s_%s.txt ' % (network, network,
            run_num)
        command += '-t %s 2> log > ' % temp
        command += './results/clusters_%s/clusters_%s_%s.txt' % (network,
            network, run_num)

        print command
        subprocess.call(command, shell=True)
    else:
        go_num = sys.argv[3]
        assert go_num in ['0', '1', '2']
        command = '%s %s 1 0 ' % (binary, num_clusters)
        command += './data/orth.txt 1 '
        command += './data/networks_%s/network_%s_%s_%s.txt ' % (network, 
            network, run_num, go_num)
        command += '-t %s 2> log > ' % temp
        command += './results/clusters_%s/clusters_%s_%s_%s.txt' % (network,
            network, run_num, go_num)

        print command
        subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))