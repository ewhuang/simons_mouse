### Author: Edward Huang

import file_operations
import subprocess
import sys
import time

### This script calls the C++ code to cluster on the network files.
### Command line arguments: temperature, number of clusters, and run number.
### Run time: 1 hour for gene edge weight threshold 0.8

def main():
    if len(sys.argv) not in [4, 5]:
        print 'Usage:python %s oclode/schaeffer/wlogv run_num go/no_go go_num <if go>' % sys.argv[0]
        exit()

    # Sort out command line arguments.
    objective_function, run_num, network = sys.argv[1], sys.argv[2], sys.argv[3]
    assert objective_function in ['oclode', 'schaeffer', 'wlogv']
    assert (network in ['go', 'no_go'])
    config_dct = file_operations.read_config_file()[run_num]
    temp = config_dct['temp']
    num_clusters = config_dct['num_clusters']

    # Get the binary associated with the desired objective function.
    if objective_function == 'oclode':
        binary = '''./OCLODE/makedir/OCLODE_efficient'''
    elif objective_function == 'schaeffer':
        binary = './SchaefferScore/makedir/SchaefferImplementNotWeighted'
    else:
        binary = './WlogV/makedir/WlogVImplement'

    # Build the command.
    command = '%s %s 1 0 ./data/orth.txt 1 ' % (binary, num_clusters)
    command += './data/networks_%s/network_%s_%s' % (network, network, run_num)

    if network == 'go':
        go_num = sys.argv[4]
        command += '_%s' % go_num

    command += '.txt -t %s 2>log > ./results/%s/clusters_%s/clusters_%s_%s' % (
        temp, objective_function, network, network, run_num)

    if network == 'go':
        command += '_%s' % go_num

    command += '.txt'

    # Execute the command in the shell.
    print command
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))