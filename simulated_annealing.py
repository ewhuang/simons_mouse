### Author: Edward Huang

import file_operations
import subprocess
import sys
import time

### This script calls the C++ code to cluster on the network files.
### Command line arguments: temperature, number of clusters, and run number.
### Run time: 1 hour for gene edge weight threshold 0.8

def main():
    if len(sys.argv) not in [5, 6]:
        print 'Usage:python %s data_type objective_function run_num go/no_go go_num <if go>' % sys.argv[0]
        exit()

    # Sort out command line arguments.
    data_type, objective_function, run_num, network = sys.argv[1:5]
    assert data_type in ['mouse', 'tcga']
    assert objective_function in ['oclode', 'schaeffer', 'wlogv']
    assert run_num.isdigit()
    assert network in ['go', 'no_go']
    
    config_dct = file_operations.read_config_file(data_type)[run_num]
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
    command += './data/%s_networks_%s/network_%s_%s' % (data_type, network,
        network, run_num)

    if network == 'go':
        go_num = sys.argv[5]
        command += '_%s' % go_num

    command += '.txt -t %s 2>log > ./%s_results/%s/clusters_%s/clusters_%s_%s' % (
        temp, data_type, objective_function, network, network, run_num)

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