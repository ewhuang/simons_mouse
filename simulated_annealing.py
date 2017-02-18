### Author: Edward Huang

import file_operations
import os
import subprocess
import sys
import time

### This script calls the C++ code to cluster on the network files.
### Command line arguments: temperature, number of clusters, and run number.
### Run time: 1 hour for gene edge weight threshold 0.8

def make_folders(data_type, obj_func):
    '''
    Construct the results folders if it's the first time running simulated
    annealing.
    '''
    results_folder = './results/%s_results' % data_type
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    obj_func_folder = '%s/%s' % (results_folder, obj_func)
    if not os.path.exists(obj_func_folder):
        os.makedirs(obj_func_folder)

def get_num_clusters(data_type):
    num_clusters = -1
    wgcna_f = open('./results/%s_results/wgcna/clus_info_go/clus_info_go_1'
        '.tsv' % data_type, 'r')
    for line in wgcna_f:
        num_clusters += 1
    wgcna_f.close()
    return num_clusters

def main():
    if len(sys.argv) not in [4, 5]:
        print ('Usage:python %s mouse/tcga/tcga_idx obj_func net_type <prosnet>'
            % sys.argv[0])
        exit()
    data_type, obj_func, net_type = sys.argv[1:4]
    assert obj_func in ['oclode', 'schaeffer', 'wlogv']
    assert net_type in ['none', 'go']

    if data_type.isdigit():
        data_type = file_operations.get_tcga_list()[int(data_type)]

    num_clusters = get_num_clusters(data_type)

    if len(sys.argv) == 5:
        assert sys.argv[4] == 'prosnet'
        data_type = 'prosnet_' + data_type

    # config_dct = file_operations.read_config_file(data_type)[run_num]
    # temp, num_clusters = config_dct['temp'], config_dct['num_clusters']

    temp = 10
    make_folders(data_type, obj_func)

    # Get the binary associated with the desired objective function.
    if obj_func == 'oclode':
        binary = '''./OCLODE/makedir/OCLODE_efficient'''
    elif obj_func == 'schaeffer':
        binary = './SchaefferScore/makedir/SchaefferImplementNotWeighted'
    else:
        binary = './WlogV/makedir/WlogVImplement'

    # Build the command.
    command = ('%s %s 1 0 "./data/%s_data/orth.txt" 1 "./data/%s_data/'
                'network_%s.txt" -t %s 2>log > "./results/%s_results/%s/'
                'clusters_%s.txt"' % (binary, num_clusters, data_type,
                    data_type, net_type, temp, data_type, obj_func, net_type))
    print command
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))