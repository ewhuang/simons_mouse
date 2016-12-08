### Author: Edward Huang

import file_operations
import os
import subprocess
import sys
import time

### This script calls the C++ code to cluster on the network files.
### Command line arguments: temperature, number of clusters, and run number.
### Run time: 1 hour for gene edge weight threshold 0.8

def make_folders():
    '''
    Construct the results folders if it's the first time running simulated
    annealing.
    '''
    results_folder = './results/%s_results/' % data_type
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    obj_func_folder = results_folder + 'wlogv/'
    if not os.path.exists(obj_func_folder):
        os.makedirs(obj_func_folder)
        for network_type in ['go', 'no_go']:
            os.makedirs(obj_func_folder + 'clus_info_%s/' % network_type)
            os.makedirs(obj_func_folder + 'cluster_enrichment_terms_%s/' % 
                network_type)
            os.makedirs(obj_func_folder + 'dbgap_enrichment_terms_%s/' %
                network_type)
            os.makedirs(obj_func_folder + 'clusters_%s/' % network_type)
            os.makedirs(obj_func_folder + 'cluster_eval_%s/' % network_type)
    plot_folder = results_folder + 'comparison_plots/'
    if not os.path.exists(plot_folder):
        os.makedirs(plot_folder)

def main():
    if len(sys.argv) != 5:
        print ('Usage:python %s data_type objective_function run_num '
                'go/no_go' % sys.argv[0])
        exit()
    global data_type
    data_type, objective_function, run_num, network = sys.argv[1:5]
    assert objective_function in ['oclode', 'schaeffer', 'wlogv']
    assert run_num.isdigit()
    assert network in ['go', 'no_go']

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]
    # Reconstruct data_type for TCGA. e.g., prosnet_0 to prosnet_carcinoma
    if 'prosnet' in data_type:
        data_type = data_type.split('_')[1]
        if data_type.isdigit():
            data_type = file_operations.get_tcga_disease_list()[int(data_type)]
        data_type = 'prosnet_' + data_type

    config_dct = file_operations.read_config_file(data_type)[run_num]
    temp, num_clusters = config_dct['temp'], config_dct['num_clusters']

    make_folders()

    # Get the binary associated with the desired objective function.
    if objective_function == 'oclode':
        binary = '''./OCLODE/makedir/OCLODE_efficient'''
    elif objective_function == 'schaeffer':
        binary = './SchaefferScore/makedir/SchaefferImplementNotWeighted'
    else:
        binary = './WlogV/makedir/WlogVImplement'

    # Build the command.
    command = ('%s %s 1 0 "./data/%s_data/orth.txt" 1 "./data/%s_data/'
                'networks_%s/network_%s_%s.txt" -t %s 2>log > "./results/'
                '%s_results/%s/clusters_%s/clusters_%s_%s.txt"' % (binary,
                    num_clusters, data_type, data_type, network, network,
                    run_num, temp, data_type, objective_function, network,
                    network, run_num))
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))