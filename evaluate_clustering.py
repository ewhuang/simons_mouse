### Author: Edward Huang

import file_operations
from multiprocessing import Pool
import os
import sys
import time

### This script calls the perl script to evaluate clusters for a given run.
### Outputs the file to a results folder.
### Run time: 4.5 minutes

def main():
    if len(sys.argv) != 4:
        print 'Usage:python %s data_type objective_function run_num' % (
            sys.argv[0])
        exit()
    data_type, objective_function, run_num = sys.argv[1:]
    assert objective_function in ('oclode', 'schaeffer', 'wlogv', 'wgcna',
        'prosnet')
    assert run_num.isdigit()

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]
    # Reconstruct data_type for TCGA. e.g., prosnet_0 to prosnet_carcinoma
    if 'prosnet' in data_type:
        data_type = data_type.split('_')[1]
        if data_type.isdigit():
            data_type = file_operations.get_tcga_disease_list()[int(data_type)]
        data_type = 'prosnet_' + data_type

    # First, get a clean cluster file without GO terms.
    file_operations.create_clean_go_file(data_type, objective_function, run_num)

    # Only for non-wgcna purposes.
    clean_dct = {'go':'clean_', 'no_go':''}

    # multiprocess for non-wgcna runs.
    if objective_function != 'wgcna':
        pool = Pool(processes=2)
    for network in ['go', 'no_go']:
        # Generate directory.
        results_folder = './results/%s_results/%s/cluster_eval_%s' % (
            data_type, objective_function, network)
        if not os.path.exists(results_folder):
            os.makedirs(results_folder)
        # Always compare against the network without GO.
        command = ('perl ./evaluate_clustering.pl "./results/%s_results/%s/'
                    'clusters_%s/clusters_%s_%s%s.txt" "./data/%s_data/'
                    'networks_no_go/real_network_no_go_%s.txt" > "%s/'
                    'cluster_eval_%s_%s.txt"' % (
                        data_type, objective_function, network, network,
                        clean_dct[network], run_num, data_type, run_num,
                        results_folder, network, run_num))
        print command
        if objective_function != 'wgcna':
            pool.apply_async(os.system, (command,))
        else:
            os.system(command)
            break # We don't want to run no_go for WGCNA.
    if objective_function != 'wgcna':
        pool.close()
        pool.join()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))