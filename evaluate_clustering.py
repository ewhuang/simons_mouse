### Author: Edward Huang

import file_operations
import subprocess
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
    assert objective_function in ['oclode', 'schaeffer', 'wlogv']
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

    clean_dct = {'go':'clean_', 'no_go':''}

    for network in ['go', 'no_go']:
        command = ('perl ./evaluate_clustering.pl "./results/%s_results/%s/'
                    'clusters_%s/clusters_%s_%s%s.txt" "./data/%s_data/'
                    'networks_%s/real_network_%s_%s.txt" > "./results/'
                    '%s_results/%s/cluster_eval_%s/cluster_eval_%s_%s.txt"' % (
                        data_type, objective_function, network, network,
                        clean_dct[network], run_num, data_type, 'no_go',
                        'no_go', run_num, data_type, objective_function,
                        network, network, run_num))
        subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))