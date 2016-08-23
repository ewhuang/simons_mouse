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
        print 'Usage:python %s data_type objective_function run_num' % sys.argv[0]
        exit()
    data_type = sys.argv[1]
    assert data_type == 'mouse' or data_type.isdigit()
    objective_function = sys.argv[2]
    assert objective_function in ['oclode', 'schaeffer', 'wlogv']
    run_num = sys.argv[3]
    assert run_num.isdigit()

    if data_type.isdigit():
        data_type = file_operations.get_tcga_diseases()[int(data_type)]

    for go_domain_num in [0]:
        command = 'perl ./evaluate_clustering.pl '
        # First, get a clean cluster file.
        file_operations.create_clean_go_file(data_type, objective_function,
            run_num, go_domain_num)
        command += '"./results/%s_results/%s/clusters_go/clusters_go_clean_%s_%d.txt" ' % (
            data_type, objective_function, run_num, go_domain_num)
        command += '"./data/%s_data/networks_go/real_network_go_%s_%d.txt" ' % (
            data_type, run_num, go_domain_num)
        command += '> "./results/%s_results/%s/cluster_eval_go/cluster_eval_go_%s_%d.txt"' % (
            data_type, objective_function, run_num, go_domain_num)
        subprocess.call(command, shell=True)

    # No GO term command.
    command = 'perl ./evaluate_clustering.pl '
    command += '"./results/%s_results/%s/clusters_no_go/clusters_no_go_%s.txt" ' % (
        data_type, objective_function, run_num)
    command += '"./data/%s_data/networks_no_go/real_network_no_go_%s.txt" ' % (
        data_type, run_num)
    command += '> "./results/%s_results/%s/cluster_eval_no_go/cluster_eval_no_go_%s.txt"' % (
        data_type, objective_function, run_num)
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))