### Author: Edward Huang

import file_operations
import subprocess
import sys
import time

### This script calls the perl script to evaluate clusters for a given run.
### Outputs the file to the results folder. File names are
### ./results/cluster_eval_no_go_RUNNUM.txt
### Run time: 4.5 minutes

def main():
    if len(sys.argv) != 3:
        print 'Usage:python %s objective_function run_num' % sys.argv[0]
        exit()
    objective_function = sys.argv[1]
    assert objective_function in ['oclode', 'schaeffer', 'wlogv']
    run_num = sys.argv[2]

    for go_domain_num in range(2):
        command = 'perl ./evaluate_clustering.pl '
        # First, get a clean cluster file.
        file_operations.create_clean_go_file(objective_function, run_num,
            go_domain_num)
        command += './results/%s/clusters_go/clusters_go_clean_%s_%d.txt ' % (
            objective_function, run_num, go_domain_num)
        command += './data/networks_go/real_network_go_%s_%d.txt ' % (run_num,
            go_domain_num)
        command += '> ./results/%s/cluster_eval_go/cluster_eval_go_%s_%d.txt' % (
            objective_function, run_num, go_domain_num)
        subprocess.call(command, shell=True)

    # No GO term command.
    command = 'perl ./evaluate_clustering.pl '
    command += './results/%s/clusters_no_go/clusters_no_go_%s.txt ' % (
        objective_function, run_num)
    command += './data/networks_no_go/real_network_no_go_%s.txt ' % (run_num)
    command += '> ./results/%s/cluster_eval_no_go/cluster_eval_no_go_%s.txt' % (
        objective_function, run_num)
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))