### Author: Edward Huang

import subprocess
import sys

### This script calls the perl script to evaluate clusters for a given run.
### Outputs the file to the results folder. File names are
### ./results/cluster_eval_no_go_RUNNUM.txt

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage:python %s run_num' % sys.argv[0]
        exit()
    run_num = sys.argv[1]

    for network in ['go', 'no_go']:

        command = 'perl ./sim_anneal/evaluate_clustering.pl '
        if network == 'go':
            command += './results/clusters_%s_clean_%s.txt ' % (network, run_num)
        else:
            command += './results/clusters_%s_%s.txt ' % (network, run_num)
        command += './data/real_network_%s_%s.txt ' % (network, run_num)
        command += '> ./results/cluster_eval_%s_%s.txt' % (network, run_num)
        subprocess.call(command, shell=True)


