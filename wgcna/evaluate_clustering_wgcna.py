### Author: Edward Huang

import subprocess
import sys

### This script calls the perl script to evaluate clusters for a given run.
### Outputs the file to the results folder.

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: %s data_type genes_only/pca/mean/median' % sys.argv[0]
        exit()
    data_type = sys.argv[1]
    assert data_type in ['mouse', 'tcga']
    go_method = sys.argv[2]
    assert go_method in ['genes_only', 'pca', 'mean', 'median']

    if go_method == 'genes_only':
        go_domain_list = ['genes_only']
    else:
        go_domain_list = ['bp']

    for go_domain in go_domain_list:
        command = 'perl ../evaluate_clustering.pl ./%s_results/%s/' % (
            data_type, go_method)
        command += 'clusters_%s.txt ../data/%s_networks_no_go/real_network_no' % (
            go_domain, data_type)
        command += '_go_1.txt > ./%s_results/%s/cluster_eval_%s.txt' % (
            data_type, go_method, go_domain)

        subprocess.call(command, shell=True)