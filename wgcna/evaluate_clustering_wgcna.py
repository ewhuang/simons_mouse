### Author: Edward Huang

import subprocess
import sys

### This script calls the perl script to evaluate clusters for a given run.
### Outputs the file to the results folder.

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: %s genes_only/pca/mean/median' % sys.argv[0]
        exit()
    go_method = sys.argv[1]
    assert go_method in ['genes_only', 'pca', 'mean', 'median']

    if go_method == 'genes_only':
        go_domain_list = ['genes_only']
    else:
        go_domain_list = ['bp', 'cc', 'mf']

    for go_domain in go_domain_list:
        command = 'perl ../sim_anneal/evaluate_clustering.pl ./%s_results/' % (
            go_method)
        command += 'clusters_%s.txt ../data/networks_no_go/real_network_no' % (
            go_domain)
        command += '_go_42.txt > ./%s_results/cluster_eval_%s.txt' % (go_method,
            go_domain)

        subprocess.call(command, shell=True)