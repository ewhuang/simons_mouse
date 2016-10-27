### Author: Edward Huang

import sys
import subprocess

### This script runs the full pipeline after running WGCNA and placing the
### resulting file in ./wgcna/data/.
### data_type is either 'mouse' or a number, corresponding to a TCGA disease.

def main():
    if len(sys.argv) != 4:
        print 'Usage:python %s data_type genes_only/pca/mean/median network_num' % sys.argv[0]
        exit()
    data_type = sys.argv[1]
    assert data_type == 'mouse' or data_type.isdigit()
    go_type = sys.argv[2]
    assert go_type in ['genes_only', 'pca', 'mean', 'median']
    network_num = sys.argv[3]
    assert network_num.isdigit()

    command = 'python clean_wgcna_module_results.py %s %s' % (data_type,
        go_type)
    subprocess.call(command, shell=True)

    command = 'python evaluate_clustering_wgcna.py %s %s %s' % (data_type,
        go_type, network_num)
    subprocess.call(command, shell=True)

    command = 'python compute_go_enrichment_wgcna.py %s %s' % (data_type,
        go_type)
    subprocess.call(command, shell=True)

    command = 'python compute_dbgap_enrichment_wgcna.py %s %s' % (data_type,
        go_type)
    subprocess.call(command, shell=True)

    command = 'python cluster_info_summary_wgcna.py %s %s %s' % (data_type,
        go_type, network_num)
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    main()