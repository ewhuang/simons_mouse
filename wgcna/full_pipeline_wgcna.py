### Author: Edward Huang

import subprocess
import sys

### This script runs the full pipeline after running WGCNA and placing the
### resulting file in ./wgcna/data/.
### data_type is either 'mouse' or a number, corresponding to a TCGA disease.

def main():
    if len(sys.argv) != 3:
        print 'Usage:python %s data_type network_num' % sys.argv[0]
        exit()
    data_type, network_num = sys.argv[1:]
    assert data_type == 'mouse' or data_type.isdigit()
    assert network_num.isdigit()

    command = 'python clean_wgcna_module_results.py %s' % data_type
    subprocess.call(command, shell=True)

    command = 'python evaluate_clustering_wgcna.py %s %s' % (data_type,
        network_num)
    subprocess.call(command, shell=True)

    command = 'python compute_label_enrichments_wgcna.py %s go' % data_type
    subprocess.call(command, shell=True)

    command = 'python compute_label_enrichments_wgcna.py %s dbgap' % data_type
    subprocess.call(command, shell=True)

    command = 'python cluster_info_summary_wgcna.py %s' % data_type
    subprocess.call(command, shell=True)

if __name__ == '__main__':
    main()