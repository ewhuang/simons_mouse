### Author: Edward Huang

import file_operations
import os
import subprocess
import sys

### This script calls the perl script to evaluate clusters for a given run.
### Outputs the file to the results folder.

def clean_module_results():
    '''
    Reformats the output of WGCNA by removing any GO terms from the clusters.
    '''
    results_folder = './results/%s_results' % data_type
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    f = open('./data/%s_module_membership.txt' % data_type, 'r')
    out = open('%s/clusters.txt' % results_folder, 'w')
    # Write dummy header line.
    out.write('dummy_header\n')
    for i, line in enumerate(f):
        # Skip header.
        if i == 0:
            continue
        # Ignore color and membership columns.
        node, module, color, membership = line.split()
        # Skip garbage module.
        if module == '0':
            continue        
        # Remove quotiation marks around the ENSMUSG ID.
        node = node.strip('"')
        # Don't add in GO terms.
        if 'GO:' in node or ('ENSMUSG' not in node and 'ENSG' not in node):
            continue
        # Write out the line.
        out.write('Species 0\tGene %s\tCluster %s\n' % (node, module))
    out.close()
    f.close()

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print 'Usage: %s data_type network_to_compare' % sys.argv[0]
        exit()
    global data_type
    data_type, comparing_network = sys.argv[1:]
    assert data_type == 'mouse' or data_type.isdigit()
    assert comparing_network.isdigit()

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    clean_module_results()

    command = ('perl ../evaluate_clustering.pl "./results/%s_results/'
        'clusters.txt" "../data/%s_data/networks_no_go/real_network_no'
        '_go_%s.txt" > "./results/%s_results/cluster_eval.txt"' % (
        data_type, data_type, comparing_network, data_type))
    subprocess.call(command, shell=True)