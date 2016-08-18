### Author: Edward Huang

import file_operations
import os
import sys

### Reformats the output of WGCNA for running through simulated annealing
### evluation scripts. Takes out any GO terms.

def main():
    if len(sys.argv) != 3:
        print 'Usage: %s data_type genes_only/pca/mean/median' % sys.argv[0]
        exit()
    data_type = sys.argv[1]
    assert data_type == 'mouse' or data_type.isdigit()
    go_method = sys.argv[2]
    assert go_method in ['genes_only', 'pca', 'mean', 'median']

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    if go_method == 'genes_only':
        domain_list = [go_method]
    else:
        domain_list = ['bp', 'cc', 'mf']

    results_folder = './results/%s_results/' % data_type
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    results_folder = './results/%s_results/%s/' % (data_type, go_method)
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)

    for go_domain in domain_list:
        f = open('./data/%s_module_membership_%s.txt' % (data_type, go_domain),
            'r')
        out = open('%sclusters_%s.txt' % (results_folder, go_domain), 'w')
        # Write dummy header line.
        out.write('dummy_header\n')

        for i, line in enumerate(f):
            # Skip header.
            if i == 0:
                continue

            # Ignore color and membershi columns.
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
    main()