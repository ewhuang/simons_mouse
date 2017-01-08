### Author: Edward Huang

import file_operations
import sys
import time

### This script pre-processes the original co-expression file to prepare for
### WGCNA. It preserves only the rows that correspond to the top 15k genes with
### the highest standard deviations across all samples' gene expression values.
### Run time: 10 for mouse, 80 seconds for tcga.

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga' % sys.argv[0]
        exit()
    category = sys.argv[1]
    assert category in ['mouse', 'tcga']

    if category == 'mouse':
        data_type_list = ['mouse']
    elif category == 'tcga':
        # Pre-processes for all TCGA cancers.
        data_type_list = file_operations.get_tcga_disease_list()

    # data_type can be 'mouse' or any of the TCGA cancers.
    for data_type in data_type_list:
        high_std_genes = file_operations.get_high_std_genes(data_type)
        f = open('../data/%s_data/expr.tsv' % data_type, 'r')
        out = open('./data/%s_expr.tsv' % data_type, 'w')
        for i, line in enumerate(f):
            # Directly write out the header file.
            if i == 0:
                out.write(line)
                continue
            split_line = line.split()
            gene, exp_vals = split_line[0], split_line[1:]
            assert 'ENSMUSG' in gene or 'ENSG' in gene            
            # Skip genes with low standard deviation.
            if gene not in high_std_genes:
                continue
            out.write(line)
        f.close()
        out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))