### Author: Edward Huang

import json
import numpy as np
from sklearn.decomposition import PCA
from scipy.stats import pearsonr
import sys
import time

### This script preprocesses the original coexpression file to prepare for
### WGCNA. It looks at the networks we have previously built, finds all the
### genes used in those networks, and then preserves only the rows in the 
### original file that matches those genes.
### Run time: 9 seconds.

def get_high_std_genes(data_type):
    '''
    Retrieves the list of genes that have high standard deviations across their
    gene expression vectors.
    '''
    high_std_genes = []
    f = open('../data/%s_data/high_std_genes.txt' % data_type, 'r')
    for line in f:
        gene = line.strip()
        assert 'ENSMUSG' in gene or 'ENSG' in gene
        high_std_genes += [gene]
    f.close()
    return high_std_genes

def main():
    if len(sys.argv) != 3:
        print 'Usage: %s data_type genes_only/pca/mean/median' % sys.argv[0]
        exit()
    data_type = sys.argv[1]
    # assert data_type in ['mouse', 'tcga']
    go_method = sys.argv[2]
    assert go_method in ['genes_only', 'pca', 'mean', 'median']

    high_std_genes = get_high_std_genes(data_type)

    if go_method == 'genes_only':
        go_domain_list = [go_method]
    else:
        go_domain_list = ['bp', 'cc', 'mf']

    # go_domain is the domain we hold out.
    for go_domain in go_domain_list:
        ## Leave one out code. This block adds two GO domains into the network,
        ## and evaluates on the domain that is not added.
        # go_domain_list_train = go_domain_list[:]
        # go_domain_list_train.remove(go_domain)

        # go_gene_dct = {}
        # for go_domain_train in go_domain_list_train:
        #     go_gene_dct.update(get_go_gene_dct(go_domain_train))

        # # This line is to both train and evaluate on the same domain.
        # if go_domain != 'genes_only':
        #     go_gene_dct = get_go_gene_dct(go_domain)

        if data_type == 'mouse':
            f = open('../data/mouse_data/mm_mrsb_log2_expression.tsv', 'r')
            out = open('./data/mm_mrsb_log2_expression_%s.tsv' % go_domain, 'w')
        else:
            f = open('../data/%s_data/expr.txt' % data_type, 'r')
            out = open('./data/%s_expr_%s.txt' % (data_type, go_domain), 'w')

        # Write out the gene expression vectors for genes.
        # gene_expression_dct = {}
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
            # exp_vals = [float(val) for val in exp_vals]
            # assert gene not in gene_expression_dct
            # gene_expression_dct[gene] = exp_vals
        f.close()

        # if go_domain == 'genes_only':
        #     break

        # # Run PCA and generate gene expression vectors for GO terms.
        # for go_term in go_gene_dct:
        #     annotated_gene_list = go_gene_dct[go_term]

        #     # Skip bad GO terms.
        #     num_annotated_gene_list = len(annotated_gene_list)
        #     if num_annotated_gene_list < 10 or num_annotated_gene_list > 1000:
        #         continue

        #     # Construct expression matrix for a GO term on the genes it
        #     # annotates.
        #     super_gene_matrix = []
        #     for annotated_gene in annotated_gene_list:
        #         ann_gene_expression = gene_expression_dct[annotated_gene]
        #         super_gene_matrix += [ann_gene_expression]
            
        #     if go_method == 'mean':
        #         # Get the average across each sample.
        #         super_gene_matrix = np.array(super_gene_matrix)
        #         super_gene = np.mean(super_gene_matrix, axis=0)
        #     elif go_method == 'pca':
        #         # Get most principal component.
        #         pca = PCA()
        #         pca.fit(super_gene_matrix)
        #         super_gene = pca.components_[0] + pca.mean_
        #     elif go_method == 'median':
        #         # Get median across each sample.
        #         super_gene_matrix = np.array(super_gene_matrix)
        #         super_gene = np.median(super_gene_matrix, axis=0)

        #     # Write out the vector to file.
        #     out.write(go_term + '\t')
        #     out.write('\t'.join(map(str, super_gene)) + '\n')

        out.close()
        f.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))