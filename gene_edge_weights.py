### Author: Edward Huang

import file_operations
import numpy as np
from scipy.stats import pearsonr, betai
import sys
import time

### This file generates the base file that we use to create the network with
### GO file and the network without GO file. Reads in the raw data, and
### finds the Pearson correlation coefficient between every pair of genes,
### making an edge in the generated network if the coefficient exceeds a 
### set threshold. Each line is "GENE_1\tGENE_2\tPEARSON_SCORE".
### Run time: 7 minutes.

# P_VALUE_THRESOLD = 0.0001

def corrcoef(matrix):
    '''
    Received code from following link:
    http://stackoverflow.com/questions/24432101/correlation-coefficients-and-p-v
    alues-for-all-pairs-of-rows-of-a-matrix
    '''
    r = np.corrcoef(matrix)
    rf = r[np.triu_indices(r.shape[0], 1)]
    df = matrix.shape[1] - 2
    ts = rf * rf * (df / (1 - rf * rf))
    pf = betai(0.5 * df, 0.5, df / (df + ts))
    p = np.zeros(shape=r.shape)
    p[np.triu_indices(p.shape[0], 1)] = pf
    p[np.tril_indices(p.shape[0], -1)] = pf
    p[np.diag_indices(p.shape[0])] = np.ones(p.shape[0])
    return r, p

def create_gene_exp_matrix(gene_exp_dct, high_std_genes):
    gene_exp_matrix = []
    for gene in high_std_genes:
        gene_exp_matrix += [gene_exp_dct[gene]]
    return np.array(gene_exp_matrix)

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga_cancer_index' % sys.argv[0]
        exit()
    data_type = sys.argv[1]    
    assert data_type == 'mouse' or data_type.isdigit()
    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    if data_type == 'mouse':
        pcc_threshold = 0.9
    else:
        # TCGA coefficients are worse.
        pcc_threshold = 0.5

    # Read in the tsv file.
    gene_exp_dct = file_operations.get_gene_expression_dct(data_type)
    high_std_genes = file_operations.get_high_std_genes(data_type)

    gene_exp_matrix = create_gene_exp_matrix(gene_exp_dct, high_std_genes)

    r, p = corrcoef(gene_exp_matrix)

    out = open('./data/%s_data/high_std_network.txt' % data_type, 'w')
    for row_idx, row in enumerate(r):
        for col_idx, pcc in enumerate(row):
            if col_idx <= row_idx or pcc < pcc_threshold or pcc == 1:
                continue
            # if p[row_idx][col_idx] > P_VALUE_THRESOLD:
            #     continue
            # Write out gene information.
            gene_a, gene_b = (high_std_genes[row_idx], high_std_genes[col_idx])
            out.write('%s\t%s\t%f\n' % (gene_a, gene_b, abs(pcc)))
    out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))