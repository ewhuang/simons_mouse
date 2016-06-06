### Author: Edward Huang

import file_operations
import numpy as np
from scipy.stats import pearsonr, betai
import time

### This file generates the base file that we use to create the network with
### GO file and the network without GO file. Reads in the raw data, and
### finds the Pearson correlation coefficient between every pair of genes,
### making an edge in the generated network if the coefficient exceeds a 
### set threshold. Each line is "GENE_1\tGENE_2\tPEARSON_SCORE".
### Run time: 7 minutes.

P_VALUE_THRESHOLD = 0.9

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
    # Read in the tsv file.
    gene_exp_dct = file_operations.get_gene_expression_dct()
    high_std_genes = file_operations.get_high_std_genes()

    gene_exp_matrix = create_gene_exp_matrix(gene_exp_dct, high_std_genes)

    r, p = corrcoef(gene_exp_matrix)

    out = open('./data/high_std_ensmusg_network.txt', 'w')
    for row_index, row in enumerate(r):
        for col_index, pcc in enumerate(row):
            if col_index <= row_index or pcc < P_VALUE_THRESHOLD or pcc == 1:
                continue
            # Write out gene information.
            gene_a, gene_b = (high_std_genes[row_index], high_std_genes[
                col_index])
            out.write('%s\t%s\t%f\n' % (gene_a, gene_b, abs(pcc)))
    out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))