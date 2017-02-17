### Author: Edward Huang

import file_operations
import numpy as np
import operator
from scipy.stats import betai
import sys
import time

# Reads in the raw data, and finds the Pearson correlation coefficient between
### every pair of genes, making an edge in the generated network if the
### coefficient exceeds a set threshold.
### Each line is "GENE_1\tGENE_2\tABS_PEARSON_SCORE".
### Run time: 50s for mouse, 100s for TCGA, 12min for all TCGA.

# def corrcoef(matrix):
#     '''
#     Computes Pearson correlation between every pair of genes.
#     Received code from following link:
#     http://stackoverflow.com/questions/24432101/correlation-coefficients-and-p
#         -values-for-all-pairs-of-rows-of-a-matrix
#     '''
#     r = np.corrcoef(matrix)
#     print 'ye'
#     rf = r[np.triu_indices(r.shape[0], 1)]
#     df = matrix.shape[1] - 2
#     ts = rf * rf * (df / (1 - rf * rf))
#     pf = betai(0.5 * df, 0.5, df / (df + ts))
#     p = np.zeros(shape=r.shape)
#     p[np.triu_indices(p.shape[0], 1)] = pf
#     p[np.tril_indices(p.shape[0], -1)] = pf
#     p[np.diag_indices(p.shape[0])] = np.ones(p.shape[0])
#     return r, p

def create_gene_exp_matrix(high_std_genes, gene_exp_dct):
    '''
    Make a 2D array, where rows are genes, columns are gene expression values.
    '''
    gene_exp_matrix = []
    for gene in high_std_genes:
        gene_exp_matrix += [gene_exp_dct[gene]]
    return np.array(gene_exp_matrix)

def write_edge_dct(high_std_genes, gene_exp_matrix, folder_name):
    '''
    Given the 2D gene expression matrix, compute the pairwise Pearson
    coefficients, and then make the edge dictionary. Returns the edges sorted
    by decreasing weight. Take only the top one million edges.
    '''
    num_genes = len(high_std_genes)
    num_edges = int(0.01 * num_genes * (num_genes - 1) / 2.0)

    edge_list = []

    # TODO: currently not using p-values.
    r = np.corrcoef(gene_exp_matrix)
    # Get the upper triangular matrix.
    r = np.triu(r, 1)
    r[r == 1] = 0

    x_idx, y_idx = np.unravel_index(np.argsort(r.ravel())[-num_edges:], r.shape)

    for x, y in zip(x_idx, y_idx):
        gene_x, gene_y = high_std_genes[x], high_std_genes[y]
        edge_list += ['%s\t%s\t%g' % (gene_x, gene_y, r[x][y])]

    out = open('./data/%s_data/high_std_network.txt' % folder_name, 'w')
    out.write('\n'.join(edge_list[::-1]))
    out.close()

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/all/tcga' % sys.argv[0]
        exit()
    gene_type = sys.argv[1]
    assert gene_type in ['mouse', 'all', 'tcga']

    if gene_type == 'all':
        folder_list = file_operations.get_tcga_disease_list()
    else:
        folder_list = [gene_type]

    for folder_name in folder_list:
        gene_exp_dct = file_operations.get_gene_expression_dct(folder_name)
        high_std_genes = file_operations.get_high_std_genes(folder_name)

        gene_exp_matrix = create_gene_exp_matrix(high_std_genes, gene_exp_dct)
        write_edge_dct(high_std_genes, gene_exp_matrix, folder_name)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))