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
### Run time: 7 minutes. An hour and 40 minutes for TCGA datasets.

def corrcoef(matrix):
    '''
    Computes Pearson correlation between every pair of genes.
    Received code from following link:
    http://stackoverflow.com/questions/24432101/correlation-coefficients-and-p
        -values-for-all-pairs-of-rows-of-a-matrix
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

def create_gene_exp_matrix(high_std_genes, gene_exp_dct):
    '''
    Make a 2D array, where rows are genes, columns are gene expression values.
    '''
    gene_exp_matrix = []
    for gene in high_std_genes:
        gene_exp_matrix += [gene_exp_dct[gene]]
    return np.array(gene_exp_matrix)

def create_sorted_edge_dct(high_std_genes, gene_exp_matrix):
    '''
    Given the 2D gene expression matrix, compute the pairwise Pearson
    coefficients, and then make the edge dictionary. Returns the edges sorted
    by decreasing weight. Take only the top one million edges.
    '''
    edge_dct = {}

    r, p = corrcoef(gene_exp_matrix)
    for row_idx, row in enumerate(r):
        gene_a = high_std_genes[row_idx]
        for col_idx, pcc in enumerate(row):
            # Skip duplicate edges. We only save one copy of each edge.
            if col_idx <= row_idx or pcc == 1 or pcc < 0.4:
                continue
            # Write out gene information.
            gene_b = high_std_genes[col_idx]
            edge_dct[(gene_a, gene_b)] = abs(pcc)
    # Get the top one million edge weights.
    # max_num_edges = int(len(edge_dct) * 0.05)
    sorted_edge_dct = sorted(edge_dct.items(), key=operator.itemgetter(1),
        reverse=True)[:int(1e6)]

    #TODO: Find p-value of worst gene pair.
    # gene_a, gene_b = sorted_edge_dct[-1][0]
    # print len(sorted_edge_dct)#, p[high_std_genes.index(gene_a)][high_std_genes.index(gene_b)]
    return sorted_edge_dct

def write_sorted_edge_dct(sorted_edge_dct, folder_name):
    '''
    Write the sorted edges out to file.
    '''
    out = open('./data/%s_data/high_std_network.txt' % folder_name, 'w')
    for (gene_a, gene_b), edge_weight in sorted_edge_dct:
        out.write('%s\t%s\t%f\n' % (gene_a, gene_b, edge_weight))
    out.close()

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga' % sys.argv[0]
        exit()
    gene_type = sys.argv[1]
    assert gene_type in ['mouse', 'tcga']

    if gene_type == 'mouse':
        folder_list = ['mouse']
    else:
        # If TCGA, make a set of GO dictionaries for each cancer type.
        folder_list = file_operations.get_tcga_disease_list()

    for folder_name in folder_list:
        gene_exp_dct = file_operations.get_gene_expression_dct(folder_name)
        high_std_genes = file_operations.get_high_std_genes(folder_name)

        gene_exp_matrix = create_gene_exp_matrix(high_std_genes, gene_exp_dct)
        sorted_edge_dct = create_sorted_edge_dct(high_std_genes,
            gene_exp_matrix)

        write_sorted_edge_dct(sorted_edge_dct, folder_name)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))