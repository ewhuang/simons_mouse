### Author: Edward Huang

import file_operations
import numpy as np
import operator
from scipy.stats import betai
import sys
import time

### Writes top 1% of edges based on Pearson correlation coefficient out to file.
### Run time: 50s for mouse, 100s for TCGA, 12min for all TCGA.

def get_high_std_genes(data_type):
    '''
    Retrieves the list of genes that have high standard deviations across their
    gene expression vectors.
    '''
    high_std_genes = []
    f = open('./data/%s_data/high_std_genes.txt' % data_type, 'r')
    for line in f:
        gene = line.strip()
        assert 'ENSMUSG' in gene or 'ENSG' in gene
        high_std_genes += [gene]
    f.close()
    return high_std_genes

def write_edges(gene_type):
    '''
    Given the 2D gene expression matrix, compute the pairwise Pearson
    coefficients, and then make the edge dictionary. Returns the edges sorted
    by decreasing weight. Take only the top one million edges.
    '''
    gene_exp_dct = file_operations.get_gene_expression_dct(gene_type)
    high_std_genes = get_high_std_genes(gene_type)

    # Convert dictionary into 2D list.
    gene_exp_matrix = []
    for gene in high_std_genes:
        gene_exp_matrix += [gene_exp_dct[gene]]

    # Compute pairwise Pearson correlation coefficients.
    r = np.corrcoef(np.array(gene_exp_matrix))
    r = np.triu(r, 1) # Get the upper triangular matrix.
    r[r == 1] = 0 # Remove identical pairs as edges.

    num_genes = len(high_std_genes)
    num_edges = int(0.01 * num_genes * (num_genes - 1) / 2.0)
    x_idx, y_idx = np.unravel_index(np.argsort(r.ravel())[-num_edges:], r.shape)

    edge_list, network_gene_set = [], set([])
    for x, y in zip(x_idx, y_idx):
        gene_x, gene_y = high_std_genes[x], high_std_genes[y]
        edge_list += ['%s\t%s\t%g' % (gene_x, gene_y, r[x][y])]
        network_gene_set.add(gene_x)
        network_gene_set.add(gene_y)

    # Write edges out to file.
    out = open('./data/%s_data/coexpression_network.txt' % gene_type, 'w')
    out.write('\n'.join(edge_list[::-1]))
    out.close()

    # Also write genes out to file.
    gene_out = open('./data/%s_data/network_genes.txt' % gene_type, 'w')
    gene_out.write('\n'.join(network_gene_set))
    gene_out.close()

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga/all' % sys.argv[0]
        exit()
    data_type = sys.argv[1]
    assert data_type in ['mouse', 'tcga', 'all']

    # Convert 'all' to list of all TCGA sub-cancers.
    gene_type_list = [data_type]
    if data_type == 'all':
        gene_type_list = file_operations.get_tcga_disease_list()

    for gene_type in gene_type_list:
        write_edges(gene_type)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))