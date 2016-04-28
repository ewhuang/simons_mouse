### Author: Edward Huang

import file_operations
from scipy.stats import pearsonr
import time

### This file generates the base file that we use to create the network with
### GO file and the network without GO file. Reads in the raw data, and
### finds the Pearson correlation coefficient between every pair of genes,
### making an edge in the generated network if the coefficient exceeds a 
### set threshold. Each line is "GENE_1\tGENE_2\tPEARSON_SCORE".
### Run time:

def main():
    # Read in the tsv file.
    gene_exp_dct = file_operations.get_gene_expression_dct()
    genes = file_operations.get_ppi_go_high_std_genes()
    genes_to_indices = file_operations.map_genes_to_indices(genes)

    # Calculate the correlations between each pair of genes.
    out = open('./data/raw_network.txt', 'w')
    for a in range(len(genes)):
        print '%f%% done...' % (float(a) / len(genes) * 100)
        gene_a = genes[a]
        exp_a = gene_exp_dct[gene_a]
        for b in range(a + 1, len(genes)):
            gene_b = genes[b]
            exp_b = gene_exp_dct[genes[b]]
            pcc, p_value = pearsonr(exp_a, exp_b)
            # Ignore genes that have PCC = 1 for the raw network.
            if p_value < 0.005 and pcc < 1.0:
                out.write('%d\t%d\t%f\n' % (genes_to_indices[gene_a],
                    genes_to_indices[gene_b], abs(pcc)))
    out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))