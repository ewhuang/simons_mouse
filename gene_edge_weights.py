### Author: Edward Huang

import file_operations
from scipy.stats import pearsonr
import time

### This file generates the base file that we use to create the network with
### GO file and the network without GO file. Reads in the raw data, and
### finds the Pearson correlation coefficient between every pair of genes,
### making an edge in the generated network if the coefficient exceeds a 
### set threshold. Each line is "GENE_1\tGENE_2\tPEARSON_SCORE".
### Run time: 1.8 hours.

def main():
    # Read in the tsv file.
    gene_exp_dct = file_operations.get_gene_expression_dct()
    high_std_genes = file_operations.get_high_std_genes()
    gene_to_index_dct = file_operations.get_gene_to_index_dct(high_std_genes)

    # Calculate the correlations between each pair of genes.
    out = open('./data/high_std_network.txt', 'w')
    for a in range(len(high_std_genes)):
        print '%f%% done...' % (float(a) / len(high_std_genes) * 100)
        gene_a = high_std_genes[a]
        exp_a = gene_exp_dct[gene_a]
        for b in range(a + 1, len(high_std_genes)):
            gene_b = high_std_genes[b]
            exp_b = gene_exp_dct[high_std_genes[b]]
            pcc, p_value = pearsonr(exp_a, exp_b)
            # Ignore genes that have PCC = 1 for the raw network.
            if p_value < 0.001 and pcc < 1.0:
                out.write('%s\t%s\t%f\n' % (gene_to_index_dct[gene_a],
                    gene_to_index_dct[gene_b], abs(pcc)))
    out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))