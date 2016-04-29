### Author: Edward Huang

import file_operations
import matplotlib.pyplot as plt
import numpy as np
import time

### This script plots a histogram of the number of genes vs. standard deviation
### of their gene expression vectors. It also writes out to file the genes
### that have high standard deviation.
### Run time: 5 seconds

def main():
    # Compute standard deviations for each gene
    std_list, high_std_genes = [], []
    gene_exp_dct = file_operations.get_gene_expression_dct()
    for gene in gene_exp_dct:
        gene_exp_vector = gene_exp_dct[gene]
        # Transform back to non-log.
        gene_exp_vector = [pow(2, val) for val in gene_exp_vector]
        gene_std = np.std(gene_exp_vector)
        if gene_std < 4.5e-15:
            continue
        elif gene_std > 0.1:
            high_std_genes += [gene]
        std_list += [gene_std]

    # Writing out genes with high standard deviation.
    out = open('./data/high_std_genes.txt', 'w')
    for gene in high_std_genes:
        out.write(gene + '\n')
    out.close()

    # Plotting GO enrichment histograms.
    bins = np.linspace(0, 5, 50)
    plt.hist(std_list, bins, alpha=0.5, label='GO')

    plt.xlabel('standard deviation')
    plt.ylabel('number of genes')
    plt.title('Number of genes for each standard deviation, 90 samples')
    plt.show()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))