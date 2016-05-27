### Author: Edward Huang

import file_operations
import matplotlib.pyplot as plt
import numpy as np
import time

### This script plots a histogram of the number of genes vs. standard deviation
### of their gene expression vectors. It also writes out to file the genes
### that have high standard deviation.
### Run time: 12 seconds

# Chosen from examining the histogram.
HIGH_STANDARD_DEVIATION = 0.1

# Plotting GO enrichment histograms.
def plot_histogram(std_list):
    bins = np.linspace(0, 5, 50)
    plt.hist(std_list, bins, alpha=0.5, label='GO')

    plt.xlabel('standard deviation')
    plt.ylabel('number of genes')
    plt.title('Number of genes for each standard deviation, 90 samples')
    plt.show()

# Writing out genes with high standard deviation.
def write_genes_to_file(high_std_genes):
    out = open('./data/high_std_genes.txt', 'w')
    for gene in high_std_genes:
        out.write(gene + '\n')
    out.close()

def main():
    gene_expression_dct = file_operations.get_gene_expression_dct()
    embedding_genes = file_operations.get_embedding_genes()

    # Compute standard deviations for each gene. std_list is used for plotting.
    std_list, high_std_genes = [], []
    for gene in gene_expression_dct:
        if gene not in embedding_genes:
            continue
        gene_exp_vector = gene_expression_dct[gene]

        # Transform back to non-log.
        gene_exp_vector = [pow(2, val) for val in gene_exp_vector]

        # Compute standard deviation.
        gene_std = np.std(gene_exp_vector)

        # 4.5e-15 happens when the genes have the same expression values for all
        # samples.
        if gene_std < 4.5e-15:
            continue
        
        if gene_std > HIGH_STANDARD_DEVIATION:
            high_std_genes += [gene]
        std_list += [gene_std]

    write_genes_to_file(high_std_genes)

    print 'Not plotting...'
    # plot_histogram(std_list)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))