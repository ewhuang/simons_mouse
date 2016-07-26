### Author: Edward Huang

import file_operations
import matplotlib
import numpy as np
import sys
import time

### This script plots a histogram of the number of genes vs. standard deviation
### of their gene expression vectors. It also writes out to file the genes
### that have high standard deviation.
### Run time: 12 seconds for mouse, 90 seconds for TCGA.

# Chosen from examining the histogram.
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab

# Plotting GO enrichment histograms.
def plot_histogram(std_list):
    bins = np.linspace(0, 5, 50)
    matplotlib.pyplot.hist(std_list, bins, alpha=0.5, label='GO')

    plt.xlabel('standard deviation')
    plt.ylabel('number of genes')
    plt.title('Number of genes for each standard deviation, 90 samples')
    plt.show()
    pylab.savefig('./data/%s_gene_standard_deviation_hist.png' % data_type)

# Writing out genes with high standard deviation.
def write_genes_to_file(high_std_genes):
    out = open('./data/%s_high_std_genes.txt' % data_type, 'w')
    for gene in high_std_genes:
        out.write(gene + '\n')
    out.close()

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga' % sys.argv[0]
        exit()
    global data_type
    data_type = sys.argv[1]
    assert data_type in ['mouse', 'tcga']

    gene_expression_dct = file_operations.get_gene_expression_dct(data_type)
    if data_type == 'mouse':
        embedding_genes = file_operations.get_embedding_genes()
        standard_deviation_lower_bound = 0.1
    else:
        # This variable is set by looking at the histogram, and finding where
        # the steepest cutoff is.
        standard_deviation_lower_bound = 0.6

    # Compute standard deviations for each gene. std_list is used for plotting.
    std_list, high_std_genes = [], []
    for gene in gene_expression_dct:
        if data_type == 'mouse' and gene not in embedding_genes:
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
        
        if gene_std > standard_deviation_lower_bound:
            high_std_genes += [gene]
        std_list += [gene_std]

    write_genes_to_file(high_std_genes)
    plot_histogram(std_list)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))