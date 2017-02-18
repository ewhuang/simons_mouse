### Author: Edward Huang

import file_operations
import matplotlib
import numpy as np
import os
import sys
import time
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab

### This script plots a histogram of the number of genes vs. standard deviation
### of their gene expression vectors. It also writes out to file the genes
### that have high standard deviation.
### Run time: 18 seconds for mouse, 110 seconds for TCGA.

# Plotting GO enrichment histograms.
def plot_histogram(gene_type, std_list):
    # Plot 50 bins.
    bins = np.linspace(0, 5, 51)
    matplotlib.pyplot.hist(std_list, bins, alpha=0.25)

    plt.xlabel('standard deviation')
    plt.ylabel('number of genes')
    plt.title('number of genes vs. STD')
    plt.show()
    pylab.savefig('./data/%s_data/gene_std_histogram.png' % gene_type)
    plt.close()

# Writing out genes with high standard deviation.
def write_genes_to_file(gene_type, high_std_genes):
    folder = './data/%s_data' % gene_type
    if not os.path.exists(folder):
        os.makedirs(folder)
    out = open('%s/high_std_genes.txt' % folder, 'w')
    out.write('\n'.join(high_std_genes))
    out.close()

def get_std_threshold(high_std_gene_dct):
    '''
    Finds a threshold of standard deviation at which to remove genes.
    '''
    # Initialize the list of bins.
    std_range_dct = {}
    std_range = [x / 10.0 for x in range(11)]
    for std in std_range:
        std_range_dct[std] = 0

    # Get a dictionary mapping STD's to the gene counts.
    for gene in high_std_gene_dct:
        gene_std = high_std_gene_dct[gene]
        if gene_std < 1.0:
            std_range_dct[np.floor(gene_std * 10.0) / 10.0] += 1

    previous_count = std_range_dct[std_range[0]]
    for std in std_range[1:]:
        # Stop when we have reached a count lower than the previous count, with
        # at > 1000 genes in the bin.
        current_count = std_range_dct[std]
        if 1000 < current_count and current_count < previous_count:
            return std
        previous_count = current_count
    exit()

def cull_low_std_genes(high_std_gene_dct, std_threshold):
    '''
    With the standard deviation threshold, remove genes that are below it.
    '''
    high_std_genes = []
    for gene in high_std_gene_dct:
        if high_std_gene_dct[gene] >= std_threshold:
            high_std_genes += [gene]
    return high_std_genes

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga/all' % sys.argv[0]
        exit()
    data_type = sys.argv[1]
    assert data_type in ['mouse', 'tcga', 'all']

    # Convert 'all' to list of all TCGA sub-cancers.
    gene_type_list = [data_type]
    if data_type == 'all':
        gene_type_list = file_operations.get_tcga_list()

    for gene_type in gene_type_list:
        gene_expression_dct = file_operations.get_gene_expression_dct(gene_type)

        # Compute standard deviations for each gene.
        std_list, high_std_gene_dct = [], {}
        for gene in gene_expression_dct:
            # Transform back to non-log.
            gene_exp_vector = [pow(2, val) for val in gene_expression_dct[gene]]

            # Compute standard deviation.
            gene_std = np.std(gene_exp_vector)

            # 4.5e-15 happens when we have the same values for all samples.
            if gene_std < 4.5e-15:
                continue
            std_list += [gene_std]
            high_std_gene_dct[gene] = gene_std

        # TODO: automated or hard-coded threshold?
        std_threshold = get_std_threshold(high_std_gene_dct)
        # if gene_type == 'mouse':
            # std_threshold = 0.1
        # elif gene_type == 'tcga':
            # std_threshold = 0.5
        # else:
        #     std_threshold = 0.5

        high_std_genes = cull_low_std_genes(high_std_gene_dct, std_threshold)
        print 'Number of high standard deviation genes: ', len(high_std_genes),
        print std_threshold, gene_type
        write_genes_to_file(gene_type, high_std_genes)
        plot_histogram(gene_type, std_list)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))