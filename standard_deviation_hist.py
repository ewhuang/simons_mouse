### Author: Edward Huang

import file_operations
import matplotlib
import numpy as np
import operator
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
def plot_histogram(data_type, std_list):
    # Plot 50 bins.
    bins = np.linspace(0, 5, 50)
    matplotlib.pyplot.hist(std_list, bins, alpha=0.25)

    plt.xlabel('standard deviation')
    plt.ylabel('number of genes')
    plt.title('number of genes vs. STD value')
    plt.show()
    pylab.savefig('./data/%s_data/gene_std_histogram.png' % data_type)
    plt.close()

# Writing out genes with high standard deviation.
def write_genes_to_file(data_type, high_std_genes):
    out = open('./data/%s_data/high_std_genes.txt' % data_type, 'w')
    for gene in high_std_genes:
        out.write(gene + '\n')
    out.close()

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga' % sys.argv[0]
        exit()
    category = sys.argv[1]
    assert category in ['mouse', 'tcga']

    if category == 'mouse':
        data_type_list = ['mouse']
    elif category == 'tcga':
         data_type_list = file_operations.get_tcga_disease_list()

    # Loop through the data type list. data_type can be 'mouse' or any of the
    # TCGA cancers.
    for data_type in data_type_list:
        gene_expression_dct = file_operations.get_gene_expression_dct(data_type)
        if data_type == 'mouse':
            embedding_genes = file_operations.get_embedding_genes()

        # Compute standard deviations for each gene.
        std_list, high_std_gene_dct = [], {}
        for gene in gene_expression_dct:
            if data_type == 'mouse' and gene not in embedding_genes:
                continue
            gene_exp_vector = gene_expression_dct[gene]

            # Transform back to non-log.
            gene_exp_vector = [pow(2, val) for val in gene_exp_vector]

            # Compute standard deviation.
            gene_std = np.std(gene_exp_vector)

            # 4.5e-15 happens when the genes have the same expression values for
            # all samples.
            if gene_std < 4.5e-15:
                continue
            
            high_std_gene_dct[gene] = gene_std
            std_list += [gene_std]

        # Get the 15k genes with the highest standard deviations.
        high_std_genes = sorted(high_std_gene_dct.items(),
            key=operator.itemgetter(1), reverse=True)[:15000]
        # Get just the gene names.
        high_std_genes = [pair[0] for pair in high_std_genes]

        write_genes_to_file(data_type, high_std_genes)
        plot_histogram(data_type, std_list)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))