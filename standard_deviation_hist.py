### Author: Edward Huang

import file_operations
import math
import matplotlib.pyplot as plt
import numpy
import time

### This script plots a histogram of the number of genes vs. standard deviation
### of their gene expression vectors.

def write_ppi_and_go_genes():
    ppi_genes = set([])
    f = open('./data/embedding_full_network.txt', 'r')
    for line in f:
        gene_a, gene_b, weight = line.split()
        ppi_genes.add(gene_a)
        ppi_genes.add(gene_b)
    f.close()

    go_genes = set([])
    f = open('./data/go_edges.txt', 'r')
    for line in f:
        gene, go_label = line.split()
        go_genes.add(gene)
    f.close()

    out = open('./data/ppi_and_go_genes.txt', 'w')
    for gene in go_genes.intersection(ppi_genes):
        out.write(gene + '\n')
    out.close()

# Returns standard deviation of a list.
def std(x):
    std = []
    for value in x:
        std.append(pow((value - (sum(x)/len(x))), 2))
    stddev = math.sqrt(sum(std)/len(std))
    mean = (sum(x)/len(x))
    return float(stddev)

def main():
    print 'Currently not writing out ppi_and_go_genes.txt...'
    # write_ppi_and_go_genes()

    # Read in the genes from the PPI enriched network and GO annotated genes.
    ppi_and_go_genes = []
    f = open('./data/ppi_and_go_genes.txt', 'r')
    for line in f:
        gene = line.strip()
        ppi_and_go_genes += [gene]
    f.close()

    # Compute standard deviations for each gene
    std_list = []
    flat_gene = []
    high_std_genes = []
    gene_exp_dct = file_operations.get_gene_expression_dct()
    for gene in gene_exp_dct:
        # gene_std = -math.log(std(gene_exp_dct[gene]), 10)
        gene_exp_vector = gene_exp_dct[gene]
        # Transform back to non-log.
        gene_exp_vector = [pow(2, val) for val in gene_exp_vector]
        gene_std = std(gene_exp_vector)
        if gene_std < 4.5e-15:
            if flat_gene == []:
                flat_gene = gene_exp_vector
            else:
                assert flat_gene == gene_exp_vector
            continue
        elif gene_std > 0.1:
            high_std_genes += [gene]
        std_list += [gene_std]

    # Writing out genes with high standard deviation.
    out = open('./data/ppi_and_go_genes_high_std.txt', 'w')
    for gene in high_std_genes:
        if gene in ppi_and_go_genes:
            out.write(gene + '\n')
    out.close()

    # Plotting GO enrichment histograms.
    bins = numpy.linspace(0, 5, 50)
    plt.hist(std_list, bins, alpha=0.5, label='GO')

    plt.xlabel('standard deviation')
    plt.ylabel('number of genes')
    plt.title('Number of genes for each standard deviation, 90 samples')
    plt.show()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))