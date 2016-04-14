### Author: Edward Huang

import file_operations
from collections import OrderedDict
import math

### This file generates the base file that we use to create the network with
### GO file and the network without GO file. Reads in the raw data, and
### finds the Pearson correlation coefficient between every pair of genes,
### making an edge in the generated network if the coefficient exceeds a 
### set threshold. Each line is "GENE_1\tGENE_2\tPEARSON_SCORE".

pearson_threshold = 0.8

def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

def pearsonr(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff
    return diffprod / math.sqrt(xdiff2 * ydiff2)

if __name__ == '__main__':
    # Read in the tsv file.
    print 'Reading in raw data...'
    gene_exp_dct = file_operations.get_gene_expression_dct()

    # Calculate the correlations between each pair of genes.
    print 'Computing Pearson correlation coefficients...'
    out = open('./data/raw_network.txt', 'w')
    out2 = open('./data/full_network.txt', 'w')
    genes = gene_exp_dct.keys()
    for a in range(len(genes)):
        print '%f%% done...' % (float(a) / len(genes) * 100)
        gene_a = genes[a]
        exp_a = gene_exp_dct[gene_a]
        for b in range(a + 1, len(genes)):
            gene_b = genes[b]
            exp_b = gene_exp_dct[genes[b]]
            pcc = abs(pearsonr(exp_a, exp_b))
            pcc = str(pcc)
            # Always write to the full network.
            out2.write(gene_a + '\t' + gene_b + '\t' + pcc + '\n')
            # Ignore genes that have PCC = 1 for the raw network.
            if pearson_threshold < pcc < 1.0:
                out.write(gene_a + '\t' + gene_b + '\t' + pcc + '\n')
    out.close()
    out2.close()