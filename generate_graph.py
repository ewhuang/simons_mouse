### Author: Edward Huang

from collections import OrderedDict
import math

### This file generates the base file that we use to create the network with
### GO file and the network without GO file. Reads in the raw data, and
### finds the Pearson correlation coefficient between every pair of genes,
### making an edge in the generated network if the coefficient exceeds a 
### set threshold. Each line is "GENE_1\tGENE_2\tPEARSON_SCORE". Each edge
### is written twice to create an undirected edge.

pearson_threshold = 0.85

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
    f = open('./data/mm_mrsb_log2_expression.tsv', 'r')
    gene_exp_dct = OrderedDict({})
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.split()
        gene, exp_vals = line[0], line[1:]
        exp_vals = [float(val) for val in exp_vals]
        gene_exp_dct[gene] = exp_vals
    f.close()

    i = 0
    # Calculate the correlations between each pair of genes.
    print 'Computing Pearson correlation coefficients...'
    out = open('./data/mouse_network.txt', 'w')
    genes = gene_exp_dct.keys()
    for a in range(len(genes)):
        gene_a = genes[a]
        exp_a = gene_exp_dct[gene_a]
        for b in range(a + 1, len(genes)):
            gene_b = genes[b]
            exp_b = gene_exp_dct[genes[b]]
            rcc = abs(pearsonr(exp_a, exp_b))
            # Ignore genes that have PCC = 1.
            if pearson_threshold < rcc < 1.0:
                rcc = str(rcc)
                out.write(gene_a + '\t' + gene_b + '\t' + rcc + '\n')
                out.write(gene_b + '\t' + gene_a + '\t' + rcc + '\n')
    out.close()