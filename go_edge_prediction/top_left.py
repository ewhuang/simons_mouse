### Author: Edward Huang

from collections import OrderedDict
import math

### Creates the top left block in the matrix to be sent to Sheng. Contains the
### gene-gene edge weights, ordered by appearance in original coexpression
### network.

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
    f = open('../data/mm_mrsb_log2_expression.tsv', 'r')
    gene_exp_dct = OrderedDict({})
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.split()
        gene, exp_vals = line[0], line[1:]
        exp_vals = [float(val) for val in exp_vals]
        gene_exp_dct[gene] = exp_vals
    f.close()

    # Calculate the correlations between each pair of genes.
    print 'Computing Pearson correlation coefficients...'
    out = open('./prediction_data/gene_gene_weights.txt', 'w')
    genes = gene_exp_dct.keys()
    for a in range(len(genes)):
        print '%f%% done...' % (float(a) / len(genes) * 100)
        gene_a = genes[a]
        exp_a = gene_exp_dct[gene_a]
        for b in range(len(genes)):
            gene_b = genes[b]
            exp_b = gene_exp_dct[gene_b]
            pcc = str(abs(pearsonr(exp_a, exp_b)))
            # If pearson is 1, then make no edge.
            if pcc == 1.0:
                pcc = 0.0
            out.write(pcc + '\t')
        out.write('\n')
    out.close()