### Author: Edward Huang

### Creates a newline separated file where each line contains a gene in the
### coexpression matrix.

if __name__ == '__main__':
    f = open('./data/mm_mrsb_log2_expression.tsv', 'r')
    gene_lst = []
    for i, line in enumerate(f):
        if i == 0:
            continue
        if line.split()[0] in gene_lst:
            print 'already here?'
        gene_lst += [line.split()[0]]
    f.close()

    out = open('./data/all_genes.txt', 'w')
    for gene in gene_lst:
        out.write('%s\n' % gene)
    out.close()