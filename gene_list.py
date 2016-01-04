### Author: Edward Huang

### Creates a newline separated file where each line contains a gene in the
### coexpression matrix.

if __name__ == '__main__':
    f = open('./data/mm_mrsb_log2_expression.tsv', 'r')
    gene_set = set([])
    for i, line in enumerate(f):
        if i == 0:
            continue
        gene_set.add(line.split()[0])
    f.close()

    out = open('./data/all_genes.txt', 'w')
    for gene in gene_set:
        out.write('%s\n' % gene)
    out.close()