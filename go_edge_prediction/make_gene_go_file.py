### Author: Edward Huang

import index_to_name

### This script reads in the gene-GO file provided by Sheng. Each line contains
### a gene, followed by a tab, followed by a GO, then a new line. These lines
### denote relationships between each gene listed and the GO.

if __name__ == '__main__':
    f = open('../data/all_genes.txt', 'r')
    all_genes = set([])
    for line in f:
        all_genes.add(line.strip())
    f.close()

    gene_index_dct = index_to_name.get_gene_index_dct()
    go_index_dct = index_to_name.get_go_index_dct()

    out = open('../data/go_edges.txt', 'w')
    f = open('./prediction_data/noisonewAnotationAllNode.txt', 'r')
    for line in f:
        gene, go = map(int, line.split())
        if gene not in gene_index_dct or go not in go_index_dct:
            continue
        genes = gene_index_dct[gene]
        go = go_index_dct[go]
        for gene in genes:
            if gene not in all_genes:
                continue
            out.write(gene + '\t' + go + '\n')
    f.close()
    out.close()