### Author: Edward Huang

import file_operations
import time

### This script reads in the gene-GO file provided by Sheng. Each line contains
### a gene, followed by a tab, followed by a GO, then a new line. These lines
### denote relationships between each gene listed and the GO.
### Run time: 3 seconds.

def main():
    all_genes = file_operations.get_all_genes()
    gene_index_dct = file_operations.get_gene_index_dct()
    go_index_dct = file_operations.get_go_index_dct()

    out = open('./data/go_edges.txt', 'w')
    fname = './go_edge_prediction/prediction_data/noisonewAnotationAllNode.txt'
    f = open(fname, 'r')
    for line in f:
        gene, go = map(int, line.split())

        # Skip genes that don't appear in our mappings.
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

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))