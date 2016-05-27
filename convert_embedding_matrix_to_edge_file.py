### Author: Edward Huang

import file_operations
import time

### This script takes the embedding network, and generates edges in the format:
### gene_a\tgene_b\tweight\n
### in ./data/embedding_edges.txt
### Run time: 178 seconds.

def main():
    # First, retrieve the genes that are available in embedding.
    embedding_genes = file_operations.get_embedding_genes()

    # Mouse.embedding is an nxn matrix, where n is the number of genes in 
    # Mouse.embedding.id.
    f = open('./Sheng/data/network/integrated_network/Mouse.embedding', 'r')
    out = open('./data/embedding_edges.txt', 'w')
    for i, line in enumerate(f):
        row_gene = embedding_genes[i]
        line = line.split()
        assert len(line) == len(embedding_genes)

        # Start at i+1, since we only want upper triangular matrix..
        for j in range(i + 1, len(line)):
            column_gene = embedding_genes[j]
            weight = line[j]
            out.write('%s\t%s\t%s\n' % (row_gene, column_gene, weight))
    out.close()
    f.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))