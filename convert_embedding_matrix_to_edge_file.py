### Author: Edward Huang

import time

### This script takes the embedding network, and generates edges in the format:
### gene_a\tgene_b\tweight\n
### in ./data/embedding_edges.txt
### Run time: 178 seconds.

data_folder = './Sheng/data/network/integrated_network/'

# Mouse.embedding.id contains a file of ENSMUSG ID's separated by newlines.
def get_embedding_genes():
    embedding_genes = []
    f = open(data_folder + 'Mouse.embedding.id', 'r')
    for i, line in enumerate(f):
        ensmusg_id = line.strip()
        assert 'ENSMUSG' in ensmusg_id
        embedding_genes += [ensmusg_id]
    f.close()
    return embedding_genes

def main():
    embedding_genes = get_embedding_genes()

    # Mouse.embedding is an nxn matrix, where n is the number of genes in 
    # Mouse.embedding.id.
    f = open(data_folder + 'Mouse.embedding', 'r')
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