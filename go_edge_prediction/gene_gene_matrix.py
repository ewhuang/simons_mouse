### Author: Edward Huang


### Creates the top left block in the matrix to be sent to Sheng. Contains the
### gene-gene edge weights, ordered by appearance in original coexpression
### network.

if __name__ == '__main__':
    # First, read in the genes.
    gene_lst = []
    f = open('../data/full_network_genes.txt', 'r')
    for line in f:
        gene = line.strip()
        gene_lst += [gene]
    f.close()

    # Read in the gene-gene edge weights.
    edge_dct = {}
    f = open('../data/full_network.txt', 'r')
    for line in f:
        gene_a, gene_b, weight = line.split()
        gene_a = gene_lst.index(gene_a)
        gene_b = gene_lst.index(gene_a)
        edge_dct[(gene_a, gene_b)] = weight
    f.close()

    out = open('gene_gene_weights.txt', 'w')
    for gene_a in gene_lst:
        for gene_b in gene_lst:
            if (gene_a, gene_b) in edge_dct:
                out.write(edge_dct[(gene_a, gene_b)] + '\t')
            else:
                out.write(edge_dct([gene_b, gene_a]) + '\t')
        out.write('\n')
    out.close()