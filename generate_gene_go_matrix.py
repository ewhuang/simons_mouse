### Author: Edward Huang


### This script creates the gene-GO matrix for predictive edge weights.


if __name__ == '__main__':
    # This is a file where each newline ccontains a gene in the network.
    f = open('./data/all_genes.txt', 'r')
    genes = []
    for line in f:
        gene = line.strip()
        genes += [gene]
    f.close()

    # Read in the file containing all of the edge weights.
    f = open('./data/raw_network.txt', 'r')
    edge_dct = {}
    for line in f:
        gene_a, gene_b, weight = line.split()
        edge_dct[(gene_a, gene_b)] = weight
    f.close()

    out = open('./data/gene_go_weight_matrix.tsv', 'w')
    for gene_a in genes:
        for gene_b in genes:
            if gene_a == gene_b:
                out.write('1\t')
            else:
                if (gene_a, gene_b) in edge_dct:
                    out.write(edge_dct[(gene_a, gene_b)])
                else:
                    out.write(edge_dct[(gene_b, gene_a)])
        out.write('\n')
        exit()

    out.close()
