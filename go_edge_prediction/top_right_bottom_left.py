### Author: Edward Huang

from collections import OrderedDict

### Creates the top right and bottom left block of the predicted gene matrix.
### These blocks are simply the transpose of each other, and contain all of the
### gene-GO edge weights.

if __name__ == '__main__':
    # Use this only to remove genes in Sheng's network that do not appear in
    # our coexpression network.
    sampled_genes = []
    f = open('../data/sampled_genes_1_pct.txt', 'r')
    for line in f:
        sampled_genes += [line.strip()]
    f.close()

    gene_index_dct = index_to_name.get_gene_index_dct()
    go_index_dct = index_to_name.get_go_index_dct()

    ### Reading in predicted weights.
    # Rows are genes, columns are GO's.
    gene_go_weight_dct = OrderedDict({}) # Keys are genes. Values are lists of edge weights.
    num_go_terms = len(go_index_dct)
    # Initiate the dictionary.
    for gene in sampled_genes:
        gene_go_weight_dct[gene] = ['0'] * num_go_terms

    f = open('./prediction_data/noisonewAnotationAllNode.txt', 'r')
    for line in f:
        gene, go = map(int, line.split())
        if gene not in gene_index_dct or go not in go_index_dct:
            continue
        # Could be multiple genes per index due to multiple gene aliases.
        genes = gene_index_dct[gene]

        for gene in genes:
            # Record a 1 if there is a relationship between a gene and a GO.
            if gene in sampled_genes:
                gene_go_weight_dct[gene][go] = '1'
    f.close()

    ### Construct gene-GO array.
    gene_GO_weights = gene_go_weight_dct.values()

    out = open('./prediction_data/gene_gene_and_gene_go_weights.txt', 'w')
    # First, write out the top right block.
    f = open('./prediction_data/gene_gene_weights.txt', 'r')
    for i, line in enumerate(f):
        gene_go = gene_GO_weights[i]
        out.write(line.strip() + '\t' + '\t'.join(gene_go) + '\n')
    f.close()

    # Then, transpose the gene-GO block and write it to the bottom left block.
    gene_GO_weights = zip(*gene_GO_weights)
    for row in gene_GO_weights:
        out.write('\t'.join(row) + '\n')
    out.close()