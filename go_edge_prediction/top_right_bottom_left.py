### Author: Edward Huang

from collections import OrderedDict
import math

### Creates the top right and bottom left block of the predicted gene matrix.
### These blocks are simply the transpose of each other, and contain all of the
### gene-GO edge weights.


if __name__ == '__main__':
    # Use this only to remove genes in Sheng's network that do not appear in
    # our coexpression network.
    sampled_genes_set = set([])
    f = open('../data/sampled_genes_1_pct.txt', 'r')
    for line in f:
        sampled_genes_set.add(line.strip())
    f.close()

    # Read in the tsv file.
    print 'Reading in raw data...'
    f = open('../data/mm_mrsb_log2_expression.tsv', 'r')
    sampled_genes = []
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.split()
        gene, exp_vals = line[0], line[1:]
        if gene not in sampled_genes_set:
            continue
        sampled_genes += [gene]
    f.close()

    ### GENE MAPPINGS.
    # Find all MGI mappings from http://www.informatics.jax.org/
    # Use Batch query.
    mgi_to_ensembl_dct = {}
    f = open('./prediction_data/mgi_to_ensembl.txt', 'r')
    for i, line in enumerate(f):
        # Skip the header line.
        if i == 0:
            continue
        line = line.split()
        if len(line) != 5 or (line[1] != 'current' and line[1] != 'old'):
            continue
        ENSMUSG = line[4]
        # Skip if a gene isn't in our coexpression network, or if it isn't up
        # to date.
        if ENSMUSG not in sampled_genes:
            continue
        mgi_id = line[0]
        if mgi_id in mgi_to_ensembl_dct:
            mgi_to_ensembl_dct[mgi_id] += [ENSMUSG]
        else:
            mgi_to_ensembl_dct[mgi_id] = [ENSMUSG]
    f.close()

    # Keys are the indices in the edge weight matrix, values are the genes.
    gene_index_dct = {}
    f = open('./prediction_data/noisogeneHash.txt', 'r')
    for line in f:
        mgi_id, row = line.split()
        if mgi_id not in mgi_to_ensembl_dct:
            continue
        # Subtract 1 to change to list indices.
        row_index = int(row) - 1
        gene_index_dct[row_index] = mgi_to_ensembl_dct[mgi_id]
    f.close()

    ### GO MAPPINGS.
    go_id_to_name_dct = {}
    f = open('./prediction_data/go_to_name.txt', 'r')
    while True:
        line = f.readline()
        if line == '':
            break
        if line.strip() == '[Term]':
            go_id = f.readline().split()[1]
            go_name = '_'.join(f.readline()[len('name: '):].split())
            go_id_to_name_dct[go_id.lower()] = go_name
            f.readline()
            next = f.readline()
            while 'alt_id' in next:
                go_id_to_name_dct[next.split()[1].lower()] = go_name
                next = f.readline()
    f.close()

    # Keys are the GO ids, values are the indices in the edge weight matrix.
    go_index_dct = {}
    f = open('./prediction_data/noisogoHash.txt', 'r')
    for line in f:
        go_id, index = line.split()
        # Subtract 1 to change to list indices.
        go_index_dct[int(index) - 1] = go_id_to_name_dct[go_id]
    f.close()

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
                gene_go_weight_dct[gene][go] = 1
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