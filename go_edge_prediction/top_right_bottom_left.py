### Author: Edward Huang

from collections import OrderedDict
import math

### Creates the top right and bottom left block of the predicted gene matrix.
### These blocks are simply the transpose of each other, and contain all of the
### gene-GO edge weights.


if __name__ == '__main__':
    # Use this only to remove genes in Sheng's network that do not appear in
    # our coexpression network.
    f = open('../data/all_genes.txt', 'r')
    gene_list = []
    for line in f:
        gene_list += [line.strip()]
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
        if ENSMUSG not in gene_list:
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
    # Find the GO terms to GO names.
    # Retrieved data from http://geneontology.org/page/download-annotations
    # go_id_to_name_dct = {}
    # f = open('./prediction_data/go_to_name.txt', 'r')
    # while True:
    #     line = f.readline()
    #     if line == '':
    #         break
    #     if line.strip() == '[Term]':
    #         go_id = f.readline().split()[1]
    #         go_name = '_'.join(f.readline()[len('name: '):].split())
    #         go_id_to_name_dct[go_id.lower()] = go_name
    #         f.readline()
    #         next = f.readline()
    #         while 'alt_id' in next:
    #             go_id_to_name_dct[next.split()[1].lower()] = go_name
    #             next = f.readline()
    # f.close()
    # # Keys are the GO ids, values are the indices in the edge weight matrix.
    # go_index_dct = {}
    # f = open('./prediction_data/noisogoHash.txt', 'r')
    # for line in f:
    #     go_id, index = line.split()
    #     # Subtract 1 to change to list indices.
    #     go_index_dct[int(index) - 1] = go_id_to_name_dct[go_id]
    # f.close()

    ### Reading in predicted weights.
    # Rows are genes, columns are GO's.
    gene_go_weight_dct = {} # Keys are genes. Values are lists of edge weights.
    num_go_terms = 0
    f = open('./prediction_data/Mouse_final_Score_matrix.txt', 'r')
    for i, line in enumerate(f):
        # Make a row of zeroes if the corresponding gene does not appear in our
        # network.
        if i not in gene_index_dct:
            continue

        line = line.split()

        for gene in gene_index_dct[i]:
            gene_go_weight_dct[gene] = line
        
        if num_go_terms == 0:
            num_go_terms = len(line)
        
        assert num_go_terms == len(line)


        # for j, weight in enumerate(line):
        #     if weight < smallest_top_weight:
        #         continue
        #     # Somtimes an MGI ID will correspond to multiple ENSMUSG's.
        #     for gene in gene_index_dct[i]:
        #         out.write('%s\t%s\n' % (gene, j))#go_index_dct[j]))
    f.close()

    ### Construct gene-GO array.
    gene_GO_weights = []
    for curr_gene in gene_list:
        if curr_gene in gene_go_weight_dct:
            gene_GO_weights += [gene_go_weight_dct[curr_gene]]
        else:
            gene_GO_weights += [['0'] * num_go_terms]

    out = open('./prediction_data/gene_gene_and_gene_go_weights.txt', 'w')
    # First, write out the top right block.
    f = open('./prediction_data/gene_gene_weights.txt', 'r')
    for i, line in enumerate(f):
        gene_gene = line.strip()
        gene_go = gene_GO_weights[i]
        out.write(gene_gene + '\t' + '\t'.join(gene_go) + '\n')
    f.close()

    # Then, transpose the gene-GO block and write it to the bottom left block.
    gene_GO_weights = zip(*gene_GO_weights)
    for row in gene_GO_weights:
        out.write('\t'.join(row) + '\n')
    out.close()