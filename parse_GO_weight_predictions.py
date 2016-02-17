### Author: Edward Huang

import math
import bisect

### This script takes in the matrix created by the GO hierarchy method in order
### to find the edge weights between GO and genes. We take the NUM_TOP_WEIGHTS
### weights from the matrix and save them as edges in our network. Writes out 
### to a file, predicted_go_edges.txt.

# NUM_TOP_WEIGHTS = int(1e06)
NUM_TOP_WEIGHTS = 33284

if __name__ == '__main__':
    # Initialize the dictionary with all of the genes in the coexpression
    # matrix.
    f = open('./data/all_genes.txt', 'r')
    all_genes = set([])
    for line in f:
        all_genes.add(line.strip())
    f.close()

    # Find all MGI mappings from http://www.informatics.jax.org/
    # Use Batch query.
    mgi_to_ensembl_dct = {}
    f = open('./go_hierarchy/mgi_to_ensembl.txt', 'r')
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
        if ENSMUSG not in all_genes:
            continue
        mgi_id = line[0]
        if mgi_id in mgi_to_ensembl_dct:
            mgi_to_ensembl_dct[mgi_id] += [ENSMUSG]
        else:
            mgi_to_ensembl_dct[mgi_id] = [ENSMUSG]
    f.close()

    # Keys are the indices in the edge weight matrix, values are the genes.
    gene_index_dct = {}
    f = open('./go_hierarchy/noisogeneHash.txt', 'r')
    for line in f:
        mgi_id, row = line.split()
        if mgi_id not in mgi_to_ensembl_dct:
            continue
        # Subtract 1 to change to list indices.
        row_index = int(row) - 1
        gene_index_dct[row_index] = mgi_to_ensembl_dct[mgi_id]
    f.close()

    # Find the GO terms to GO names.
    # Retrieved data from http://geneontology.org/page/download-annotations
    # go_id_to_name_dct = {}
    # f = open('./go_hierarchy/go_to_name.txt', 'r')
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
    # f = open('./go_hierarchy/noisogoHash.txt', 'r')
    # for line in f:
    #     go_id, index = line.split()
    #     # Subtract 1 to change to list indices.
    #     go_index_dct[int(index) - 1] = go_id_to_name_dct[go_id]
    # f.close()

    # Run through the matrix file a first time, and gather the top edge weights.
    top_weights = [0] * NUM_TOP_WEIGHTS
    f = open('./go_hierarchy/Mouse_final_Score_matrix.txt', 'r')
    for i, line in enumerate(f):
        # Skip a row if the corresponding gene does not appear in our network.
        if i not in gene_index_dct:
            continue
        line = [float(weight) for weight in line.split()]
        for j, weight in enumerate(line):
            if j < i or weight <= 0:
                continue
            # Insert the weight into our sorted list of weights.
            bisect.insort(top_weights, weight)
            # Pop the first element of the list to retain size of the list.
            top_weights.pop(0)
    f.close()

    assert (len(top_weights) == NUM_TOP_WEIGHTS)
    smallest_top_weight = top_weights[0]

    # Run through the matrix file a second time, and write out the edges.
    out = open('./data/predicted_go_edges.txt', 'w')
    # Rows are genes, columns are GO's.
    f = open('./go_hierarchy/Mouse_final_Score_matrix.txt', 'r')
    for i, line in enumerate(f):
        # Skip a row if the corresponding gene does not appear in our network.
        if i not in gene_index_dct:
            continue
        line = [float(weight) for weight in line.split()]
        for j, weight in enumerate(line):
            if weight < smallest_top_weight:
                continue
            # Somtimes an MGI ID will correspond to multiple ENSMUSG's.
            for gene in gene_index_dct[i]:
                out.write('%s\t%s\n' % (gene, j))#go_index_dct[j]))
    f.close()
    out.close()