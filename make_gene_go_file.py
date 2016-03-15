### Author: Edward Huang


### This script reads in the gene-GO file provided by Sheng. Each line contains
### a gene, followed by a tab, followed by a GO, then a new line. These lines
### denote relationships between each gene listed and the GO.


if __name__ == '__main__':

    f = open('./data/all_genes.txt', 'r')
    all_genes = set([])
    for line in f:
        all_genes.add(line.strip())
    f.close()

    mgi_to_ensembl_dct = {}
    f = open('./go_edge_prediction/prediction_data/mgi_to_ensembl.txt', 'r')
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
    f = open('./go_edge_prediction/prediction_data/noisogeneHash.txt', 'r')
    for line in f:
        mgi_id, row = line.split()
        if mgi_id not in mgi_to_ensembl_dct:
            continue
        # Subtract 1 to change to list indices.
        row_index = int(row) - 1
        gene_index_dct[row_index] = mgi_to_ensembl_dct[mgi_id]
    f.close()

    go_id_to_name_dct = {}
    f = open('./go_edge_prediction/prediction_data/go_to_name.txt', 'r')
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
    f = open('./go_edge_prediction/prediction_data/noisogoHash.txt', 'r')
    for line in f:
        go_id, index = line.split()
        # Subtract 1 to change to list indices.
        go_index_dct[int(index) - 1] = go_id_to_name_dct[go_id]
    f.close()

    out = open('./data/go_edges.txt', 'w')
    f = open('./go_edge_prediction/prediction_data/noisonewAnotationAllNode.txt', 'r')
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