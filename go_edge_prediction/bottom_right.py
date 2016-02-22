### Author: Edward Huang

from collections import OrderedDict
import math

### Creates the bottom right block of the predicted gene matrix. This block
### contains the GO-GO edge weights.

if __name__ == '__main__':
    ### GO MAPPINGS.
    # Find the GO terms to GO names.
    # Retrieved data from http://geneontology.org/page/download-annotations
    # Keys are the GO terms, values are lists of other GO terms the keys have
    # an edge to.
    go_go_edge_dct = {}
    f = open('./prediction_data/go_to_name.txt', 'r')
    while True:
        line = f.readline()
        if line == '':
            break
        if line.strip() == '[Term]':
            # Keep track of the current term's GO id's.
            current_go_id_list = []

            go_id = f.readline().split()[1].lower()
            current_go_id_list += [go_id]
            go_name = f.readline()
            assert ('name:' in go_name)
            namespace = f.readline()
            assert ('namespace' in namespace)
            next = f.readline()
            while 'alt_id' in next:
                alt_id = next.split()[1].lower()
                current_go_id_list += [alt_id]
                next = f.readline()
            # Keep reading until we see an is_a relationship.
            neighbors = []
            while next.strip() != '':
                if 'is_a:' in next or 'part_of' in next:
                    if 'is_a' in next:
                        neighbor_GO = next.split()[1].lower()
                    else:
                        neighbor_GO = next.split()[2].lower()
                    neighbors += [neighbor_GO]
                next = f.readline()


            for go_id in current_go_id_list:
                if go_id in go_go_edge_dct:
                    go_go_edge_dct[go_id] += neighbors
                else:
                    go_go_edge_dct[go_id] = neighbors[:]

            for neighbor_id in neighbors:
                if neighbor_id in go_go_edge_dct:
                    go_go_edge_dct[neighbor_id] += current_go_id_list
                else:
                    go_go_edge_dct[neighbor_id] = current_go_id_list[:]
    f.close()

    # Keys are the GO ids, values are the indices in the edge weight matrix.
    go_index_dct = {}
    f = open('./prediction_data/noisogoHash.txt', 'r')
    for line in f:
        go_id, index = line.split()
        # Subtract 1 to change to list indices.
        go_index_dct[int(index) - 1] = go_id
    f.close()

    go_go_matrix = []
    for i in range(len(go_index_dct.keys())):
        go_id_i = go_index_dct[i]
        current_row = []
        for j in range(len(go_index_dct.keys())):
            # GO nodes have no edge to themselves.
            if i == j:
                current_row += ['0']
                continue
            go_id_j = go_index_dct[j]
            if go_id_i in go_go_edge_dct:
                if go_id_j in go_go_edge_dct[go_id_i]:
                    current_row += ['1']
                else:
                    current_row += ['0']
            else:
                current_row += ['0']
        go_go_matrix += [current_row]


    f = open('./prediction_data/gene_gene_and_gene_go_weights.txt', 'r')
    out = open('./prediction_data/full_gene_go_matrix.txt', 'w')
    len_prev_line = 0
    go_index = 0
    for line in f:
        splitline = line.split()
        if len(splitline) >= len_prev_line:
            len_prev_line = len(splitline)
            out.write(line)
            continue
        assert(len_prev_line - len(splitline) == len(go_index_dct))
        out.write(line.strip() + '\t'.join(go_go_matrix[go_index]) + '\n')
        go_index += 1
    out.close()
    f.close()