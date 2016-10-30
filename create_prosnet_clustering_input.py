### Author: Edward Huang

import file_operations
import numpy as np
import operator
from sklearn.metrics.pairwise import cosine_similarity
import sys
import time

### This script takes the low-dimensional vector representations created by
### prosnet and then computes pairwise cosine similarites between every pair of
### nodes. Then, write out to a file similar to those written by create_
### clustering_input.py.
### Run time: 40 seconds.

prosnet_folder = './Sheng/prosnet/result'

def get_gene_matrix(network_type):
    '''
    Reads the output of the Prosnet script. Makes a 2d list where each row 
    corresponds to a gene. List of genes returned as a second list.
    '''
    node_list, node_vec_rep_list = [], []

    f = open('%s/dca_genes_%s_%s.vec' % (prosnet_folder, network_type, run_num))
    for i, line in enumerate(f):
        # Skip the header line.
        if i == 0:
            continue
        line = line.split()
        node, node_vec_rep = line[0], line[1:]

        node_list += [node]
        # Convert the low-dimensional vector representation to floats.
        node_vec_rep = map(float, node_vec_rep)
        node_vec_rep_list += [node_vec_rep]

    f.close()
    assert len(node_list) == len(node_vec_rep_list)
    return node_list, np.array(node_vec_rep_list)

def get_edge_dct(network_type):
    '''
    Returns a dictionary of (gene, gene) -> cosine similarity. Also returns the
    list of unique nodes.
    '''
    node_list, node_vec_rep_list = get_gene_matrix(network_type)

    # Compute the pairwise cosine similarity values.
    cos_sim_matrix = np.abs(cosine_similarity(node_vec_rep_list))

    edge_dct, edge_genes = {}, set([])
    num_nodes = len(cos_sim_matrix)
    # Loop through the pairwise cosine similarity matrix.
    for row_idx, row in enumerate(cos_sim_matrix):
        gene_a = node_list[row_idx]
        for col_idx in range(row_idx + 1, num_nodes):
            cos = row[col_idx]
            # TODO: maybe make this threshold higher.
            if cos < 0.5:
                continue
            gene_b = node_list[col_idx]
            edge_dct[(gene_a, gene_b)] = cos
            edge_genes.add(gene_a)
            edge_genes.add(gene_b)
    # Get the top million genes.
    sorted_edge_dct = sorted(edge_dct.items(), key=operator.itemgetter(1),
        reverse=True)[:int(1e6)]
    return sorted_edge_dct, edge_genes

def write_file(network_type):
    '''
    Writes the networks of edges.
    '''
    sorted_edge_dct, edge_genes = get_edge_dct(network_type)

    out_folder = './data/%s_data/networks_%s' % (data_type, network_type)

    # Regular network file for clustering.
    out = open('%s/network_%s_%s.txt' % (out_folder, network_type, run_num),
        'w')
    out.write('0\n%d\n' % len(edge_genes))

    # Real network file for cluster evaluation.
    real_out = open('%s/real_network_%s_%s.txt' % (out_folder, network_type,
        run_num), 'w')
    real_out.write('Real network\n')

    for (gene_a, gene_b), edge_weight in sorted_edge_dct:
        # Write in each edge twice to make it undirected.
        out.write('%s\t%s\t%g\n%s\t%s\t%g\n' % (gene_a, gene_b, edge_weight,
            gene_b, gene_a, edge_weight))
        real_out.write('0\t%s\t%s\t%g\n0\t%s\t%s\t%g\n' % (gene_a, gene_b,
            edge_weight, gene_b, gene_a, edge_weight))
    out.close()
    real_out.close()

    # Write out the orth file.
    orth_out = open('./data/%s_data/orth.txt' % data_type, 'w')
    sample_gene = edge_genes.pop()
    orth_out.write('0\t0\t%s\t%s' % (sample_gene, sample_gene))
    orth_out.close()

def main():
    if len(sys.argv) != 3:
        print 'Usage:python %s data_type run_num' % sys.argv[0]
        exit()
    global data_type, run_num

    data_type, run_num = sys.argv[1:]
    assert data_type == 'mouse' or data_type.isdigit()
    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    data_type = 'prosnet_' + data_type

    for network_type in ['go', 'no_go']:
        write_file(network_type)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))