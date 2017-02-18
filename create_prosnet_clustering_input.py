### Author: Edward Huang

import file_operations
import numpy as np
import operator
import os
from sklearn.metrics.pairwise import cosine_similarity
import sys
import time

### This script takes the low-dimensional vector representations created by
### prosnet and then computes pairwise cosine similarites between every pair of
### nodes. Then, write out to a file similar to those written by create_
### clustering_input.py.
### Run time: 40 seconds.

def get_gene_matrix(network_type):
    '''
    Reads the output of the Prosnet script. Makes a 2d list where each row 
    corresponds to a gene. List of genes returned as a second list.
    '''
    node_list, node_vec_rep_list = [], []
    prosnet_folder = './Sheng/prosnet/result/%s' % data_type
    f = open('%s/dca_genes_%s_%s.vec' % (prosnet_folder, network_type, run_num))
    for i, line in enumerate(f):
        # Skip the header line.
        if i == 0 or ('ENSMUSG' not in line and 'ENSG' not in line):
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

def get_genes_from_edges(edge_set):
    gene_set = set([])
    for gene_a, gene_b in edge_set:
        gene_set.add(gene_a)
        gene_set.add(gene_b)
    return gene_set

def get_edge_dct(network_type, base_edge_set):
    '''
    Returns a dictionary of (gene, gene) -> cosine similarity. Also returns the
    list of unique nodes.
    '''
    node_list, node_vec_rep_list = get_gene_matrix(network_type)

    # Compute the pairwise cosine similarity values. TODO: absolute value?
    # cos_sim_matrix = np.abs(cosine_similarity(node_vec_rep_list))
    # cos_sim_matrix = cosine_similarity(node_vec_rep_list)
    cos_sim_matrix = np.corrcoef(node_vec_rep_list)

    edge_set, edge_genes = set([]), set([])
    edge_dct = {}
    num_nodes = len(cos_sim_matrix)
    # Loop through the pairwise cosine similarity matrix.
    for row_idx, row in enumerate(cos_sim_matrix):
        gene_a = node_list[row_idx]
        for col_idx in range(row_idx + 1, num_nodes):
            cos = row[col_idx]
            # TODO: maybe make this threshold higher.
            # if cos > 0.95:
            #     gene_b = node_list[col_idx]
            #     edge_set.add((gene_a, gene_b))
            #     edge_set.add((gene_b, gene_a))
            #     # edge_dct[(gene_a, gene_b)] = cos
            #     edge_genes.add(gene_a)
            #     edge_genes.add(gene_b)

            if cos < 0.3:
                continue
            gene_b = node_list[col_idx]
            if (gene_a, gene_b) in base_edge_set:
                continue
            edge_dct[(gene_a, gene_b)] = cos
    # Get the top million genes. TODO: setting the number of edges.
    sorted_edge_dct = sorted(edge_dct.items(), key=operator.itemgetter(1),
        reverse=True)[:int(200000)]
    edge_set = set([edge for edge, val in sorted_edge_dct])
    # return sorted_edge_dct, edge_genes
    edge_genes = get_genes_from_edges(edge_set)
    return edge_set, edge_genes

def write_file(network_type, base_edge_set, base_genes):
    '''
    Writes the networks of edges.
    '''
    # sorted_edge_dct, edge_genes = get_edge_dct(network_type)
    edge_set, edge_genes = get_edge_dct(network_type, base_edge_set)

    out_folder = './data/%s_data/networks_%s' % (data_type, network_type)
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)

    # Regular network file for clustering.
    out = open('%s/network_%s_%s.txt' % (out_folder, network_type, run_num),
        'w')
    out.write('0\n%d\n' % len(edge_genes.union(base_genes)))

    # Real network file for cluster evaluation.
    real_out = open('%s/real_network_%s_%s.txt' % (out_folder, network_type,
        run_num), 'w')
    real_out.write('Real network\n')

    # combined_edges = edge_set.union(base_edge_set)
    assert len(base_edge_set) % 2 == 0 # Should have an even number of edges.
    for (gene_a, gene_b) in base_edge_set:
        out.write('%s\t%s\t%d\n' % (gene_a, gene_b, 1))
        real_out.write('0\t%s\t%s\t%d\n' % (gene_a, gene_b, 1))

    # Have to write twice for edge set.
    for (gene_a, gene_b) in edge_set:
        out.write('%s\t%s\t%d\n%s\t%s\t%d\n' % (gene_a, gene_b, 1, gene_b,
            gene_a, 1))
        real_out.write('0\t%s\t%s\t%d\n0\t%s\t%s\t%d\n' % (gene_a, gene_b, 1,
            gene_b, gene_a, 1))
    out.close()
    real_out.close()

    # Write out the orth file.
    orth_out = open('./data/%s_data/orth.txt' % data_type, 'w')
    sample_gene = edge_genes.pop()
    orth_out.write('0\t0\t%s\t%s' % (sample_gene, sample_gene))
    orth_out.close()

def read_base_network(data_type):
    assert 'prosnet' not in data_type
    base_edge_set, base_genes = set([]), set([])
    f = open('./data/%s_data/networks_no_go/network_no_go_%s.txt' % (data_type,
        run_num), 'r')
    for i, line in enumerate(f):
        if i == 0: # Skip header line.
            continue
        if i == 1:
            num_genes = int(line.strip())
            continue
        gene_a, gene_b, weight = line.split()
        base_edge_set.add((gene_a, gene_b))
        base_genes.add(gene_a)
        base_genes.add(gene_b)
    f.close()
    assert len(base_genes) == num_genes
    return base_edge_set, base_genes

def main():
    if len(sys.argv) != 3:
        print 'Usage:python %s data_type run_num' % sys.argv[0]
        exit()
    global data_type, run_num
    data_type, run_num = sys.argv[1:]
    assert data_type == 'mouse' or data_type.isdigit()
    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    base_edge_set, base_genes = read_base_network(data_type)

    data_type = 'prosnet_' + data_type

    for network_type in ['go', 'no_go']:
        write_file(network_type, base_edge_set, base_genes)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))