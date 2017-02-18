### Author: Edward Huang

import file_operations
import json
import os
import sys
import time

### Creates files for the network with the GO terms and the network without.

def get_mf_dct():
    with open('./data/%s_data/mf_dct.json' % gene_type, 'r') as fp:
        mf_dct = json.load(fp)
    fp.close()
    return mf_dct

def get_go_go_edges():
    with open('./data/mf_go_go_dct.json', 'r') as fp:
        mf_go_go_dct = json.load(fp)
    fp.close()
    return mf_go_go_dct

def write_orth_file(gene_example):
    '''
    Writes out the orthology file.
    '''
    # Write out the orth file.
    orth_out = open('./data/%s_data/orth.txt' % gene_type, 'w')
    orth_out.write('0\t0\t%s\t%s' % (gene_example, gene_example))
    orth_out.close()

def write_clustering_inputs(edge_dct, network_genes, net_type):
    data_folder = './data/%s_data' % gene_type

    num_nodes = len(network_genes)
    if net_type == 'go':
        mf_dct, mf_go_go_dct = get_mf_dct(), get_go_go_edges()
        num_nodes += len(mf_dct)

    # Write out network genes to DCA file.
    dca_genes_out = open('%s/dca_genes_%s.txt' % (data_folder, net_type), 'w')
    for gene in network_genes:
        dca_genes_out.write('%s\tA\n' % gene)

    net_out = open('%s/network_%s.txt' % (data_folder, net_type), 'w')
    net_out.write('0\n%d\n' % num_nodes)
    real_out = open('%s/real_network_%s.txt' % (data_folder, net_type), 'w')
    real_out.write('Real network\n')
    dca_edges_out = open('%s/dca_edges_%s.txt' % (data_folder, net_type), 'w')
    # Write out GO-GO edges if this is the with-GO method. Only for DCA network.
    if net_type == 'go':
        for go in mf_go_go_dct:
            if go not in mf_dct:
                continue
            go_neighbor_list = mf_go_go_dct[go]
            for go_neighbor in go_neighbor_list:
                if go_neighbor not in mf_dct:
                    continue
                dca_edges_out.write('%s\t%s\t1\tg\n' % (go, go_neighbor))
        # Write gene-GO edges.
        for go in mf_dct:
            dca_genes_out.write('%s\tB\n' % go)

            for gene in mf_dct[go]:
                if gene not in network_genes:
                    continue
                # Write in each edge twice to make it undirected.
                net_out.write('%s\t%s\t1\n%s\t%s\t1\n' % (gene, go, go, gene))
                real_out.write('0\t%s\t%s\t1\n0\t%s\t%s\t1\n' % (gene, go, go,
                    gene))
                dca_edges_out.write('%s\t%s\t1\to\n' % (gene, go))

    # Write out co-expression edges.
    for gene_a, gene_b in edge_dct:
        net_out.write('%s\t%s\t1\n%s\t%s\t1\n' % (gene_a, gene_b, gene_b,
            gene_a))
        real_out.write('0\t%s\t%s\t1\n0\t%s\t%s\t1\n' % (gene_a, gene_b, gene_b,
            gene_a))
        dca_edges_out.write('%s\t%s\t1\tn\n' % (gene_a, gene_b))

    dca_genes_out.close()
    dca_edges_out.close()
    real_out.close()
    net_out.close()

    write_orth_file(gene_a)

def main():
    if len(sys.argv) not in [2, 3]:
        print ('Usage:python %s mouse/tcga/tcga_idx prosnet<optional>' %
            sys.argv[0])
        exit()
    global gene_type
    gene_type = sys.argv[1]
    assert gene_type in ['mouse', 'tcga'] or gene_type.isdigit()

    if gene_type.isdigit():
        gene_type = file_operations.get_tcga_list()[int(gene_type)]

    if len(sys.argv) == 3:
        assert sys.argv[2] == 'prosnet'
        gene_type = 'prosnet_' + gene_type

    edge_dct = file_operations.get_edge_dct(gene_type)
    network_genes = file_operations.get_network_genes(gene_type)

    for net_type in ['go', 'none']:
        write_clustering_inputs(edge_dct, network_genes, net_type)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))