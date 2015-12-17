### Author: Edward Huang

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiment.
### To perform on clusters without GO, simply change the name of the input
### and output files.

if __name__ == '__main__':
    print 'Extracting cluster data with GO...'
    # We use the "dirty" clusters with GO to analyze gene-GO edges.
    f = open('./results/clusters_go_1.txt', 'r')
    # Dictionary for clusters with GO. Keys are cluster indexes, values are the
    # list of genes in those clusters.
    clst_go_dct = {}
    for i, line in enumerate(f):
        if i == 0:
            continue
        newline = line.split()
        gene, cluster = newline[3], newline[5]
        if cluster in clst_go_dct:
            clst_go_dct[cluster] += [gene]
        else:
            clst_go_dct[cluster] = [gene]
    f.close()

    # We then read in the network to figure out how many GO terms went into the
    # clustering.
    print 'Extracting network data with GO...'
    num_genes_net = 0
    num_gg_net = 0
    num_ggo_net = 0
    f = open('./data/network_go_1.txt', 'r')
    edge_list_go = []
    for i, line in enumerate(f):
        print i
        if i == 0:
            continue
        if i == 1:
            num_genes_net = line.strip()
            continue
        # If not first two lines, find out if each edge is G-GO or G-G.
        node_a, node_b, weight = line.split()
        if 'ENSMUSG' not in node_a or 'ENSMUSG' not in node_b:
            num_ggo_net += 1
        else:
            num_gg_net += 1
        edge_list_go += [(node_a, node_b)]
    f.close()
    # Divide the two numbers by two to account for each edge in twice.
    assert(num_ggo_net % 2 == 0)
    assert(num_ggo_net % 2 == 0)
    num_ggo_net /= 2
    num_gg_net /= 2

    print 'Writing out information with GO...'
    # Write out to file.
    out = open('./results/clus_info_go_1.txt', 'w')
    out.write('num_genes_in_net\tnum_g_g_net\tnum_g_go_net\n')
    out.write('%s\t%d\t%d\n' % (num_genes_net, num_gg_net, num_ggo_net))
    out.write('num_genes\tnum_go_terms_in\tnum_g_g_edges\tnum_g_go_edges\n')
    for cid in clst_go_dct:
        clus = clst_go_dct[cid]
        num_genes = len(clus)
        num_go = 0
        num_gg = 0
        num_ggo = 0
        GO_terms = set([])
        for node in clus:
            if 'ENSMUSG' not in node:
                num_go += 1
        for edge in edge_list_go:
            node_a, node_b = edge
            if node_a in clus and node_b in clus:
                if 'ENSMUSG' not in node_a:
                    num_ggo += 1
                    GO_terms.add(node_a)
                elif 'ENSMUSG' not in node_b:
                    num_ggo += 1
                    GO_terms.add(node_b)
                else:
                    num_gg += 1
        assert(num_gg % 2 == 0)
        assert(num_ggo % 2 == 0)
        num_gg /= 2
        num_ggo /= 2
        out.write('%d\t%d\t%d\t%d\t' % (num_genes, num_go, num_gg, num_ggo))
        out.write('\t'.join(list(GO_terms)) + '\n')
    out.close()