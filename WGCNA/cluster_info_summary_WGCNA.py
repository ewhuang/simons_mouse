### Author: Edward Huang

import sys

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiment.

if __name__ == '__main__':
    # Find the GO terms to GO names.
    # Retrieved data from http://geneontology.org/page/download-annotations
    go_id_to_name_dct = {}
    f = open('../go_edge_prediction/prediction_data/go_to_name.txt', 'r')
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
    f = open('../go_edge_prediction/prediction_data/noisogoHash.txt', 'r')
    for line in f:
        go_id, index = line.split()
        # Subtract 1 to change to list indices.
        go_index_dct[int(index) - 1] = go_id_to_name_dct[go_id]
    f.close()

    print 'Extracting cluster data...'
    # We use the "dirty" clusters with GO to analyze gene-GO edges.
    f = open('../results/WGCNA_results/WGCNA_clusters_all_genes.txt', 'r')
    # Dictionary for clusters with GO. Keys are cluster indexes, values are
    # the list of genes in those clusters.
    clst_go_dct = {}
    for i, line in enumerate(f):
        if i == 0:
            continue
        newline = line.strip().split('\t')
        gene = newline[1][len('Gene '):]
        cluster = newline[2][len('Cluster '):]
        if cluster == '0':
            continue
        if cluster in clst_go_dct:
            clst_go_dct[cluster] += [gene]
        else:
            clst_go_dct[cluster] = [gene]
    f.close()

    # We then read in the network to figure out how many GO terms went into
    # the clustering.
    print 'Extracting network data...'
    num_genes_net = 0
    num_gg_net = 0
    num_ggo_net = 0
    f = open('../data/network_no_go_20.txt', 'r')
    edge_list_go = []
    for i, line in enumerate(f):
        if i == 0:
            continue
        if i == 1:
            num_genes_net = line.strip()
            continue
        # If not first two lines, find out if each edge is G-GO or G-G.
        node_a, node_b, weight = line.strip().split('\t')
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

    # Get the output files of the Perl script evaluate_clustering.pl and
    # find the in-density and out-density of each cluster.
    dens_dct = {}
    f = open('../results/WGCNA_results/WGCNA_cluster_eval.txt', 'r')
    for i, line in enumerate(f):
        if line[:7] != 'Cluster':
            continue
        line = line.split()
        clus_id, in_density, out_density = line[1], line[7], line[9]
        dens_dct[clus_id] = (float(in_density), float(out_density))
    f.close()

    # Find the best p-value GO enrichments for each cluster.
    best_enrichment_dct = {}
    f = open('../results/WGCNA_results/cluster_enrichment_terms_wgcna.txt', 'r')
    while True:
        line = f.readline()
        if line == '':
            break
        line = line.split()
        if line[0] == 'Cluster':
            cid = line[1]
            # Skip two lines, and read in the top p-value.
            line = f.readline()
            line = f.readline().split()
            best_enrichment_dct[cid] = line[0]
    f.close()

    print 'Writing out clustering summary...'
    # Write out to file.
    out = open('../results/WGCNA_results/clus_info_WGCNA.txt', 'w')
    out.write('num_genes_in_net\tnum_g_g_net\tnum_g_go_net\n')
    out.write('%s\t%d\t%d\n' % (num_genes_net, num_gg_net, num_ggo_net))
    out.write('cluster_number\tin_dens\tout_dens\tin/(in+out)\t')
    out.write('num_genes\tnum_go_terms_in\tnum_g_g_edges\tnum_g_go_edges\t')
    out.write('top_enrichment_p\n')
    for i in range(len(clst_go_dct)):
        cid = str(i + 1)
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
                    # This means we have a GO term.
                    num_ggo += 1
                    if node_a.isdigit():
                        # This means that we need to convert the GO index
                        # to the GO term name.
                        GO_terms.add(go_index_dct[int(node_a)])
                    else:
                        GO_terms.add(node_a)
                elif 'ENSMUSG' not in node_b:
                    num_ggo += 1
                    if node_b.isdigit():
                        GO_terms.add(go_index_dct[int(node_b)])
                    else:
                        GO_terms.add(node_b)
                else:
                    # This means that we have a gene-gene edge.
                    num_gg += 1
        assert(num_gg % 2 == 0)
        assert(num_ggo % 2 == 0)
        # Divide each of these values by two because each edge is written
        # twice in the network.
        num_gg /= 2
        num_ggo /= 2
        in_dens, out_dens = dens_dct[cid]
        if in_dens == 0:
            ratio = float('inf')
        else:
            ratio = in_dens / (in_dens + out_dens)
        out.write('%s\t%g\t%g\t%g\t' % (cid, in_dens, out_dens, ratio))
        out.write('%d\t%d\t%d\t%d\t' % (num_genes, num_go, num_gg, num_ggo))
        out.write('%s\t' % best_enrichment_dct[cid])
        out.write('\t'.join(list(GO_terms)) + '\n')
    out.close()