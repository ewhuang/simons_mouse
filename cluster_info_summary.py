### Author: Edward Huang

import file_operations
import sys

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiment.

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage:python %s run_num' % sys.argv[0]
        exit()
    run_num = sys.argv[1]

    config_dct = file_operations.read_config_file()[run_num]
    num_clusters = int(config_dct['num_clusters'])
    # Find the GO terms to GO names.
    go_index_dct = file_operations.get_go_index_dct()

    for mode in ['go', 'no_go']:
        # We use the "dirty" clusters with GO to analyze gene-GO edges.
        cluster_fname = './results/clusters_%s_%s.txt' % (mode, run_num)
        clst_go_dct = file_operations.create_cluster_dct(cluster_fname)
        assert len(clst_go_dct) == num_clusters

        network_fname = './data/network_%s_%s.txt' % (mode, run_num)
        (num_genes_net, num_gg_net, num_ggo_net,
            edge_list_go) = file_operations.get_network_stats(network_fname)

        eval_fname = './results/cluster_eval_%s_%s.txt' % (mode, run_num)
        dens_dct = file_operations.get_cluster_evaluation_densities(eval_fname)

        # Find the best p-value GO enrichments for each cluster.
        enrichment_fname = './results/cluster_enrichment_terms_%s_%s.txt' % (
            mode, run_num)
        best_enrichment_dct = file_operations.get_best_enrichment_dct(
            enrichment_fname)

        # Write out to file.
        out = open('./results/clus_info_%s_%s.txt' % (mode, run_num), 'w')
        out.write('num_genes_in_net\tnum_g_g_net\tnum_g_go_net\n')
        out.write('%s\t%d\t%d\n' % (num_genes_net, num_gg_net, num_ggo_net))
        out.write('cluster_number\tin_dens\tout_dens\tin/(in+out)\t')
        out.write('num_genes\tnum_go_terms_in\tnum_g_g_edges\tnum_g_go_edges\t')
        out.write('top_enrichment_p\n')
        for i in range(len(clst_go_dct)):
            cid = str(i + 1)
            clus = clst_go_dct[cid]
            num_genes = len(clus)
            num_go, num_gg, num_ggo = 0, 0, 0
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
            ratio = in_dens / (in_dens + out_dens)
            out.write('%s\t%g\t%g\t%g\t' % (cid, in_dens, out_dens, ratio))
            out.write('%d\t%d\t%d\t%d\t' % (num_genes, num_go, num_gg, num_ggo))
            out.write('%s\t' % best_enrichment_dct[cid])
            out.write('\t'.join(list(GO_terms)) + '\n')
        out.close()