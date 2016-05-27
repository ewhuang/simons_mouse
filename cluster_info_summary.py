### Author: Edward Huang

import file_operations
import sys

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiment.

def write_summary(run_num, mode):
    # Get cluster dictionary. Values are lists of genes.
    fname = './results/clusters_%s_%s.txt' % (mode, run_num)
    clst_go_dct = file_operations.get_cluster_dictionary(fname)

    # Get the number of gene-gene and gene-GO edges from the network.
    fname = './data/network_%s_%s.txt' % (mode, run_num)
    num_gg_net, num_ggo_net = file_operations.get_network_stats(fname)

    # Get density dictionary. Values are in-density and out-density.
    fname = './results/cluster_eval_%s_%s.txt' % (mode, run_num)
    density_dct = file_operations.get_cluster_densities(fname)

    # Find the best p-value GO enrichments for each cluster.
    fname = './results/cluster_enrichment_terms_%s_%s.txt' % (mode, run_num)
    enrichment_dct = file_operations.get_enrichment_dct(fname)

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
        for node in clus:
            if 'ENSMUSG' not in node:
                num_go += 1
        for edge in edge_list_go:
            node_a, node_b = edge
            if node_a in clus and node_b in clus:
                if 'ENSMUSG' not in node_a:
                    # This means we have a GO term.
                    num_ggo += 1
                elif 'ENSMUSG' not in node_b:
                    num_ggo += 1
                else:
                    # This means that we have a gene-gene edge.
                    num_gg += 1
        assert(num_gg % 2 == 0)
        assert(num_ggo % 2 == 0)
        # Divide each of these values by two because each edge is written
        # twice in the network.
        num_gg /= 2
        num_ggo /= 2
        in_dens, out_dens = density_dct[cid]
        ratio = in_dens / (in_dens + out_dens)
        out.write('%s\t%g\t%g\t%g\t' % (cid, in_dens, out_dens, ratio))
        out.write('%d\t%d\t%d\t%d\t' % (num_genes, num_go, num_gg, num_ggo))
        out.write('%s\n' % best_enrichment_dct[cid])
    out.close()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage:python %s run_num' % sys.argv[0]
        exit()
    run_num = sys.argv[1]

    for mode in ['go', 'no_go']:
        write_summary(run_num, mode)