### Author: Edward Huang

import file_operations
import sys
import time

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiment.

def write_summary(run_num, clus_fname, net_fname, eval_fname, enrich_fname,
    out_fname):
    # Get cluster dictionary. Values are lists of genes.
    clst_go_dct = file_operations.get_cluster_dictionary(clus_fname)

    # Get the number of gene-gene and gene-GO edges from the network.
    num_genes_net, num_gg_net, num_ggo_net = file_operations.get_network_stats(
        net_fname)

    # Get density dictionary. Values are in-density and out-density.
    density_dct = file_operations.get_cluster_densities(eval_fname)

    # Find the best p-value GO enrichments for each cluster.
    enrichment_dct = file_operations.get_enrichment_dct(enrich_fname)

    # Write out to file.
    out = open(out_fname, 'w')
    out.write('num_genes_in_net\tnum_g_g_net\tnum_g_go_net\n')
    out.write('%s\t%d\t%d\n' % (num_genes_net, num_gg_net, num_ggo_net))
    out.write('cluster_number\tin_dens\tout_dens\t')
    out.write('num_genes\tnum_go_terms_in\tnum_g_g_edges\tnum_g_go_edges\t')
    out.write('top_enrichment_p\n')
    for cid in density_dct:
        # cid = str(i + 1)
        clus = clst_go_dct[cid]
        num_genes = len(clus)
        num_go, num_gg, num_ggo = 0, 0, 0
        for node in clus:
            if 'ENSMUSG' not in node:
                num_go += 1
        # for edge in edge_list_go:
        #     node_a, node_b = edge
        #     if node_a in clus and node_b in clus:
        #         if 'ENSMUSG' not in node_a:
        #             # This means we have a GO term.
        #             num_ggo += 1
        #         elif 'ENSMUSG' not in node_b:
        #             num_ggo += 1
        #         else:
        #             # This means that we have a gene-gene edge.
        #             num_gg += 1
        # assert(num_gg % 2 == 0)
        # assert(num_ggo % 2 == 0)
        # # Divide each of these values by two because each edge is written
        # # twice in the network.
        # num_gg /= 2
        # num_ggo /= 2
        in_dens, out_dens, ratio = density_dct[cid]

        # ratio = in_dens / (in_dens + out_dens)
        out.write('%s\t%g\t%g\t%g\t' % (cid, in_dens, out_dens, ratio))
        out.write('%d\t%d\t%d\t%d\t' % (num_genes, num_go, num_gg, num_ggo))
        out.write('%s\n' % enrichment_dct[cid])
    out.close()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage:python %s run_num' % sys.argv[0]
        exit()
    run_num = sys.argv[1]

    start_time = time.time()

    # Get no_go summary.
    clus_fname = './results/clusters_no_go/clusters_no_go_%s.txt' % run_num
    net_fname = './data/networks_no_go/network_no_go_%s.txt' % run_num
    eval_fname = './results/cluster_eval_no_go/cluster_eval_no_go_%s.txt' % (
        run_num)
    enrich_fname = './results/cluster_enrichment_terms_no_go/'
    enrich_fname += 'cluster_enrichment_terms_no_go_%s.txt' % run_num
    out_fname = './results/clus_info_no_go/clus_info_no_go_%s.txt' % run_num
    write_summary(run_num, clus_fname, net_fname, eval_fname, enrich_fname,
        out_fname)

    # Get GO summaries.
    for domain_index in range(3):
        clus_fname = './results/clusters_go/clusters_go_%s_%d.txt' % (run_num,
            domain_index)
        net_fname = './data/networks_go/network_go_%s_%d.txt' % (run_num,
            domain_index)
        eval_fname = './results/cluster_eval_go/cluster_eval_go_%s_%d.txt' % (
            run_num, domain_index)
        enrich_fname = './results/cluster_enrichment_terms_go/'
        enrich_fname += 'cluster_enrichment_terms_go_%s_%d.txt' % (run_num,
            domain_index)
        out_fname = './results/clus_info_go/clus_info_go_%s_%d.txt' % (run_num,
            domain_index)
        write_summary(run_num, clus_fname, net_fname, eval_fname, enrich_fname,
            out_fname)

    print("--- %s seconds ---" % (time.time() - start_time))