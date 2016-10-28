### Author: Edward Huang

import file_operations
import sys
import time

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiments.

def write_summary(clus_fname, net_fname, eval_fname, go_enrich_fname, 
    dbgap_enrich_fname, out_fname):
    '''
    Reads each of the files, and creates a summary table suitable for reading.
    '''
    # Get cluster dictionary. Values are lists of genes.
    clst_go_dct = file_operations.get_cluster_dictionary(clus_fname)
    # Get the number of gene-gene and gene-GO edges from the network.
    num_genes_net, num_gg_net, num_ggo_net = file_operations.get_network_stats(
        net_fname)
    # Get density dictionary. Values are in-density and out-density.
    density_dct = file_operations.get_cluster_densities(eval_fname)
    # Find the best p-value GO enrichments for each cluster.
    go_enrichment_dct = file_operations.get_enrichment_dct(go_enrich_fname)
    dbgap_enrichment_dct = file_operations.get_enrichment_dct(
        dbgap_enrich_fname)

    # Write out to file.
    out = open(out_fname, 'w')
    out.write('num_genes_in_net\tnum_g_g_net\tnum_g_go_net\n')
    out.write('%s\t%d\t%d\n' % (num_genes_net, num_gg_net, num_ggo_net))
    out.write('cluster_number\tin_dens\tout_dens\tratio\tnum_genes\t'
                'num_go_terms_in\tnum_g_g_edges\tnum_g_go_edges\t'
                'top_go_enrich_p\ttop_dbgap_enrich_p\tgo_terms_in\n')
    for cid in density_dct:
        cluster_go_terms = []
        clus = clst_go_dct[cid]
        num_genes = len(clus)
        num_go, num_gg, num_ggo = 0, 0, 0
        for node in clus:
            if ('ENSMUSG' not in node and 'ENSG' not in node):
                num_go += 1
                cluster_go_terms += [node]
        in_dens, out_dens, ratio = density_dct[cid]

        out.write('%s\t%g\t%g\t%g\t' % (cid, in_dens, out_dens, ratio)
                    '%d\t%d\t%d\t%d\t' % (num_genes, num_go, num_gg, num_ggo)
                    '%s\t%s\t%s\n' % (go_enrichment_dct[cid], 
                        dbgap_enrichment_dct[cid], '\t'.join(cluster_go_terms)))
    out.close()

def generate_filenames():
    '''
    Generates the filenames for each of the files we read from to summarize 
    the clustering performances of a particular run.
    '''
    results_folder = './results/%s_results/%s' % (data_type, objective_function)

    for network_type in ['go', 'no_go']:
        # Generate the format string.
        format_str = (results_folder, network_type, network_type, run_num)
        # Simulated annealing cluster results filename.
        clus_fname = '%s/clusters_%s/clusters_%s_%s.txt' % format_str
        # Input network filename.
        net_fname = './data/%s_data/networks_%s/network_%s_%s.txt' % format_str
        # Perl script evaluation results filename.
        eval_fname = '%s/cluster_eval_%s/cluster_eval_%s_%s.txt' % format_str

        # GO enrichment results filename.
        go_enrich_fname = ('%s/cluster_enrichment_terms_%s/cluster_enrichment'
                            '_terms_%s_%s.txt') % format_str

        # DBGAP enrichment results filename.
        dbgap_enrich_fname = ('%s/dbgap_enrichment_terms_%s/dbgap_enrichment_'
                                'terms_%s_%s.txt') % format_str

        # Finally, the output filename.
        out_fname = '%s/clus_info_%s/clus_info_%s_%s.txt' % format_str

        write_summary(clus_fname, net_fname, eval_fname, go_enrich_fname,
            dbgap_enrich_fname, out_fname)

def main():
    if len(sys.argv) != 4:
        print 'Usage:python %s data_type objective_function run_num' % (
            sys.argv[0])
        exit()
    global data_type, objective_function, run_num
    data_type, objective_function, run_num = sys.argv[1:]
    assert data_type in ['mouse', 'prosnet_mouse'] or data_type.isdigit()
    assert objective_function in ['oclode', 'schaeffer', 'wlogv']
    assert run_num.isdigit()

    # Integer data_type means it's a TCGA cancer.
    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    generate_filenames()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))