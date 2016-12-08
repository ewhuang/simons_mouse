### Author: Edward Huang

import file_operations
import sys
import time

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiments.

def read_cheat_evaluation():
    '''
    Reads the cheating evaluation file. Returns a dictionary.
    Key: cluster ID -> str
    Value: p-value of t-test between in/(in + out) ratios of genes labeled
        by best GO term vs. ratios of genes not labeled -> float
    '''
    cheat_eval_dct = {}
    f = open('./results/%s_results/%s/cheating_evaluation/cheat_eval_%s.txt' % (
        data_type, objective_function, run_num), 'r')
    for line in f:
        # (cluster_id, p_value, labeled_mean, labeled_size, unlabeled_mean,
        #     unlabeled_size) = line.split()
        line = line.split()
        cheat_eval_dct[line[0]] = tuple(map(float, line[1:]))
    f.close()
    return cheat_eval_dct

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

    if 'prosnet' not in clus_fname and 'no_go' not in clus_fname:
        cheat_eval_dct = read_cheat_evaluation()

    # Write out to file.
    out = open(out_fname, 'w')
    # out.write('num_genes_in_net\tnum_g_g_net\tnum_g_go_net\n')
    # out.write('%s\t%d\t%d\n' % (num_genes_net, num_gg_net, num_ggo_net))
    
    out.write('Cluster ID\tIn-Density\tOut-Density\tNumber of genes\t'
        'Top GO enrichment p-value\tTop DBGAP enrichment p-value\t')
    # If it's a regular run with GO, add in the cheating evaluation.
    if 'prosnet' not in clus_fname and 'no_go' not in clus_fname:
        out.write('t-test p-value\tLabeled mean\tLabeled variance\t'
            'Labeled size\tUnlabeled mean\tUnlabeled variance\t'
            'Unlabeled size\t')
    if 'no_go' not in clus_fname:
        out.write('Number of cluster GO members\tCluster GO members')
    out.write('\n')

    for cid in sorted(density_dct.keys(), key=lambda x: int(x)):
        cluster_go_terms = []
        clus = clst_go_dct[cid]
        num_genes = len(clus)
        for node in clus:
            if 'ENSMUSG' not in node and 'ENSG' not in node:
                cluster_go_terms += [node]
        in_dens, out_dens = density_dct[cid]

        out.write('%s\t%g\t%g\t%d\t%s\t%s\t' % (cid, in_dens, out_dens,
            num_genes, go_enrichment_dct[cid], dbgap_enrichment_dct[cid]))

        if 'prosnet' not in clus_fname and 'no_go' not in clus_fname:
            out.write('%f\t%f\t%f\t%d\t%f\t%f\t%d\t' % cheat_eval_dct[cid])
        if 'no_go' not in clus_fname:
            out.write('%d\t%s' % (len(cluster_go_terms),
                '\t'.join(cluster_go_terms)))
        out.write('\n')
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
        net_fname = './data/%s_data/networks_%s/network_%s_%s.txt' % (data_type,
            network_type, network_type, run_num)
        # Perl script evaluation results filename.
        eval_fname = '%s/cluster_eval_%s/cluster_eval_%s_%s.txt' % format_str

        # GO enrichment results filename.
        go_enrich_fname = ('%s/cluster_enrichment_terms_%s/cluster_enrichment'
                            '_terms_%s_%s.txt') % format_str

        # DBGAP enrichment results filename.
        dbgap_enrich_fname = ('%s/dbgap_enrichment_terms_%s/dbgap_enrichment_'
                                'terms_%s_%s.txt') % format_str

        # Finally, the output filename.
        out_fname = '%s/clus_info_%s/clus_info_%s_%s.tsv' % format_str

        write_summary(clus_fname, net_fname, eval_fname, go_enrich_fname,
            dbgap_enrich_fname, out_fname)

def main():
    if len(sys.argv) != 4:
        print 'Usage:python %s data_type objective_function run_num' % (
            sys.argv[0])
        exit()
    global data_type, objective_function, run_num
    data_type, objective_function, run_num = sys.argv[1:]
    assert objective_function in ['oclode', 'schaeffer', 'wlogv']
    assert run_num.isdigit()

    # Integer data_type means it's a TCGA cancer.
    if 'prosnet' in data_type:
        data_type = data_type.split('_')[1]
        if data_type.isdigit():
            data_type = file_operations.get_tcga_disease_list()[int(data_type)]
        data_type = 'prosnet_' + data_type
    elif data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    generate_filenames()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))