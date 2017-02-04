### Author: Edward Huang

import file_operations
import os
import sys
import time

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiments.

subfolder = './results/plots_and_tables'
if not os.path.exists(subfolder):
    os.makedirs(subfolder)

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
    dbgap_enrich_fname, gwas_enrich_fname, kegg_enrich_fname, ctd_enrich_fname,
    out_fname):
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
    gwas_enrichment_dct = file_operations.get_enrichment_dct(gwas_enrich_fname)
    kegg_enrichment_dct = file_operations.get_enrichment_dct(kegg_enrich_fname)
    ctd_enrichment_dct = file_operations.get_enrichment_dct(ctd_enrich_fname)

    if ('prosnet' not in clus_fname and 'no_go' not in clus_fname and
        objective_function != 'wgcna'):
        cheat_eval_dct = read_cheat_evaluation()

    # Write out to file.
    out = open(out_fname, 'w')
    
    out.write('Cluster ID\tIn-Density\tOut-Density\tIn/(In+Out)\t'
        'Number of genes\t1st GO p-value\t1st GO term\t2nd GO p-value\t'
        '2nd GO term\t3rd GO p-value\t3rd GO term\t4th GO p-value\t4th GO term'
        '\t5th GO p-value\t5th GO term\tTop DBGAP p-value\tTop DBGAP term\t'
        'Top GWAS p-value\tTop GWAS term\tTop KEGG p-value\tTop KEGG term\t'
        'Top CTD p-value\tTop CTD term\t')
    # If it's a regular run with GO, add in the cheating evaluation.
    if ('prosnet' not in clus_fname and 'no_go' not in clus_fname and
        objective_function != 'wgcna'):
        out.write('t-test p-value\tLabeled mean\tLabeled variance\t'
            'Labeled size\tUnlabeled mean\tUnlabeled variance\t'
            'Unlabeled size\t')
    if 'no_go' not in clus_fname and objective_function != 'wgcna':
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

        interleaved_terms_and_p = '\t'.join([val for pair in zip(
            go_enrichment_dct[cid][1], go_enrichment_dct[cid][0]
            ) for val in pair])

        best_dbgap_p = dbgap_enrichment_dct[cid][1][0]
        best_dbgap_term = dbgap_enrichment_dct[cid][0][0]

        best_gwas_p = gwas_enrichment_dct[cid][1][0]
        best_gwas_term = gwas_enrichment_dct[cid][0][0]

        best_kegg_p = kegg_enrichment_dct[cid][1][0]
        best_kegg_term = kegg_enrichment_dct[cid][0][0]

        best_ctd_p = ctd_enrichment_dct[cid][1][0]
        best_ctd_term = ctd_enrichment_dct[cid][0][0]

        out.write('%s\t%g\t%g\t%g\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (
            cid, in_dens, out_dens, in_dens / (in_dens + out_dens), num_genes,
            interleaved_terms_and_p, best_dbgap_p, best_dbgap_term, best_gwas_p,
            best_gwas_term, best_kegg_p, best_kegg_term, best_ctd_p,
            best_ctd_term))

        if ('prosnet' not in clus_fname and 'no_go' not in clus_fname and
            objective_function != 'wgcna'):
            out.write('%f\t%f\t%f\t%d\t%f\t%f\t%d\t' % cheat_eval_dct[cid])
        if 'no_go' not in clus_fname and objective_function != 'wgcna':
            out.write('%d\t%s' % (len(cluster_go_terms), '\t'.join(
                cluster_go_terms)))
        out.write('\n')
    out.close()

def generate_filenames():
    '''
    Generates the filenames for each of the files we read from to summarize 
    the clustering performances of a particular run.
    '''
    results_folder = './results/%s_results/%s' % (data_type, objective_function)

    for network_type in ['go', 'no_go']:
        # WGCNA name convention follow 'go' conventions.
        if objective_function == 'wgcna' and network_type == 'no_go':
            continue
        # Generate the format string.
        format_str = (results_folder, network_type, network_type, run_num)
        # Simulated annealing cluster results filename.
        if objective_function == 'wgcna':
            clus_fname = '%s/clusters_%s/clusters_%s_clean_%s.txt' % format_str
        else:
            clus_fname = '%s/clusters_%s/clusters_%s_%s.txt' % format_str
        # Input network filename.
        net_fname = './data/%s_data/networks_%s/network_%s_%s.txt' % (data_type,
            network_type, network_type, run_num)
        # Perl script evaluation results filename.
        eval_fname = '%s/cluster_eval_%s/cluster_eval_%s_%s.txt' % format_str

        # GO enrichment results filename.
        go_enrich_fname = ('%s/go_enrichment_terms_%s/go_enrichment_terms_%s_'
                            '%s.txt') % format_str

        # DBGAP enrichment results filename.
        dbgap_enrich_fname = ('%s/dbgap_enrichment_terms_%s/dbgap_enrichment_'
                                'terms_%s_%s.txt') % format_str

        # GWAS enrichment results filename.
        gwas_enrich_fname = ('%s/gwas_enrichment_terms_%s/gwas_enrichment_'
                                'terms_%s_%s.txt') % format_str

        # KEGG enrichment results filename.
        kegg_enrich_fname = ('%s/kegg_enrichment_terms_%s/kegg_enrichment_'
                                'terms_%s_%s.txt') % format_str

        # CTD disease enrichment results filename.
        ctd_enrich_fname = ('%s/ctd_enrichment_terms_%s/ctd_enrichment_'
                                'terms_%s_%s.txt') % format_str

        # Finally, the output filename.
        info_folder = '%s/clus_info_%s' % (results_folder, network_type)
        if not os.path.exists(info_folder):
            os.makedirs(info_folder)
        out_fname = '%s/clus_info_%s/clus_info_%s_%s.tsv' % format_str

        write_summary(clus_fname, net_fname, eval_fname, go_enrich_fname,
            dbgap_enrich_fname, gwas_enrich_fname, kegg_enrich_fname,
            ctd_enrich_fname, out_fname)
        # plots_and_tables filename.
        if objective_function == 'wgcna':
            subfolder_fname = '%s/%s_wgcna.tsv' % (subfolder, data_type)
        else:
            subfolder_fname = '%s/%s_%s.tsv' % (subfolder, data_type,
                network_type)
        write_summary(clus_fname, net_fname, eval_fname, go_enrich_fname,
            dbgap_enrich_fname, gwas_enrich_fname, kegg_enrich_fname,
            ctd_enrich_fname, subfolder_fname)

def main():
    if len(sys.argv) != 4:
        print 'Usage:python %s data_type objective_function run_num' % (
            sys.argv[0])
        exit()
    global data_type, objective_function, run_num
    data_type, objective_function, run_num = sys.argv[1:]
    assert objective_function in ['oclode', 'schaeffer', 'wlogv', 'wgcna']
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