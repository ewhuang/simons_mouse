### Author: Edward Huang

import file_operations
import sys
import time

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiment.
### Run time: 1 minute.

def get_cluster_dictionary():
    '''
    Returns a dictionary.
    Key: cluster ID -> str
    Value: list of genes in the corresponding clusters -> list(str)
    '''
    cluster_wgcna_dct = {}
    f = open('./results/%s_results/clusters.txt' % data_type, 'r')
    # Read in the cluster file to create the cluster dictionary.
    for i, line in enumerate(f):
        if i == 0:
            continue
        newline = line.strip().split('\t')
        cluster = newline[2][len('Cluster '):]
        # Skip trashcan clusters.
        if cluster == '0':
            continue
        gene = newline[1][len('Gene '):]
        if cluster not in cluster_wgcna_dct:
            cluster_wgcna_dct[cluster] = []
        cluster_wgcna_dct[cluster] += [gene]
    f.close()
    return cluster_wgcna_dct

def get_density_dct():
    '''
    Get the output files of the Perl script evaluate_clustering.pl and find
    the in-density and out-density of each cluster.
    Key: cluster ID -> str
    Value: (in-density, out-density) -> (float, float)
    '''
    density_dct = {}
    f = open('./results/%s_results/cluster_eval.txt' % data_type, 'r')
    for i, line in enumerate(f):
        if line[:7] != 'Cluster':
            continue
        line = line.split()
        clus_id, in_density, out_density = line[1], line[7], line[9]
        density_dct[clus_id] = (float(in_density), float(out_density))
    f.close()
    return density_dct

def get_enrichment_dct(label_type):
    '''
    Find the best p-value GO enrichments for each cluster.
    '''
    enrichment_dct = {}
    f = open('./results/%s_results/%s_enrichment_terms.txt' % (data_type,
        label_type), 'r')
    while True:
        line = f.readline()
        if line == '':
            break
        line = line.split()
        if line[0] == 'Cluster':
            cid = line[1]
            # Skip two lines, and read in the top p-value.
            line = f.readline()
            term_list = line.split()
            line = f.readline()
            p_list = line.split()
            enrichment_dct[cid] = (term_list, p_list)
    f.close()
    return enrichment_dct

def write_summary(out_fname, density_dct, cluster_wgcna_dct, go_enrichment_dct,
    dbgap_enrichment_dct, gwas_enrichment_dct):
    # Write out to file.
    out = open(out_fname, 'w')
    out.write('Cluster ID\tIn-Density\tOut-Density\tIn/(In+Out)\t'
        'Number of genes\t1st GO p-value\t1st GO term\t2nd GO p-value\t'
        '2nd GO term\t3rd GO p-value\t3rd GO term\t4th GO p-value\t4th GO term'
        '\t5th GO p-value\t5th GO term\tTop DBGAP p-value\tTop DBGAP term\t'
        'Top GWAS p-value\tTop GWAS term\n')
    # Loop through the clusters of genes.
    for cid in sorted(cluster_wgcna_dct.keys(), key=lambda x: int(x)):
        clus = cluster_wgcna_dct[cid]
        in_dens, out_dens = density_dct[cid]

        interleaved_terms_and_p = '\t'.join([val for pair in zip(
            go_enrichment_dct[cid][1], go_enrichment_dct[cid][0]
            ) for val in pair])

        best_dbgap_p = dbgap_enrichment_dct[cid][1][0]
        best_dbgap_term = dbgap_enrichment_dct[cid][0][0]

        best_gwas_p = gwas_enrichment_dct[cid][1][0]
        best_gwas_term = gwas_enrichment_dct[cid][0][0]

        out.write('%s\t%g\t%g\t%g\t%d\t%s\t%s\t%s\t%s\t%s\n' % (cid, in_dens, out_dens,
            in_dens / (in_dens + out_dens), len(clus), interleaved_terms_and_p,
            best_dbgap_p, best_dbgap_term, best_gwas_p, best_gwas_term))
    out.close()

def main():
    if len(sys.argv) != 2:
        print 'Usage: %s data_type' % sys.argv[0]
        exit()
    global data_type
    data_type = sys.argv[1]
    assert data_type == 'mouse' or data_type.isdigit()

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    go_enrichment_dct = get_enrichment_dct('go')
    dbgap_enrichment_dct = get_enrichment_dct('dbgap')
    gwas_enrichment_dct = get_enrichment_dct('gwas')
    
    cluster_wgcna_dct = get_cluster_dictionary()
    density_dct = get_density_dct()

    write_summary('./results/%s_results/clus_info.tsv' % data_type, density_dct,
        cluster_wgcna_dct, go_enrichment_dct, dbgap_enrichment_dct,
        gwas_enrichment_dct)
    write_summary('../results/plots_and_tables/%s_wgcna.tsv' % data_type,
        density_dct, cluster_wgcna_dct, go_enrichment_dct,
        dbgap_enrichment_dct, gwas_enrichment_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))