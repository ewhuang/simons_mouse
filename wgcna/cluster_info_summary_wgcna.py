### Author: Edward Huang

import file_operations
import sys
import time

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiment.
### Run time: 1 minute.

def get_cluster_dictionary(go_method, go_domain):
    '''
    Returns a dictionary, keys=cluster ID's, values=lists of genes in the
    corresponding clusters.
    '''
    cluster_wgcna_dct = {}
    f = open('./results/%s_results/%s/clusters_%s.txt' % (data_type, go_method,
        go_domain), 'r')
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
            cluster_wgcna_dct[cluster] = [gene]
        else:
            cluster_wgcna_dct[cluster] += [gene]
    f.close()
    return cluster_wgcna_dct

def get_density_dct(go_method, go_domain):
    '''
    Get the output files of the Perl script evaluate_clustering.pl and find
    the in-density and out-density of each cluster.
    Keys are cluster ID's, and values are tuples.
    Values: (in-density, out-density, number of in-cluster edges)
    '''
    density_dct = {}
    f = open('./results/%s_results/%s/cluster_eval_%s.txt' % (data_type, go_method,
        go_domain), 'r')
    for i, line in enumerate(f):
        if line[:7] != 'Cluster':
            continue
        line = line.split()
        clus_id, in_density, out_density = line[1], line[7], line[9]
        density_dct[clus_id] = (float(in_density), float(out_density))
    f.close()
    return density_dct

def get_enrichment_dct(go_method, go_domain):
    '''
    Find the best p-value GO enrichments for each cluster.
    '''
    enrichment_dct = {}
    f = open('./results/%s_results/%s/cluster_enrichment_terms_%s_%s.txt' % (
        data_type, go_method, go_method, go_domain), 'r')
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
            enrichment_dct[cid] = line[0]
    f.close()
    return enrichment_dct

def get_dbgap_enrichment_dct(go_method):
    '''
    Find the best p-value DBGAP enrichments for each cluster.
    '''
    enrichment_dct = {}
    f = open('./results/%s_results/%s/dbgap_enrichment_terms_%s.txt' % (
        data_type, go_method, go_method), 'r')
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
            enrichment_dct[cid] = line[0]
    f.close()
    return enrichment_dct

def main():
    if len(sys.argv) != 4:
        print ('Usage: %s data_type genes_only/pca/mean/median network_to_'
            'compare' % sys.argv[0])
        exit()
    global data_type
    data_type = sys.argv[1]
    assert data_type == 'mouse' or data_type.isdigit()
    go_method = sys.argv[2]
    assert go_method in ['genes_only', 'pca', 'mean', 'median']
    comparing_network = sys.argv[3]
    assert comparing_network.isdigit()

    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    domain_list = ['bp']

    dbgap_enrichment_dct = get_dbgap_enrichment_dct(go_method)
    
    for go_domain in domain_list:
        cluster_wgcna_dct = get_cluster_dictionary(go_method, go_method)
        density_dct = get_density_dct(go_method, go_method)

        enrichment_dct = get_enrichment_dct(go_method, go_domain)
        out = open('./results/%s_results/%s/clus_info_%s_%s.tsv' % (data_type, 
            go_method, go_method, go_domain), 'w')

        # Write out to file.
        out.write('Cluster ID\tIn-Density\tOut-Density\tNumber of genes\t'
            'Top GO enrichment p-value\tTop DBGAP enrichment p-value\n')

        # Loop through the clusters of genes.
        for cid in sorted(cluster_wgcna_dct.keys(), key=lambda x: int(x)):
            clus = cluster_wgcna_dct[cid]

            in_dens, out_dens = density_dct[cid]

            out.write('%s\t%g\t%g\t%d\t%s\t%s\n' % (cid, in_dens, out_dens,
                len(clus), enrichment_dct[cid], dbgap_enrichment_dct[cid]))
        out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))