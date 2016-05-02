### Author: Edward Huang

import sys
import time

### Generates readable, tab-separated file to provide information on the
### clusters generated from the simulated annealing experiment.
### Run time: 1 minute.

def get_cluster_dictionary(go_domain):
    '''
    Returns a dictionary, keys=cluster ID's, values=lists of genes in the
    corresponding clusters.
    '''
    cluster_wgcna_dct = {}
    f = open('./results/clusters_%s.txt' % go_domain, 'r')
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

def get_density_dct(go_domain):
    '''
    Get the output files of the Perl script evaluate_clustering.pl and find
    the in-density and out-density of each cluster.
    Keys are cluster ID's, and values are tuples.
    Values: (in-density, out-density, number of in-cluster edges)
    '''
    density_dct = {}
    f = open('./results/cluster_eval_%s.txt' % go_domain, 'r')
    for i, line in enumerate(f):
        if line[:7] != 'Cluster':
            continue
        line = line.split()
        clus_id, in_density, out_density = line[1], line[7], line[9]
        num_gg_edges = float(line[25])
        density_dct[clus_id] = (float(in_density), float(out_density),
            num_gg_edges)
    f.close()
    return density_dct

def get_network_statistics(run_num):
    '''
    Reads a particular network we generated, and returns the number of genes in
    the network and the number of gene-gene edges.
    '''
    num_genes_net = 0
    num_gg_net = 0
    f = open('../data/networks_no_go/network_no_go_%d.txt' % run_num, 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue
        if i == 1:
            num_genes_net = line.strip()
            continue
        line = line.split()
        num_gg_net += 1
    f.close()
    # Divide the by two to account for each edge listed twice.
    num_gg_net /= 2
    return num_genes_net, num_gg_net

def get_enrichment_dct(go_domain):
    '''
    Find the best p-value GO enrichments for each cluster.
    '''
    enrichment_dct = {}
    f = open('./results/cluster_enrichment_terms_%s.txt' % go_domain, 'r')
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
    num_genes_net, num_gg_net = get_network_statistics(42)

    for go_domain in ['bp', 'cc', 'mf']:
        cluster_wgcna_dct = get_cluster_dictionary(go_domain)
        density_dct = get_density_dct(go_domain)
        enrichment_dct = get_enrichment_dct(go_domain)

        # Write out to file.
        out = open('./results/clus_info_%s.txt' % go_domain, 'w')
        out.write('num_genes_in_net\tnum_g_g_net\tnum_g_go_net\n')
        # Automatically write 0 for the last column, since there are no gene-GO
        # edges for WGCNA.
        out.write('%s\t%d\t0\n' % (num_genes_net, num_gg_net))
        out.write('cluster_number\tin_dens\tout_dens\tin/(in+out)\t')
        out.write('num_genes\tnum_go_terms_in\tnum_g_g_edges\tnum_g_go_edges\t')
        out.write('top_enrichment_p\n')

        # Loop through the clusters of genes.
        for i in range(len(cluster_wgcna_dct)):
            print float(i) / len(cluster_wgcna_dct) * 100
            cid = str(i + 1)
            clus = cluster_wgcna_dct[cid]

            num_genes = len(clus)

            in_dens, out_dens, num_gg_edges = density_dct[cid]
            # if in_dens == 0:
            #     ratio = float('inf')
            # else:
            ratio = in_dens / (in_dens + out_dens)
            out.write('%s\t%g\t%g\t%g\t' % (cid, in_dens, out_dens, ratio))
            out.write('%d\t0\t%d\t0\t' % (num_genes, num_gg_edges))
            out.write('%s\n' % enrichment_dct[cid])
        out.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))