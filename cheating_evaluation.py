### Author: Edward Huang

import file_operations
import json
import numpy as np
import os
from scipy import stats
import sys
import time

### This script evaluates the resulting clusters. For each node in each cluster,
### first find the in/(in + out) ratio. The nodes labeled by the best enriched
### GO term should have roughly the same ratio as those not labeled by the term.
### Run time: 10 minutes.

def generate_folders():
    res_dir = './results/%s_results/%s/cheating_evaluation' % (data_type,
        objective_function)
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)

def read_cluster_file():
    '''
    Generates the cluster filename, and then returns the dictionary of clusters.
    '''
    cluster_fname = ('./results/%s_results/%s/clusters_go/clusters_go_clean'
                        '_%s.txt' % (data_type, objective_function, run_num))
    cluster_dct = file_operations.get_cluster_dictionary(cluster_fname)
    return cluster_dct

def get_in_out_edge_dct(cluster_dct):
    '''
    For each gene, find out how many edges it has leading to other clusters,
    and how many edges lead to other genes in its own cluster.
    Key: ENSMUSG gene ID -> str
    Value: [number of edges within cluster, number of edges leaving] ->
        list(int)
    '''
    def increment_edges(gene_density_pair_dct, gene_a, gene_b, in_out_index):
        '''
        Increments a gene pair's in-cluster or out-cluster edge count.
        in_out_index == 0 means in-cluster count. 1 means out-cluster count.
        '''
        if gene_a not in gene_density_pair_dct:
            gene_density_pair_dct[gene_a] = [0, 0]
        gene_density_pair_dct[gene_a][in_out_index] += 1
        if gene_b not in gene_density_pair_dct:
            gene_density_pair_dct[gene_b] = [0, 0]
        gene_density_pair_dct[gene_b][in_out_index] += 1

    gene_density_pair_dct = {}
    # Read through the network file.
    f = open('./data/%s_data/networks_no_go/network_no_go_%s.txt' % (data_type,
        run_num))
    for i, line in enumerate(f):
        # Skip header and every other line.
        if i < 2 or i % 2 == 1:
            continue
        gene_a, gene_b, weight = line.split()
        # Loop through the cluster dictionary, and determine the edge type.
        # for cluster in cluster_dct:
        #     cluster_genes = cluster_dct[cluster]
        for cluster_genes in cluster_dct.values():
            if gene_a in cluster_genes:
                if gene_b in cluster_genes:
                    # This is an in-cluster edge.
                    increment_edges(gene_density_pair_dct, gene_a, gene_b, 0)
                else:
                    # This is an out-cluster edge.
                    increment_edges(gene_density_pair_dct, gene_a, gene_b, 1)

            elif gene_b in cluster_genes:
                # This is also an out-cluster edge.
                increment_edges(gene_density_pair_dct, gene_a, gene_b, 1)
    f.close()
    return gene_density_pair_dct

def get_bp_dct():
    '''
    Fetch the BP dictionary, and remove GO terms that overlap too much with MF.
    '''
    with open('./data/%s_data/bp_dct.json' % data_type, 'r') as fp:
        bp_dct = json.load(fp)
    fp.close()
    return bp_dct

def get_best_go_per_cluster():
    '''
    Reads the GO enrichment file and determines what the best GO term for each
    cluster is. Returns a dictionary.
    Key: Cluster ID number -> str
    Value: GO term best enriched in the cluster -> str
    '''
    best_go_dct = {}
    f = open('./results/%s_results/%s/go_enrichment_terms_go/go_enrichment_'
        'terms_go_%s.txt' % (data_type, objective_function, run_num))
    line = f.readline()
    while line != '':
        cluster_number = line.split()[1]
        line = f.readline()
        best_go = line.strip().split()[0]
        line = f.readline()
        line = f.readline()
        best_go_dct[cluster_number] = best_go
    f.close()
    return best_go_dct

def evaluate_clusters(cluster_dct):
    '''
    Goes through each cluster, and checks whether the best GO term-labeled 
    nodes have roughly the same in/(in + out) as those that aren't directly
    labeled.
    '''
    gene_density_pair_dct = get_in_out_edge_dct(cluster_dct)
    # Normalize the densities of each cluster.
    num_total_genes = float(len(gene_density_pair_dct))
    for cluster_gene_list in cluster_dct.values():
        num_cluster_genes = float(len(cluster_gene_list))
        out_denominator = num_total_genes - num_cluster_genes
        for gene in cluster_gene_list:
            # Update in-count to in-density. Subtract 1 for self edge.
            in_count, out_count = gene_density_pair_dct[gene]
            in_count /= (num_cluster_genes - 1)
            out_count /= out_denominator
            # Convert the pair to an in/(in + out) ratio.
            gene_density_pair_dct[gene] = in_count / (in_count + out_count)

    best_go_dct = get_best_go_per_cluster()
    print best_go_dct
    bp_dct = get_bp_dct()

    # Write out the results to file.
    out = open('./results/%s_results/%s/cheating_evaluation/cheat_eval_%s.txt'
        % (data_type, objective_function, run_num), 'w')
    for cluster in cluster_dct:
        labeled_densities, unlabeled_densities = [], []
        best_go_genes = bp_dct[best_go_dct[cluster]]
        # Determine which genes are labeled by the best enriched GO term.
        for gene in cluster_dct[cluster]:
            if gene in best_go_genes:
                labeled_densities += [gene_density_pair_dct[gene]]
            else:
                unlabeled_densities += [gene_density_pair_dct[gene]]
        # Compare the two density lists.
        t_test = stats.ttest_ind(labeled_densities, unlabeled_densities)
        out.write('%s\t%g\t%f\t%f\t%d\t%f\t%f\t%d\n' % (cluster, t_test[1],
            np.mean(labeled_densities), np.var(labeled_densities),
            len(labeled_densities), np.mean(unlabeled_densities),
            np.var(unlabeled_densities), len(unlabeled_densities)))
    out.close()

def main():
    if len(sys.argv) != 4:
        print ('Usage:python %s data_type objective_function run_num' % 
            sys.argv[0])
        exit()
    global data_type, objective_function, run_num
    data_type, objective_function, run_num = sys.argv[1:]
    assert objective_function in ['oclode', 'schaeffer', 'wlogv', 'prosnet']
    assert run_num.isdigit()
    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    generate_folders()
    # Evaluate only on clusters with GO terms.
    cluster_dct = read_cluster_file()
    evaluate_clusters(cluster_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))