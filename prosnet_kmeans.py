### Author: Edward Huang

import os
from sklearn.cluster import KMeans, SpectralClustering
import sys
import time

### This script takes the low-dimensional vector representations of genes and
### GO terms as calculated by PROSNET on the DCA networks, runs k-means on them,
### then prepares a cluster file similar to that produced by simulated
### annealing.
### Run time:

def get_num_wgcna_clusters():
    '''
    This function reads the WGCNA equivalent, and gets the same number of 
    clusters that WGCNA made.
    '''
    wgcna_folder = './wgcna/results'
    f = open('%s/%s_results/genes_only/clus_info_genes_only_bp.txt' % (
        wgcna_folder, data_type), 'r')
    num_lines = 0
    for line in f:
        line = line.split()
        if len(line) < 4 or not line[0].isdigit():
            continue
        num_lines += 1
    f.close()
    return num_lines

def get_prosnet_feature_matrix(network_type):
    '''
    Given a network type (go/no_go), returns the feature matrix, read for
    k-means.
    '''
    assert network_type in ['go', 'no_go']
    prosnet_folder = './Sheng/prosnet/result'
    f = open('%s/dca_genes_%s_%s.vec' % (prosnet_folder, network_type, run_num))

    node_list, node_vec_rep_list = [], []
    for i, line in enumerate(f):
        # Skip the header line.
        if i == 0:
            continue
        line = line.split()
        node, node_vec_rep = line[0], line[1:]

        node_list += [node]
        # Convert the low-dimensional vector representation to floats.
        node_vec_rep = map(float, node_vec_rep)
        node_vec_rep_list += [node_vec_rep]
    f.close()
    assert len(node_list) == len(node_vec_rep_list)
    return node_list, node_vec_rep_list

def run_k_means(node_vec_rep_list):
    '''
    Runs k-means on the 2D list of low-dimensional vector representations,
    given a network from either GO or no-GO.
    '''
    num_clusters = get_num_wgcna_clusters()
    # TODO: currently running spectral.
    clf = KMeans(n_clusters=num_clusters, init='k-means++')
    # clf = SpectralClustering(n_clusters=num_clusters,
    #     affinity='nearest_neighbors')
    return clf.fit_predict(node_vec_rep_list)
    # labels = clf.labels_
    # return labels

def write_cluster_file(network_type, node_list, node_labels):
    '''
    Given a clustering result, write it out to file, in the same format that
    simulated annealing outputs.
    '''
    assert network_type in ['go', 'no_go']

    # Create the folders if they don't exist.
    results_folder = './results/%s_results/prosnet/' % data_type
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
        os.makedirs(results_folder + 'clus_info_no_go/')
        os.makedirs(results_folder + 'clus_info_go/')
        os.makedirs(results_folder + 'cluster_enrichment_terms_go/')
        os.makedirs(results_folder + 'cluster_enrichment_terms_no_go/')
        os.makedirs(results_folder + 'cluster_eval_no_go/')
        os.makedirs(results_folder + 'clusters_go/')
        os.makedirs(results_folder + 'cluster_eval_go/')
        os.makedirs(results_folder + 'clusters_no_go/')

    if network_type == 'go':
        out = open('%sclusters_%s/clusters_%s_%s_0.txt' % (results_folder,
            network_type, network_type, run_num), 'w')
    else:
        out = open('%sclusters_%s/clusters_%s_%s.txt' % (results_folder,
            network_type, network_type, run_num), 'w')
    out.write('dummy_header\n')
    for node_idx, label in enumerate(node_labels):
        out.write('Species 0\tGene %s\tCluster %d\n' % (node_list[node_idx],
            label+1))
    out.close()

def main():
    if len(sys.argv) != 3:
        print 'Usage:python %s mouse/tcga_cancer_index run_num' % sys.argv[0]
        exit()
    global data_type, run_num
    data_type = sys.argv[1]
    assert data_type == 'mouse' or data_type.isdigit()
    run_num = sys.argv[2]
    assert run_num.isdigit()

    # Get the list of vectors and node lists for each network type.
    go_node_list, go_node_vec_rep_list = get_prosnet_feature_matrix('go')
    no_go_node_list, no_go_node_vec_rep_list = get_prosnet_feature_matrix(
        'no_go')

    go_node_labels = run_k_means(go_node_vec_rep_list)
    assert len(go_node_labels) == len(go_node_list)

    no_go_node_labels = run_k_means(no_go_node_vec_rep_list)
    assert len(no_go_node_labels) == len(no_go_node_list)

    write_cluster_file('go', go_node_list, go_node_labels)
    write_cluster_file('no_go', no_go_node_list, no_go_node_labels)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))