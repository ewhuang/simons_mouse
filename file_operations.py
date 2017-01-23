### Author: Edward Huang

from collections import OrderedDict
import math
import os

# standard_deviation_hist.py
# gene_edge_weights.py
def get_gene_expression_dct(data_type):
    '''
    Returns dictionary where keys are genes, and values are gene expression
    vectors. data_type can be either 'mouse' or any on of the TCGA diseases.
    '''
    f = open('./data/%s_data/expr.tsv' % data_type, 'r')
    gene_expression_dct = OrderedDict({})
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.split()
        gene, exp_vals = line[0], map(float, line[1:])
        assert 'ENSMUSG' in gene or 'ENSG' in gene
        assert gene not in gene_expression_dct
        gene_expression_dct[gene] = exp_vals
    f.close()
    return gene_expression_dct

# convert_embedding_matrix_to_edge_file.py
def get_embedding_genes():
    '''
    Mouse.embedding.id contains a file of ENSMUSG ID's separated by newlines.
    Returns the list of genes.
    '''
    embedding_genes = []
    f = open('./Sheng/data/network/integrated_network/Mouse.embedding.id', 'r')
    for i, line in enumerate(f):
        ensmusg_id = line.strip()
        assert 'ENSMUSG' in ensmusg_id
        embedding_genes += [ensmusg_id]
    f.close()
    return embedding_genes

# # compute_go_enrichment.py
def read_go_overlap(data_type):
    '''
    Gets the BP and MF GO terms that are too similar to each other. We exclude
    them from training, but still use them to evaluate.
    '''
    overlap_list = []
    if 'mouse' in data_type:
        f = open('./data/mouse_data/overlapping_bp_mf_go_labels.txt', 'r')
    else:
        f = open('./data/tcga_data/overlapping_bp_mf_go_labels.txt', 'r')
    for line in f:
        bp_label, mf_label, p_value = line.strip().split('\t')
        overlap_list += [(bp_label, mf_label)]
    f.close()
    return overlap_list

# evaluate_clustering.py
def create_clean_go_file(data_type, objective_function, run_num):
    '''
    Takes a clustered file from simulated annealing and outputs a clean file
    without the GO nodes. Reformats the file for WGCNA outputs.
    '''
    results_folder = './results/%s_results/%s' % (data_type, objective_function)
    if not os.path.exists(results_folder):
        os.makedirs(results_folder)
    cluster_folder = '%s/clusters_go' % results_folder
    if not os.path.exists(cluster_folder):
        os.makedirs(cluster_folder)

    if objective_function != 'wgcna':
        fname = '%s/clusters_go_%s.txt' % (cluster_folder, run_num)
    else:
        fname = './data/wgcna_data/%s_module_membership.txt' % data_type

    f = open(fname, 'r')
    out = open('%s/clusters_go_clean_%s.txt' % (cluster_folder, run_num), 'w')
    for i, line in enumerate(f):
        if i == 0:
            out.write(line)
            continue
        # Skip GO terms in clusters.
        if 'ENSMUSG' not in line and 'ENSG' not in line:
            continue
        # Have to convert line from WGCNA output to simulated annealing output.
        if objective_function == 'wgcna':
            # Ignore color and membership columns.
            node, module, color, membership = line.split()
            # Skip garbage module.
            if module == '0':
                continue
            # Remove quotiation marks around the ENSMUSG ID.
            node = node.strip('"')
            line = 'Species 0\tGene %s\tCluster %s\n' % (node, module)
        out.write(line)
    out.close()
    f.close()

# compute_go_enrichment.py
# cluster_info_summary.py
def get_cluster_dictionary(filename):
    '''
    Returns a dictionary of clusters.
    Key: cluster ID -> str
    Value: lists of genes in the cluster-> list(str)
    '''
    cluster_dct = {}
    f = open(filename, 'r')
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
        if cluster not in cluster_dct:
            cluster_dct[cluster] = []
        cluster_dct[cluster] += [gene]
    f.close()
    return cluster_dct

# cluster_info_summary.txt
def get_cluster_densities(eval_fname):
    '''
    Get the output files of the Perl script evaluate_clustering.pl and find
    the in-density and out-density of each cluster.
    Key: Cluster Index
    Value: (in-density, out-density)
    '''
    density_dct = {}
    f = open(eval_fname, 'r')
    for i, line in enumerate(f):
        if line[:7] != 'Cluster':
            continue
        line = line.split()
        if len(line) != 35:
            continue
        clus_id, in_density, out_density = line[1], line[7], line[9]
        density_dct[clus_id] = (float(in_density), float(out_density))
    f.close()
    return density_dct

# cluster_info_summary.py
def get_enrichment_dct(enrichment_fname):
    '''
    Find the best p-value GO enrichments for each cluster.
    '''
    enrichment_dct = {}
    f = open(enrichment_fname, 'r')
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

# cluster_info_summary.py
def get_network_stats(network_fname):
    '''
    Reads a network and counts the number of genes, gene-gene edges, and gene-GO
    edges in the network.
    '''
    num_gg_net, num_ggo_net = 0, 0
    
    f = open(network_fname, 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue
        if i == 1:
            num_genes_net = line.strip()
            continue

        if 'GO:' in line:
            num_ggo_net += 1
        else:
            num_gg_net += 1
    f.close()
    # Divide the two numbers by two to account for each edge in twice.
    return num_genes_net, num_gg_net / 2, num_ggo_net / 2

# dump_go_dictionary_files.py
# gene_edge_weights.py
# compute_go_enrichment.py
def get_high_std_genes(data_type):
    '''
    Retrieves the list of genes that have high standard deviations across their
    gene expression vectors.
    '''
    high_std_genes = []
    f = open('./data/%s_data/high_std_genes.txt' % data_type, 'r')
    for line in f:
        gene = line.strip()
        assert 'ENSMUSG' in gene or 'ENSG' in gene
        high_std_genes += [gene]
    f.close()
    return high_std_genes

# create_clustering_input.py
def get_high_std_edge_dct(data_type):
    '''
    Returns edges between high standard deviation genes that have very
    correlated gene expression vectors.
    '''
    high_std_edge_dct = {}
    f = open('./data/%s_data/high_std_network.txt' % data_type, 'r')
    for i, line in enumerate(f):
        gene_a, gene_b, pcc = line.split()
        high_std_edge_dct[(gene_a, gene_b)] = pcc
    f.close()
    return high_std_edge_dct

# create_clustering_input.py
# simulated_annealing.py
def read_config_file(data_type):
    '''
    Returns dictionary. Keys are run_num strings, values are dicts of config
    options. Each dct has key of subgraph_decimal, temp, min_go_size,
    max_go_size, pearson/embedding edges, lambda, num_clusters.
    '''
    num_options = 6
    config_dct = {}
    
    # Get the number of clusters from the WGCNA runs.
    num_clusters = -1
    try:
        wgcna_f = open('./results/%s_results/wgcna/clus_info_go/clus_info_go_1'
            '.tsv' % data_type, 'r')
        for line in wgcna_f:
            num_clusters += 1
        wgcna_f.close()
    except:
        pass

    if 'prosnet' in data_type:
        data_type = '_'.join(data_type.split('_')[1:])
    if 'mouse' in data_type:
        f = open('mouse_config.txt', 'r')
    else:
        f = open('tcga_config.txt', 'r')
    for i, line in enumerate(f):
        config_num = i % (num_options + 1)
        if config_num == 0:
            run_num = line.split(' ')[1].strip()
            config_dct[run_num] = {}
        elif config_num == 1:
            temp = line.split()[2].strip()
            config_dct[run_num]['temp'] = temp
        elif config_num == 2:
            min_go_size = line.split()[2].strip()
            config_dct[run_num]['min_go_size'] = int(min_go_size)
        elif config_num == 3:
            max_go_size = line.split()[2].strip()
            config_dct[run_num]['max_go_size'] = int(max_go_size)
        elif config_num == 4:
            edge_method = line.strip()
            config_dct[run_num]['edge_method'] = edge_method
        elif config_num == 5:
            lamb = line.split()[2].strip()
            config_dct[run_num]['lamb'] = float(lamb)
        config_dct[run_num]['num_clusters'] = num_clusters
    f.close()
    return config_dct

# standard_deviation_hist.py
# full_pipeline.py
def get_tcga_disease_list():
    '''
    Get the list of valid diseases from the TCGA dataset.
    '''
    tcga_disease_list = []
    f = open('./data/tcga_data/tcga_diseases.txt', 'r')
    for line in f:
        tcga_disease_list += [line.strip()]
    f.close()
    return tcga_disease_list