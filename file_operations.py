### Author: Edward Huang

from collections import OrderedDict
import math

# standard_deviation_hist.py
# gene_edge_weights.py
def get_gene_expression_dct(data_type):
    '''
    Returns dictionary where keys are genes, and values are gene expression
    vectors.
    '''
    if data_type == 'mouse':
        f = open('./data/mm_mrsb_log2_expression.tsv', 'r')
    else:
        f = open('./data/tcga_expr.txt', 'r')
    gene_expression_dct = OrderedDict({})
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.split()
        gene, exp_vals = line[0], line[1:]
        assert 'ENSMUSG' in gene or 'ENSG' in gene
        exp_vals = [float(val) for val in exp_vals]
        assert gene not in gene_expression_dct
        gene_expression_dct[gene] = exp_vals
    f.close()
    return gene_expression_dct

# standard_deviation_hist.py
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

# evaluate_clustering.py
def create_clean_go_file(data_type, objective_function, run_num, go_domain_num):
    '''
    Takes a clustered file from simulated annealing and outputs a clean file
    without the GO nodes.
    '''
    f = open('./%s_results/%s/clusters_go/clusters_go_%s_%d.txt' % (
        data_type, objective_function, run_num, go_domain_num), 'r')
    out = open('./%s_results/%s/clusters_go/clusters_go_clean_%s_%d.txt' % (
        data_type, objective_function, run_num, go_domain_num), 'w')
    for i, line in enumerate(f):
        if i == 0:
            out.write(line)
            continue
        # Skip GO terms in clusters.
        if 'ENSMUSG' not in line or 'ENSG' not in line:
            continue
        out.write(line)
    out.close()
    f.close()

# compute_go_enrichment.py
def get_cluster_dictionary(filename):
    '''
    Returns a dictionary, keys=cluster ID's, values=lists of genes in the
    corresponding clusters.
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
            cluster_dct[cluster] = [gene]
        else:
            cluster_dct[cluster] += [gene]
    f.close()
    return cluster_dct

# make_gene_go_file.py
def get_coexpression_genes():
    '''
    Reads the original co-expression matrix, and returns the set of genes that
    appear in the matrix.
    '''
    coexpression_genes = set([])
    f = open('./data/mm_mrsb_log2_expression.tsv', 'r')
    gene_lst = []
    for i, line in enumerate(f):
        if i == 0:
            continue
        gene = line.split()[0]
        assert 'ENSMUSG' in gene
        coexpression_genes.add(gene)
    f.close()
    return coexpression_genes

# make_gene_go_file.py
def get_gene_index_dct():
    '''
    Returns a dictionary that maps gene indices to their ENSMUSG id's.
    '''
    # We only want coexpression genes in gene-GO relationships.
    coexpression_genes = get_coexpression_genes()

    mgi_to_ensembl_dct = {}
    f = open('./go_edge_prediction/prediction_data/mgi_to_ensembl.txt', 'r')
    for i, line in enumerate(f):
        # Skip the header line.
        if i == 0:
            continue
        line = line.split()
        if len(line) != 5 or line[1] != 'current':
            continue
        ensmusg_id = line[4]
        assert 'ENSMUSG' in ensmusg_id
        # Skip if a gene isn't in our coexpression network, or if it isn't up
        # to date.
        if ensmusg_id not in coexpression_genes:
            continue
        mgi_id = line[0]
        if mgi_id in mgi_to_ensembl_dct:
            mgi_to_ensembl_dct[mgi_id] += [ensmusg_id]
        else:
            mgi_to_ensembl_dct[mgi_id] = [ensmusg_id]
    f.close()

    # Keys are the indices in the edge weight matrix, values are the genes.
    gene_index_dct = {}
    f = open('./go_edge_prediction/prediction_data/noisogeneHash.txt', 'r')
    for line in f:
        mgi_id, row = line.split()
        if mgi_id not in mgi_to_ensembl_dct:
            continue
        # Subtract 1 to change to list indices.
        row_index = int(row) - 1
        gene_index_dct[row_index] = mgi_to_ensembl_dct[mgi_id]
    f.close()
    # Make sure we indeed subtract 1 from the row indices.
    assert -1 not in gene_index_dct
    return gene_index_dct

# create_clustering_input.py
def get_go_id_to_name_dct():
    '''
    Returns a dictionary that maps GO ID's to their English names.
    '''
    go_id_to_name_dct = {}
    f = open('./go_edge_prediction/prediction_data/go_to_name.txt', 'r')
    while True:
        line = f.readline()
        if line == '':
            break
        if line.strip() == '[Term]':
            go_id = f.readline().split()[1].lower()
            go_name = '_'.join(f.readline()[len('name: '):].split())

            assert go_id not in go_id_to_name_dct
            go_id_to_name_dct[go_id] = go_name
            f.readline()
            next = f.readline()
            while 'alt_id' in next:
                go_alt_id = next.split()[1].lower()
                assert go_alt_id not in go_id_to_name_dct
                go_id_to_name_dct[go_alt_id] = go_name
                next = f.readline()
    f.close()
    # Manual tests.
    assert go_id_to_name_dct['go:0000001'] == 'mitochondrion_inheritance'
    assert go_id_to_name_dct['go:0000011'] == 'vacuole_inheritance'
    return go_id_to_name_dct

# make_gene_go_file.py
def get_go_index_to_id_dct():
    '''
    This function returns a dictionary where keys are GO id's, and values are
    their corresponding indices.
    '''
    # Keys are the GO ids, values are the indices in the edge weight matrix.
    go_index_dct = {}
    f = open('./go_edge_prediction/prediction_data/noisogoHash.txt', 'r')
    for line in f:
        go_id, index = line.split()
        # Subtract 1 to change to list indices.
        go_index_dct[int(index) - 1] = go_id
    f.close()
    # Make sure we indeed subtract 1 from the indices.
    assert -1 not in go_index_dct
    return go_index_dct

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
        ratio = float(line[34])
        density_dct[clus_id] = (float(in_density), float(out_density), ratio)
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
            line = f.readline().split()
            enrichment_dct[cid] = line[0]
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

# # Keyword is either 'full' or 'sampled'. Gets the respective embedding network.
# def get_embedding_edge_dct(keyword):
#     assert keyword in ['full', 'sampled']
#     embedding_edge_dct = {}
#     f = open('./data/embedding_%s_network.txt' % keyword, 'r')
#     for line in f:
#         gene_a, gene_b, weight = line.split()
#         embedding_edge_dct[(gene_a, gene_b)] = weight
#     f.close()
#     return embedding_edge_dct

# gene_edge_weights.py
# dump_go_dictionary_files.py
# compute_go_enrichment.py
def get_high_std_genes(data_type):
    '''
    Retrieves the list of genes that have high standard deviations across their
    gene expression vectors.
    '''
    high_std_genes = []
    f = open('./data/%s_high_std_genes.txt' % data_type, 'r')
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
    f = open('./data/high_std_%s_network.txt' % data_type, 'r')
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
    num_options = 8
    config_dct = {}
    f = open('%s_config.txt' % data_type, 'r')
    for i, line in enumerate(f):
        config_num = i % (num_options + 1)
        if config_num == 0:
            run_num = line.split(' ')[1].strip()
            config_dct[run_num] = {}
        elif config_num == 1:
            subgraph_decimal = line.split()[0].strip()
            config_dct[run_num]['subgraph_decimal'] = subgraph_decimal
        elif config_num == 2:
            temp = line.split()[2].strip()
            config_dct[run_num]['temp'] = temp
        elif config_num == 3:
            min_go_size = line.split()[2].strip()
            config_dct[run_num]['min_go_size'] = int(min_go_size)
        elif config_num == 4:
            max_go_size = line.split()[2].strip()
            config_dct[run_num]['max_go_size'] = int(max_go_size)
        elif config_num == 5:
            edge_method = line.strip()
            config_dct[run_num]['edge_method'] = edge_method
        elif config_num == 6:
            lamb = line.split()[2].strip()
            config_dct[run_num]['lamb'] = float(lamb)
        elif config_num == 7:
            num_clusters = line.split()[1].strip()
            config_dct[run_num]['num_clusters'] = num_clusters
    f.close()
    return config_dct
