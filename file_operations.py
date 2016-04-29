### Author: Edward Huang

from collections import OrderedDict
import math

### This file contains functions that parse the data files and return the 
### data objects that we work with in our scripts.

# Returns dictionary where keys are genes, and values are gene expression
# vectors.
def get_gene_expression_dct():
    f = open('./data/mm_mrsb_log2_expression.tsv', 'r')
    gene_exp_dct = OrderedDict({})
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.split()
        gene, exp_vals = line[0], line[1:]
        exp_vals = [float(val) for val in exp_vals]
        gene_exp_dct[gene] = exp_vals
    f.close()
    return gene_exp_dct

# create_clustering_input.py
def get_high_std_edge_dct():
    '''
    Returns a dictionary of genes with high standard deviation. The keys are
    pairs of genes, and the values are the Pearson correlation weights between
    them.
    '''
    edge_dct = {}
    f = open('./data/high_std_network.txt', 'r')
    for i, line in enumerate(f):
        gene_a, gene_b, pcc = line.split()
        edge_dct[(gene_a, gene_b)] = pcc
    f.close()
    return edge_dct

# gene_edge_weights.py
def get_gene_to_index_dct(genes):
    '''
    Takes a list of genes, and returns a dictionary where keys are the genes,
    and the values are their corresponding indices in the lists.
    '''
    gene_to_index_dct = {}
    for index, gene in enumerate(genes):
        gene_to_index_dct[gene] = str(index)
    return gene_to_index_dct

# This function returns a dictionary, with keys as the names of GO annotations
# and values as lists of genes annotated by the keys. Genes are represented by
# their indices in the high_std_genes.txt file.
def get_go_labels(gene_set):
    high_std_genes = get_high_std_genes()
    gene_to_index_dct = get_gene_to_index_dct(high_std_genes)

    go_dct = {}
    f = open('./data/go_edges.txt', 'r')
    for i, line in enumerate(f):
        gene, go_label = line.split()
        
        if gene not in high_std_genes:
            continue
        # Map the name of the gene to the index in the array of high std genes.
        gene = gene_to_index_dct[gene]
        if gene not in gene_set:
            continue

        if go_label not in go_dct:
            go_dct[go_label] = [gene]
        else:
            go_dct[go_label] += [gene]
    f.close()
    return go_dct

# Takes a clustered file from simulated annealing and outputs a clean file
# without the GO nodes.
def create_clean_go_file(run_num):
    f = open('./results/clusters_go_%s.txt' % run_num, 'r')
    out = open('./results/clusters_go_clean_%s.txt' % run_num, 'w')
    for i, line in enumerate(f):
        if i == 0:
            out.write(line)
            continue
        newline = line.split()
        if 'ENSMUSG' not in newline[3]:
            continue
        out.write(line)
    out.close()
    f.close()

# Returns a dictionary where keys are cluster ID's, and values are lists of
# genes in the keys' clusters.
def create_cluster_dct(filename):
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
def get_all_genes():
    '''
    Reads the original co-expression matrix, and returns the set of genes that
    appear in the matrix.
    '''
    all_genes = set([])
    f = open('./data/mm_mrsb_log2_expression.tsv', 'r')
    gene_lst = []
    for i, line in enumerate(f):
        if i == 0:
            continue
        gene = line.split()[0]
        all_genes.add(gene)
    f.close()
    return all_genes

# make_gene_go_file.py
def get_gene_index_dct():
    '''
    Returns a dictionary that maps gene indices to their ENSMUSG id's.
    '''
    all_genes = get_all_genes()
    mgi_to_ensembl_dct = {}
    f = open('./go_edge_prediction/prediction_data/mgi_to_ensembl.txt', 'r')
    for i, line in enumerate(f):
        # Skip the header line.
        if i == 0:
            continue
        line = line.split()
        if len(line) != 5 or (line[1] != 'current' and line[1] != 'old'):
            continue
        ENSMUSG = line[4]
        # Skip if a gene isn't in our coexpression network, or if it isn't up
        # to date.
        if ENSMUSG not in all_genes:
            continue
        mgi_id = line[0]
        if mgi_id in mgi_to_ensembl_dct:
            mgi_to_ensembl_dct[mgi_id] += [ENSMUSG]
        else:
            mgi_to_ensembl_dct[mgi_id] = [ENSMUSG]
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
    return gene_index_dct

# Takes a set of GO terms and then separates them into their respective
# categories: biological process, cellular component, or molecular function.
def chunkify_go_terms(go_terms):
    go_id_to_name_dct = get_go_id_to_name_dct()

    # One chunk for each category.
    go_chunks = [[], [], []]
    f = open('./data/GO.namespace', 'r')
    for line in f:
        go_id, category = line.split()
        go_name = go_id_to_name_dct[go_id.lower()]
        if go_name not in go_terms:
            continue
        category = int(category) - 1
        go_chunks[category] += [go_name]
    f.close()
    return go_chunks

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
            go_id = f.readline().split()[1]
            go_name = '_'.join(f.readline()[len('name: '):].split())
            go_id_to_name_dct[go_id.lower()] = go_name
            f.readline()
            next = f.readline()
            while 'alt_id' in next:
                go_id_to_name_dct[next.split()[1].lower()] = go_name
                next = f.readline()
    f.close()
    return go_id_to_name_dct

# make_gene_go_file.py
def get_go_index_dct():
    '''
    This function returns a dictionary where keys are GO id's, and values are
    their corresponding indices.
    '''
    go_id_to_name_dct = get_go_id_to_name_dct()
    # Keys are the GO ids, values are the indices in the edge weight matrix.
    go_index_dct = {}
    f = open('./go_edge_prediction/prediction_data/noisogoHash.txt', 'r')
    for line in f:
        go_id, index = line.split()
        # Subtract 1 to change to list indices.
        go_index_dct[int(index) - 1] = go_id_to_name_dct[go_id]
    f.close()
    return go_index_dct

# Get the output files of the Perl script evaluate_clustering.pl and find the
# in-density and out-density of each cluster.
def get_cluster_evaluation_densities(eval_fname):
    dens_dct = {}
    f = open(eval_fname, 'r')
    for i, line in enumerate(f):
        if line[:7] != 'Cluster':
            continue
        line = line.split()
        clus_id, in_density, out_density = line[1], line[7], line[9]
        dens_dct[clus_id] = (float(in_density), float(out_density))
    f.close()
    return dens_dct

# Find the best p-value GO enrichments for each cluster. Returns a dictionary
# where keys are cluster ID's, and values are the best p-values.
def get_best_enrichment_dct(enrichment_fname):
    best_enrichment_dct = {}
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
            best_enrichment_dct[cid] = line[0]
    f.close()
    return best_enrichment_dct

# Read in the network to figure out how many GO terms went into the clustering.
def get_network_stats(network_fname):
    num_genes_net = 0
    num_gg_net = 0
    num_ggo_net = 0
    f = open(network_fname, 'r')
    edge_list_go = []
    for i, line in enumerate(f):
        if i == 0:
            continue
        if i == 1:
            num_genes_net = line.strip()
            continue
        # If not first two lines, find out if each edge is G-GO or G-G.
        node_a, node_b, weight = line.strip().split('\t')
        if 'ENSMUSG' not in node_a or 'ENSMUSG' not in node_b:
            num_ggo_net += 1
        else:
            num_gg_net += 1
        edge_list_go += [(node_a, node_b)]
    f.close()
    # Divide the two numbers by two to account for each edge in twice.
    assert(num_ggo_net % 2 == 0)
    assert(num_ggo_net % 2 == 0)
    num_ggo_net /= 2
    num_gg_net /= 2
    return num_genes_net, num_gg_net, num_ggo_net, edge_list_go

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
def get_high_std_genes():
    '''
    Retrieves the list of genes that have high standard deviations across their
    gene expression vectors.
    '''
    high_std_genes = []
    f = open('./data/high_std_genes.txt', 'r')
    for line in f:
        high_std_genes += [line.strip()]
    f.close()
    return high_std_genes

# Returns dictionary. Keys are run_num strings, values are dicts of config
# options. Each dct has key of subgraph_decimal, temp, min_go_size,
# max_go_size, pearson/embedding edges, lambda, num_clusters.
def read_config_file():
    NUM_OPTIONS = 8
    config_dct = {}
    f = open('config.txt', 'r')
    for i, line in enumerate(f):
        config_num = i % (NUM_OPTIONS + 1)
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
            config_dct[run_num]['min_go_size'] = min_go_size
        elif config_num == 4:
            max_go_size = line.split()[2].strip()
            config_dct[run_num]['max_go_size'] = max_go_size
        elif config_num == 5:
            edge_method = line.strip()
            config_dct[run_num]['edge_method'] = edge_method
        elif config_num == 6:
            lamb = line.split()[2].strip()
            config_dct[run_num]['lamb'] = lamb
        elif config_num == 7:
            num_clusters = line.split()[1].strip()
            config_dct[run_num]['num_clusters'] = num_clusters
    f.close()
    return config_dct

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