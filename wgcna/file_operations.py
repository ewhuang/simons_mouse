# preprocess_wgcna.py
# clean_wgcna_module_results.py
# evaluate_clustering_wgcna.py
# compute_label_enrichments.py
def get_tcga_disease_list():
    disease_list = []
    f = open('../data/tcga_data/tcga_diseases.txt', 'r')
    for line in f:
        disease_list += [line.strip()]
    f.close()
    return disease_list

# preprocess_wgcna.py
# compute_label_enrichments.py
def get_high_std_genes(data_type):
    '''
    Retrieves the list of genes that have high standard deviations across their
    gene expression vectors.
    '''
    high_std_genes = []
    f = open('../data/%s_data/high_std_genes.txt' % data_type, 'r')
    for line in f:
        gene = line.strip()
        assert 'ENSMUSG' in gene or 'ENSG' in gene
        high_std_genes += [gene]
    f.close()
    return high_std_genes