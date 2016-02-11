### Author: Edward Huang

### This script runs through the outputs of the WGCNA scripts and tries to
### output files that are in the same format as the simulated annealing cluster
### results.


if __name__ == '__main__':
    # Find out what genes to keep.
    f = open('./data/WGCNA_data/good_gene_booleans_WGCNA.txt', 'r')
    good_gene_indices = []
    for line in f:
        good_gene_indices += [line.strip()]
    f.close()

    # Get the actual genes corresponding to the genes kept by WGCNA.
    good_genes_lst = []
    f = open('./data/mm_mrsb_log2_expression.tsv', 'r')
    for i, line in enumerate(f):
        # Skip header.
        if i == 0:
            continue
        # Minus 1 because we didn't have a header for gg file.
        gg_lst_index = i - 1
        gg_bool = good_gene_indices[gg_lst_index]
        if gg_bool == 'TRUE':
            current_gene = line.split()[0]
            good_genes_lst += [current_gene]
    f.close()

    # Find the module memberships for the relevant genes.
    gene_cluster_dct = {}
    module_lst = []
    f = open('./data/WGCNA_data/module_membership_WGCNA.txt', 'r')
    for i, line in enumerate(f):
        # Skip header.
        if i == 0:
            continue
        # Columns are module, color, and module_membership. We ignore the last
        # two.
        module, color, membership = line.split()
        if module == '0':
            continue
        # Because we have a header.
        good_gene_index = i - 1
        gene = good_genes_lst[good_gene_index]
        # Make sure each gene is only in one cluster at most.
        assert gene not in gene_cluster_dct
        gene_cluster_dct[gene] = module
        if module not in module_lst:
            module_lst += [module]
    f.close()

    # Extract all of the genes that appear in our randomly sampled network.
    sampled_network_genes = set([])
    f = open('./data/network_no_go_20.txt', 'r')
    for i, line in enumerate(f):
        if i < 2:
            continue
        gene_a, gene_b, weight = line.split()
        sampled_network_genes.add(gene_a)
        sampled_network_genes.add(gene_b)
    f.close()

    # Write an output file in the style of the simulated annealing.
    out = open('./results/WGCNA_results/WGCNA_clusters_all_genes.txt', 'w')
    out.write('header\n')
    for gene in gene_cluster_dct:
        if gene not in sampled_network_genes:
            continue
        cluster_num = module_lst.index(gene_cluster_dct[gene]) + 1
        out.write('Species 0\tGene %s\tCluster %s\n' % (gene, cluster_num))
    out.close()