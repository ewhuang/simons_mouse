### Author: Edward Huang

import sys

### Reformats the output of WGCNA for running through simulated annealing
### evluation scripts. Takes out any GO terms.

def get_gene_to_index_dct(genes):
    '''
    Takes a list of genes, and returns a dictionary where keys are the genes,
    and the values are their corresponding indices in the lists.
    '''
    gene_to_index_dct = {}
    for index, gene in enumerate(genes):
        gene_to_index_dct[gene] = str(index)
    return gene_to_index_dct

def get_high_std_genes():
    '''
    Retrieves the list of genes that have high standard deviations across their
    gene expression vectors.
    '''
    high_std_genes = []
    f = open('../data/high_std_genes.txt', 'r')
    for line in f:
        gene = line.strip()        
        assert 'ENSMUSG' in gene
        high_std_genes += [gene]
    f.close()
    return high_std_genes

def main():
    if len(sys.argv) != 2:
        print 'Usage: %s genes_only/pca/mean/median' % sys.argv[0]
        exit()
    go_method = sys.argv[1]
    assert go_method in ['genes_only', 'pca', 'mean', 'median']

    high_std_genes = get_high_std_genes()
    # Map the gene ENSMUSG ID's to indices.
    gene_to_index_dct = get_gene_to_index_dct(high_std_genes)

    if go_method == 'genes_only':
        domain_list = ['genes_only']
    else:
        domain_list = ['bp', 'cc', 'mf']

    for go_domain in domain_list:
        f = open('./%s_results/module_membership_%s.txt' % (go_method,
            go_domain), 'r')
        out = open('./%s_results/clusters_%s.txt' % (go_method,
            go_domain), 'w')
        # Write dummy header line.
        out.write('dummy_header\n')

        for i, line in enumerate(f):
            # Skip header.
            if i == 0:
                continue

            # Ignore color and membershi columns.
            node, module, color, membership = line.split()

            # Skip garbage module.
            if module == '0':
                continue
            
            # Remove quotiation marks around the ENSMUSG ID.
            node = node.strip('"')

            # Don't add in GO terms.
            if 'GO:' in node:
                continue

            # Convert ENSMUSG ID to index.
            gene_index = gene_to_index_dct[node]

            # Write out the line.
            out.write('Species 0\tGene %s\tCluster %s\n' % (gene_index, module))
        out.close()
        f.close()

if __name__ == '__main__':
    main()