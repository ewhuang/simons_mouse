### Author: Edward Huang


### This script preprocesses the original coexpression file to prepare for
### WGCNA. It looks at the networks we have previously built, finds all the
### genes used in those networks, and then preserves only the rows in the 
### original file that matches those genes.


if __name__ == '__main__':
    # First, extract all of the genes that appear in our randomly sampled
    # network.
    network_genes = set([])
    f = open('./data/network_no_go_20.txt', 'r')
    for i, line in enumerate(f):
        if i < 2:
            continue
        gene_a, gene_b, weight = line.split()
        network_genes.add(gene_a)
        network_genes.add(gene_b)
    f.close()

    # Next, write out the gene rows that appear in the network.
    f = open('./data/mm_mrsb_log2_expression.tsv', 'r')
    out = open('./data/mm_mrsb_log2_expression_sampled.tsv', 'w')
    for i, line in enumerate(f):
        if i == 0:
            out.write(line)
            continue
        gene = line.split()[0]
        if gene not in network_genes:
            continue
        out.write(line)
    out.close()
    f.close()