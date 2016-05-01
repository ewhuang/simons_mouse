### Author: Edward Huang

### This script preprocesses the original coexpression file to prepare for
### WGCNA. It looks at the networks we have previously built, finds all the
### genes used in those networks, and then preserves only the rows in the 
### original file that matches those genes.
### Run time: 10 seconds.

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
    # First, get the high standard deviation genes.
    high_std_genes = get_high_std_genes()

    # Next, write out the gene rows that appear in the network.
    f = open('../data/mm_mrsb_log2_expression.tsv', 'r')
    out = open('./data/mm_mrsb_log2_expression_high_std.tsv', 'w')
    for i, line in enumerate(f):
        # Directly write out the header file.
        if i == 0:
            out.write(line)
            continue
        # Skip genes not high_std_genes.
        gene = line.split()[0]
        if gene not in high_std_genes:
            continue

        out.write(line)
    out.close()
    f.close()

if __name__ == '__main__':
    main()