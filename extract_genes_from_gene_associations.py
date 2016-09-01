### Author: Edward Huang

import sys
import time

### This script reads the gene association file from either TCGA or MGI data
### and then writes out a list for the gene sets for BP terms and for MF terms.

def read_gene_association_file():
    '''
    Reads the gene association file, whether it's human or mouse, and outputs
    a list of the unique ID's that appear in the association file. Returns
    one list for the genes associated with MF terms, and one for BP terms.
    Returns MGI ID's for mouse genes.
    '''
    gene_list = []
    if gene_type == 'mouse':
        f = open('./data/mouse_data/gene_association.mgi', 'r')
    for line in f:
        line = line.strip().split('\t')
        # Ninth column is which ontology the GO term belongs to.
        aspect = line[8]
        assert aspect in ['P', 'F', 'C']
        # Skip cellular process associations.
        if aspect == 'C':
            continue
        # Second column is the MGI unique identifier.
        mgi_gene_id = line[1]
        assert 'MGI:' in mgi_gene_id
        # Skip cellular component associations.
        gene_list += [mgi_gene_id]

    f.close()
    return gene_list

def write_gene_list(gene_list):
    '''
    Takes in a list of genes, and the aspect (ontology) to which the genes 
    belong. Writes them out to file.
    '''
    out = open('./data/%s_data/consortium_gene_list.txt' % gene_type, 'w')
    for gene in gene_list:
        out.write('%s\n' % gene)
    out.close()

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga' % sys.argv[0]
        exit()
    global gene_type
    gene_type = sys.argv[1]
    assert gene_type in ['mouse', 'tcga']

    gene_list = read_gene_association_file()
    write_gene_list(gene_list)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))