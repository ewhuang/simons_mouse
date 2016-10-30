### Author: Edward Huang

import sys
import time

### This script reads the gene association file from either TCGA or MGI data
### and then writes out a file of genes annotated by either BP terms MF terms.
### The output is uploaded to find ensembl analogs.
### Run time: a second.

def read_gene_association_file():
    '''
    Reads the gene association file, whether it's human or mouse, and outputs
    a list of the unique ID's that appear in the association file. Returns
    one list for the genes associated with MF terms, and one for BP terms.
    Returns MGI IDs for mouse genes and UniProtKB IDs for human genes.
    '''
    gene_list = []
    if gene_type == 'mouse':
        f = open('./data/mouse_data/gene_association.mgi', 'r')
    elif gene_type == 'tcga':
        f = open('./data/tcga_data/goa_human.gaf', 'r')
    for line in f:
        line = line.strip().split('\t')
        # Ninth column is which ontology the GO term belongs to.
        aspect = line[8]
        assert aspect in ['P', 'F', 'C']
        # Skip cellular component associations.
        if aspect == 'C':
            continue
        gene_id = line[1]
        gene_list += [gene_id]

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