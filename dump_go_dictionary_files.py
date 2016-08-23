### Author: Edward Huang

import file_operations
import json
import operator
from scipy.stats import fisher_exact
import sys
import time

### Go to http://pantherdb.org/
### Upload genes with high standard deviation file. Select mus musculus.
### Customize Gene List: Add GO Database CC/BP/MF Complete to "Columns sorted
### using user preferences". Send list to file. It will download as
### pantherGeneList.txt. Rename with a _mouse or _tcga suffix.
### Run time: 3 seconds.

def process_go_term_list(go_term_list):
    '''
    Takes as input a string of GO terms, and processes them into a list of GO
    ID's.
    '''
    go_id_list = []
    
    # GO terms are separated by semicolons.
    go_term_list = go_term_list.strip().split(';')
    
    # Some genes aren't labeled by a specific domain.
    if go_term_list == ['']:
        return go_id_list

    for go_term in go_term_list:
        assert 'GO:' in go_term
        
        # Extract the GO ID between parentheses.
        go_term = go_term[:go_term.index('(')]
        go_term = '_'.join(go_term.split())
        go_id_list += [go_term]

    return go_id_list

def add_to_dictionary(gene, term_list, go_dct):
    '''
    Adds GO terms as keys, with the ENSMUSG ID as values to the dictionary.
    '''
    # Either this is the first time we add the gene to the dictionary, or there
    # is nothing meaningful to add.
    for term in term_list:
        if term in go_dct:
            go_dct[term] += [gene]
        else:
            go_dct[term] = [gene]

def get_go_domains(data_type):
    '''
    Returns three dictionaries, one corresponding to each GO domain. Keys are
    genes, and values are GO terms.
    '''
    high_std_genes = file_operations.get_high_std_genes(data_type)

    # bp = biological process, mf = molecular function.
    bp_go_gene_dct = {}
    mf_go_gene_dct = {}

    if data_type == 'mouse':
        f = open('./data/mouse_data/pantherGeneList_mouse.txt', 'r')
    else:
        f = open('./data/tcga_data/pantherGeneList_tcga.txt', 'r')

    for line in f:
        (gene_id, ensembl_id_list, ortholog, family, protein_class, species,
            bp_term_list, mf_term_list) = line.split('\t')

        # Process the GO terms.
        bp_term_list = process_go_term_list(bp_term_list)
        mf_term_list = process_go_term_list(mf_term_list)

        # Associate each GO term list with the corresponding gene.
        ensembl_id_list = ensembl_id_list.split(',')
        for gene in ensembl_id_list:
            assert ('ENSMUSG' in gene or 'ENSG' in gene)
            if gene not in high_std_genes:
                continue
            add_to_dictionary(gene, bp_term_list, bp_go_gene_dct)
            add_to_dictionary(gene, mf_term_list, mf_go_gene_dct)
    f.close()

    return bp_go_gene_dct, mf_go_gene_dct

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga_cancers' % sys.argv[0]
        exit()
    data_type = sys.argv[1]    
    assert data_type == 'mouse' or data_type.isdigit()
    if data_type.isdigit():
        data_type = file_operations.get_tcga_disease_list()[int(data_type)]

    bp_go_gene_dct, mf_go_gene_dct = get_go_domains(data_type)

    # Delete empty keys.
    if '' in bp_go_gene_dct:
        del bp_go_gene_dct['']
    if '' in mf_go_gene_dct:
        del mf_go_gene_dct['']

    # Dump each dictionary out to file.
    with open('./data/%s_data/bp_dct.json' % data_type, 'w') as fp:
        json.dump(bp_go_gene_dct, fp)
    fp.close()

    with open('./data/%s_data/mf_dct.json' % data_type, 'w') as fp:
        json.dump(mf_go_gene_dct, fp)
    fp.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))