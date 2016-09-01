### Author: Edward Huang

import json
import operator
from scipy.stats import fisher_exact
import sys
import time

### This script writes out a file showing overlapping GO terms as computed by
### Fisher's test. We write out all pairs of BP and MF terms that have a
### Fisher's p-value of less than 1e-10. We can write a file for mouse or TCGA
### data.
### Run time: 82 minutes.

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

def read_mouse_go_dictionaries():
    assert data_type == 'mouse'
    with open('./data/%s_data/bp_dct.json' % data_type, 'r') as fp:
        bp_go_gene_dct = json.load(fp)
    fp.close()

    with open('./data/%s_data/mf_dct.json' % data_type, 'r') as fp:
        mf_go_gene_dct = json.load(fp)
    fp.close()

    return (bp_go_gene_dct, mf_go_gene_dct)

def get_tcga_go_domains():
    '''
    Returns three dictionaries, one corresponding to each GO domain. Keys are
    genes, and values are GO terms.
    '''
    # bp = biological process, mf = molecular function.
    bp_go_gene_dct = {}
    mf_go_gene_dct = {}

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
            assert 'ENSG' in gene
            add_to_dictionary(gene, bp_term_list, bp_go_gene_dct)
            add_to_dictionary(gene, mf_term_list, mf_go_gene_dct)
    f.close()

    return bp_go_gene_dct, mf_go_gene_dct

def write_go_overlap(bp_go_gene_dct, mf_go_gene_dct):
    '''
    Given two dictionaries of GO terms, find the overlapping terms' Fisher's
    test p-values.
    '''
    gene_universe = set([item for sublist in bp_go_gene_dct.values() for item
        in sublist]).union([item for sublist in mf_go_gene_dct.values() for
        item in sublist])

    fisher_dct = {}
    for bp_label in bp_go_gene_dct:
        bp_genes = set(bp_go_gene_dct[bp_label])

        for mf_label in mf_go_gene_dct:
            mf_genes = set(mf_go_gene_dct[mf_label])

            # Compute the four sets for Fisher's test.
            bp_and_mf = len(bp_genes.intersection(mf_genes))
            bp_not_mf = len(bp_genes.difference(mf_genes))
            mf_not_bp = len(mf_genes.difference(bp_genes))
            neither = len(gene_universe) - len(mf_genes.union(bp_genes))

            # Run Fisher's test.
            f_table = ([[bp_and_mf, bp_not_mf], [mf_not_bp, neither]])
            o_r, p_value = fisher_exact(f_table, alternative='greater')

            # Handle overflow issues.
            p_value = max(p_value, 1e-300)

            fisher_dct[(bp_label, mf_label)] = p_value

    # Sort the fisher dictionary by p-value. Break when we hit a p-value of
    # 1e-10.
    fisher_dct = sorted(fisher_dct.items(), key=operator.itemgetter(1))
    out = open('./data/%s_data/overlapping_bp_mf_go_labels.txt' % data_type,
        'w')
    # First column BP, second column MF, third column p-value.
    for (bp_label, mf_label), p_value in fisher_dct:
        out.write(bp_label + '\t' + mf_label + '\t' + str(p_value) + '\n')
        if p_value > 1e-10:
            break
    out.close()

def main():    
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga' % sys.argv[0]
        exit()
    global data_type
    data_type = sys.argv[1]
    assert data_type in ['mouse', 'tcga']

    # If the data is mouse, then we just simply load the JSON files in the
    # directory.
    if data_type == 'mouse':
        bp_go_gene_dct, mf_go_gene_dct = read_mouse_go_dictionaries()
    # If the data type is TCGA, then we must create new dictionaries for it,
    # since we didn't previously create dictionaries for TCGA, but rather for
    # its subset cancers.
    elif data_type == 'tcga':
        bp_go_gene_dct, mf_go_gene_dct = get_tcga_go_domains()
    write_go_overlap(bp_go_gene_dct, mf_go_gene_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))