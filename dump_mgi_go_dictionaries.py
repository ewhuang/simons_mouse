### Author: Edward Huang

import json
import sys
import time

### This script writes out new json dictionary files for GO terms, as given
### by the official GO consortium. We must first map from MGI to ENSMUSG.
### Run time: 4 seconds.

def get_mgi_to_ensembl_dct():
    '''
    Reads the mapping file as created by informatics.jax.org and then returns
    a dictionary. An MGI might map to multiple ensembl ID's.
    Key: MGI ID (MGI:xxxx) -> str
    Value: list of ENSMUSG ID's (ENSMUGXXXXXX) -> list(str)
    '''
    mgi_to_ensembl_dct = {}
    if gene_type == 'mouse':
        f = open('./data/mouse_data/mgi_to_ensembl.txt', 'r')
    for i, line in enumerate(f):
        # Skip header line.
        if i == 0:
            continue
        line = line.split()

        # Skip lines that do not have ensembl mappings.
        if len(line) != 4:
            continue
        
        mgi_gene = line[2]
        assert 'MGI:' in mgi_gene
        ensembl_gene = line[3]
        assert 'ENSMUSG' in ensembl_gene

        if mgi_gene in mgi_to_ensembl_dct:
            mgi_to_ensembl_dct[mgi_gene] += [ensembl_gene]
        else:
            mgi_to_ensembl_dct[mgi_gene] = [ensembl_gene]
    f.close()

    return mgi_to_ensembl_dct

def add_to_dictionary(go_dct, go_term, gene_id_list):
    '''
    Adds GO terms as keys, with the ENSMUSG ID as values to the dictionary.
    '''
    if go_term in go_dct:
        go_dct[go_term] = go_dct[go_term].union(gene_id_list)
    else:
        go_dct[go_term] = set(gene_id_list)

def read_gene_association_file():
    '''
    Reads the gene association file, whether it's human or mouse, and outputs
    two dictionaries: one for BP and one for MF.
    Key: GO term name (GO:xxxxxx) -> str
    Value: list of ENSMUSG ID's (ENSMUSGXXXXXX) -> list(str)
    '''
    mgi_to_ensembl_dct = get_mgi_to_ensembl_dct()

    bp_go_dct, mf_go_dct = {}, {}
    if gene_type == 'mouse':
        f = open('./data/mouse_data/gene_association.mgi', 'r')
    for i, line in enumerate(f):
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
        if mgi_gene_id not in mgi_to_ensembl_dct:
            continue
        ensembl_gene_id_list = mgi_to_ensembl_dct[mgi_gene_id][:]
        go_term = line[4]
        assert 'GO:' in go_term
        # Skip cellular component associations.
        if aspect == 'P':
            add_to_dictionary(bp_go_dct, go_term, ensembl_gene_id_list)
        elif aspect =='F':
            add_to_dictionary(mf_go_dct, go_term, ensembl_gene_id_list)
    f.close()
    return bp_go_dct, mf_go_dct

def convert_set_values_to_lists(go_dct):
    '''
    We convert the genes sets to lists because sets are not JSON serializable.
    '''
    for go_term in go_dct:
        go_dct[go_term] = list(go_dct[go_term])

def dump_go_dct(go_domain, go_dct):
    '''
    Dump each dictionary out to file.
    '''
    with open('./data/%s_data/%s_dct.json' % (gene_type, go_domain),
        'w') as fp:
        json.dump(go_dct, fp)
    fp.close()

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga' % sys.argv[0]
        exit()
    global gene_type
    gene_type = sys.argv[1]
    assert gene_type in ['mouse', 'tcga']

    bp_go_dct, mf_go_dct = read_gene_association_file()

    go_go_edge_dct = get_go_go_edge_dct()

    dump_go_dct('mf_go_go', go_go_edge_dct)

    convert_set_values_to_lists(bp_go_dct)
    convert_set_values_to_lists(mf_go_dct)
    
    dump_go_dct('bp', bp_go_dct)
    dump_go_dct('mf', mf_go_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))