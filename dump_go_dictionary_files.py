### Author: Edward Huang

import file_operations
import json
import time

### Go to http://pantherdb.org/
### Upload genes with high standard deviation file. Select mus musculus.
### Customize Gene List: Add GO Database CC/BP/MF Complete to "Columns sorted
### using user preferences". Send list to file. It will download as
### pantherGeneList.txt.

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
        go_term = go_term[go_term.rfind("(") + 1:go_term.rfind(")")]
        assert len(go_term) == 10
        
        go_id_list += [go_term]

    return go_id_list

def add_to_dictionary(ensmusg_id, term_list, go_dct):
    '''
    Adds GO terms as keys, with the ENSMUSG ID as values to the dictionary.
    '''
    # Either this is the first time we add the gene to the dictionary, or there
    # is nothing meaningful to add.
    for term in term_list:
        if term in go_dct:
            go_dct[term] += [ensmusg_id]
        else:
            go_dct[term] = [ensmusg_id]

def get_go_domains():
    '''
    Returns three dictionaries, one corresponding to each GO domain. Keys are
    genes, and values are GO terms.
    '''
    high_std_genes = file_operations.get_high_std_genes()

    # cc = cellular component, bp = biological process, mf = molecular function.
    cc_go_gene_dct = {}
    bp_go_gene_dct = {}
    mf_go_gene_dct = {}

    f = open('./data/pantherGeneList.txt', 'r')
    for line in f:
        (gene_id, ensmusg_id_list, ortholog, cc_term_list, bp_term_list,
            mf_term_list) = line.split('\t')

        # Process the GO terms.
        cc_term_list = process_go_term_list(cc_term_list)
        bp_term_list = process_go_term_list(bp_term_list)
        mf_term_list = process_go_term_list(mf_term_list)

        # Associate each GO term list with the corresponding gene.
        ensmusg_id_list = ensmusg_id_list.split(',')
        for ensmusg_id in ensmusg_id_list:
            assert 'ENSMUSG' in ensmusg_id and ensmusg_id in high_std_genes
            add_to_dictionary(ensmusg_id, cc_term_list, cc_go_gene_dct)
            add_to_dictionary(ensmusg_id, bp_term_list, bp_go_gene_dct)
            add_to_dictionary(ensmusg_id, mf_term_list, mf_go_gene_dct)
    f.close()

    return cc_go_gene_dct, bp_go_gene_dct, mf_go_gene_dct

def main():
    cc_go_gene_dct, bp_go_gene_dct, mf_go_gene_dct = get_go_domains()

    # Dump each dictionary out to file.
    with open('./data/cellular_component_go.json', 'w') as fp:
        json.dump(cc_go_gene_dct, fp)
    fp.close()

    with open('./data/biological_process.json', 'w') as fp:
        json.dump(bp_go_gene_dct, fp)
    fp.close()

    with open('./data/molecular_function.json', 'w') as fp:
        json.dump(mf_go_gene_dct, fp)
    fp.close()

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))