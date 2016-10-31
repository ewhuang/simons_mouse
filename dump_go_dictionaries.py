### Author: Edward Huang

import file_operations
import json
import sys
import time

### Dumps two dictionaries, a BP and an MF dictionary, as json files.
### Run time: 50 seconds for mouse, 30 minutes for TCGA. MF GO-GO dictionary
### takes less than a second.

def get_go_dictionaries(folder_name):
    '''
    folder_name is either 'mouse' or a TCGA cancer name.
    Reads the file downloaded from biomart, and creates a dictionary mapping
    GO terms to genes. Only considers genes that have high standard deviations
    across gene expression values.
    '''
    bp_dct, mf_dct = {}, {}

    def add_to_dictionary(go_dct, go_term, gene):
        # Augment a dictionary with a (go_term, gene) pair.
        if go_term in go_dct:
            if gene not in go_dct[go_term]:
                go_dct[go_term] += [gene]
        else:
            go_dct[go_term] = [gene]

    high_std_genes = file_operations.get_high_std_genes(folder_name)
    # Open the GO annotation file.
    if gene_type == 'mouse':
        f = open('./data/mouse_data/ensmusg_to_go.txt', 'r')
    else:
        f = open('./data/tcga_data/ensg_to_go.txt', 'r')
    for i, line in enumerate(f):
        # Skip header line.
        if i == 0:
            continue
        line = line.strip().split('\t')
        if len(line) != 3:
            continue
        ensembl_gene_id, go_term_name, go_domain = line
        # Skip cellular component GO labels.
        if (go_domain == 'cellular_component' or ensembl_gene_id not in
            high_std_genes):
            continue
        go_term_name = '_'.join(go_term_name.split())
        assert 'ENS' in ensembl_gene_id

        # Add the GO-gene relationship to the appropriate dictionary.
        if go_domain == 'biological_process':
            add_to_dictionary(bp_dct, go_term_name, ensembl_gene_id)
        else:
            assert go_domain == 'molecular_function'
            add_to_dictionary(mf_dct, go_term_name, ensembl_gene_id)
    f.close()
    return bp_dct, mf_dct

def get_go_go_edge_dct():
    '''
    Returns a dictionary mapping GO terms to their synonyms.
    '''    
    go_go_edge_dct = {}
    f = open('./data/go-basic.obo.txt', 'r')
    while True:
        line = f.readline()
        if line == '':
            break
        if line.strip() == '[Term]':
            # Keep track of the current term's GO id's.
            current_go_id_list = []

            go_id = f.readline().split()[1]
            go_term = '_'.join(f.readline().split()[1:])
            current_go_id_list += [go_term]
            namespace = f.readline()
            if ('cellular_component' in namespace) or ('biological_process'
                in namespace):
                continue
            assert ('namespace' in namespace)
            next = f.readline()
            while 'alt_id' in next:
                alt_id = next.split()[1]
                # current_go_id_list += [alt_id]
                next = f.readline()
            # Keep reading until we see an is_a relationship.
            neighbors = []
            while next.strip() != '':
                if 'is_a:' in next or 'part_of' in next:
                    if 'is_a' in next:
                        neighbor_GO = '_'.join(next.split()[3:])
                    else:
                        neighbor_GO = '_'.join(next.split()[4:])
                    neighbors += [neighbor_GO]
                next = f.readline()

            if neighbors == []:
                continue

            for go_id in current_go_id_list:
                if go_id in go_go_edge_dct:
                    go_go_edge_dct[go_id] += neighbors
                else:
                    go_go_edge_dct[go_id] = neighbors[:]
    f.close()
    return go_go_edge_dct

def dump_go_dct(file_name, go_dct):
    with open(file_name, 'w') as fp:
        json.dump(go_dct, fp)
    fp.close()

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga/mf_go_go' % sys.argv[0]
        exit()
    global gene_type
    gene_type = sys.argv[1]
    assert gene_type in ['mouse', 'tcga', 'mf_go_go']

    # Construct just the MF GO-GO edge dictionary.
    if gene_type == 'mf_go_go':
        go_go_edge_dct = get_go_go_edge_dct()
        dump_go_dct('./data/mf_go_go_dct.json', go_go_edge_dct)
    else:
        # Otherwise, make the normal GO dictionaries.
        if gene_type == 'mouse':
            folder_list = ['mouse']
        else:
            # If TCGA, make a set of GO dictionaries for each cancer type.
            folder_list = file_operations.get_tcga_disease_list()

        for folder_name in folder_list:
            bp_dct, mf_dct = get_go_dictionaries(folder_name)
            folder = './data/%s_data/' % folder_name
            dump_go_dct(folder + 'bp_dct.json', bp_dct)
            dump_go_dct(folder + 'mf_dct.json', mf_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))