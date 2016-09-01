### Author: Edward Huang

import json
import sys
import time

### This script creates the GO-GO dictionary. The GO terms are not GO ID's,
### but the actual names.

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

def dump_go_dct(go_domain, go_dct):
    '''
    Dump each dictionary out to file.
    '''
    with open('./data/%s_data/%s_dct.json' % (gene_type, go_domain),
        'w') as fp:
        json.dump(go_dct, fp)
    fp.close()

def test_go_go_edge_dct(go_go_edge_dct):
    '''
    Runs basic tests on the resulting GO-GO dictionary.
    '''
    assert ('organelle_inheritance' in
        go_go_edge_dct['mitochondrion_(inheritance'])
    assert ('_'.join(('response to light intensity').split()) in 
        go_go_edge_dct['_'.join(('response to absence of light').split())])
    assert ('S-acetyltransferase_activity' in go_go_edge_dct[
        'thioethanolamine_S-acetyltransferase_activity'])

def main():
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga' % sys.argv[0]
        exit()
    global gene_type
    gene_type = sys.argv[1]
    assert gene_type in ['mouse', 'tcga']

    go_go_edge_dct = get_go_go_edge_dct()
    dump_go_dct('mf_go_go', go_go_edge_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))