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
### Run time: 80 minutes.

# This variable determines two GO terms that are "too correlated".
P_THRESHOLD = 1e-10
gene_universe = set([])

def get_go_dictionaries():
    '''
    This function should be identical to the one in dump_go_dictionaries.py,
    except for the fact that this one does __not__ exclude genes not in the
    high standard deviation list. Another small difference is the recording of
    gene_universe.
    '''
    global gene_universe
    bp_dct, mf_dct = {}, {}

    def add_to_dictionary(go_dct, go_term, gene):
        # Augment a dictionary with a (go_term, gene) pair.
        if go_term in go_dct:
            go_dct[go_term] += [gene]
        else:
            go_dct[go_term] = [gene]

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
        if go_domain == 'cellular_component':
            continue
        gene_universe.add(ensembl_gene_id)
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

def write_go_overlap(bp_dct, mf_dct):
    '''
    Write out a file containing BP terms that overlap with MF terms.
    Format:
    BP\tMF\tFisher's p-value\n
    '''
    fisher_dct = {}
    for bp_label in bp_dct:
        bp_genes = set(bp_dct[bp_label])

        for mf_label in mf_dct:
            mf_genes = set(mf_dct[mf_label])

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

            if p_value <= P_THRESHOLD:
                fisher_dct[(bp_label, mf_label)] = p_value

    # Sort the fisher dictionary by p-value. Break when we hit a p-value of
    # 1e-10.
    fisher_dct = sorted(fisher_dct.items(), key=operator.itemgetter(1))
    out = open('./data/%s_data/overlapping_bp_mf_go_labels.txt' % gene_type,
        'w')

    # First column BP, second column MF, third column p-value.
    for (bp_label, mf_label), p_value in fisher_dct:
        out.write('%s\t%s\t%g\n' % (bp_label, mf_label, p_value))
    out.close()

def main():    
    if len(sys.argv) != 2:
        print 'Usage:python %s mouse/tcga' % sys.argv[0]
        exit()
    global gene_type
    gene_type = sys.argv[1]
    assert gene_type in ['mouse', 'tcga']

    # Get the unfiltered GO dictionaries from the desired dataset.
    bp_dct, mf_dct = get_go_dictionaries()

    write_go_overlap(bp_dct, mf_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))