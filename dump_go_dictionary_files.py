### Author: Edward Huang

import file_operations
import json
import operator
from scipy.stats import fisher_exact
import time

### Go to http://pantherdb.org/
### Upload genes with high standard deviation file. Select mus musculus.
### Customize Gene List: Add GO Database CC/BP/MF Complete to "Columns sorted
### using user preferences". Send list to file. It will download as
### pantherGeneList.txt.
### Run time: 4 seconds.

def get_gene_to_index_dct(genes):
    '''
    Takes a list of genes, and returns a dictionary where keys are the genes,
    and the values are their corresponding indices in the lists.
    '''
    gene_to_index_dct = {}
    for index, gene in enumerate(genes):
        gene_to_index_dct[gene] = str(index)
    return gene_to_index_dct

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

def get_go_domains(gene_type):
    '''
    Returns three dictionaries, one corresponding to each GO domain. Keys are
    genes, and values are GO terms.
    '''
    high_std_genes = file_operations.get_high_std_genes()
    # Map the gene ENSMUSG ID's to indices.
    gene_to_index_dct = get_gene_to_index_dct(high_std_genes)

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
        for gene in ensmusg_id_list:
            assert 'ENSMUSG' in gene and gene in high_std_genes
            if gene_type == 'index':
                gene = gene_to_index_dct[gene]
            add_to_dictionary(gene, cc_term_list, cc_go_gene_dct)
            add_to_dictionary(gene, bp_term_list, bp_go_gene_dct)
            add_to_dictionary(gene, mf_term_list, mf_go_gene_dct)
    f.close()

    return cc_go_gene_dct, bp_go_gene_dct, mf_go_gene_dct

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
        if len(bp_genes) > 1000 or len(bp_genes) < 10:
            continue

        for mf_label in mf_go_gene_dct:
            mf_genes = set(mf_go_gene_dct[mf_label])
            if len(mf_genes) > 1000 or len(mf_genes) < 10:
                continue

            # Compute the four sets for Fisher's test.
            bp_and_mf = len(bp_genes.intersection(mf_genes))
            bp_not_mf = len(bp_genes.difference(mf_genes))
            mf_not_bp = len(mf_genes.difference(bp_genes))
            neither = len(gene_universe) - len(mf_genes.union(bp_genes))

            # Run Fisher's test.
            f_table = ([[bp_and_mf, bp_not_mf], [mf_not_bp, neither]])
            o_r, p_value = fisher_exact(f_table)

            # Handle overflow issues.
            p_value = max(p_value, 1e-300)

            fisher_dct[(bp_label, mf_label)] = p_value

    fisher_dct = sorted(fisher_dct.items(), key=operator.itemgetter(1))
    out = open('./data/overlapping_bp_mf_go_labels.txt', 'w')
    for (bp_label, mf_label), p_value in fisher_dct:
        out.write(bp_label + '\t' + mf_label + '\t' + str(p_value) + '\n')
        if p_value > 1e-10:
            break
    out.close()

def main():
    # Write out two types of values: ENSMUSG ID's and their indices in the
    # high standard deviation genes.
    for gene_type in ['index', 'ensmusg']:
        cc_go_gene_dct, bp_go_gene_dct, mf_go_gene_dct = get_go_domains(
            gene_type)

        # Dump each dictionary out to file.
        with open('./data/cc_%s.json' % gene_type, 'w') as fp:
            json.dump(cc_go_gene_dct, fp)
        fp.close()

        with open('./data/bp_%s.json' % gene_type, 'w') as fp:
            json.dump(bp_go_gene_dct, fp)
        fp.close()

        with open('./data/mf_%s.json' % gene_type, 'w') as fp:
            json.dump(mf_go_gene_dct, fp)
        fp.close()

    write_go_overlap(bp_go_gene_dct, mf_go_gene_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))