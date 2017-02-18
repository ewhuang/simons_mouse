### Author: Edward Huang

import file_operations
import json
import sys
import time

### Dumps two dictionaries, a BP and an MF dictionary, as json files.
### Run time: 50 seconds for mouse, 15 minutes for TCGA. MF GO-GO dictionary
### takes 5 seconds. Takes less time for DBGAP dictionaries.

def get_go_dictionaries(gene_type, network_genes):
    '''
    gene_type is either 'mouse' or a TCGA cancer name.
    Reads the file downloaded from biomart, and creates a dictionary mapping
    GO terms to genes. Only considers genes that are in the network.
    '''
    bp_dct, mf_dct = {}, {}

    def add_to_dictionary(go_dct, go_term, gene):
        # Augment a dictionary with a (go_term, gene) pair.
        if go_term in go_dct:
            if gene not in go_dct[go_term]:
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
        if (go_domain == 'cellular_component' or ensembl_gene_id not in
            network_genes):
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
    # Remove terms that are too large or small from the MF dictionary.
    bad_mf_terms = []
    for term in mf_dct:
        term_size = len(mf_dct[term])
        if term_size > 1000 or term_size < 10:
            bad_mf_terms += [term]
    for term in bad_mf_terms:
        del mf_dct[term]
    return bp_dct, mf_dct

def get_go_go_edge_dct():
    '''
    Returns a dictionary mapping GO terms to their synonyms. Get the file from
    http://geneontology.org/page/download-ontology
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

def get_ensg_to_ensmusg_dct():
    '''
    Gets a dictionary mapping ENSG ID's to their mouse homologs.
    Key: ENSG ID -> str
    Value: list of ENSMUSG IDs -> list(str)
    '''
    ensg_to_ensmusg_dct = {}
    f = open('./data/mouse_data/esng_to_ensmusg.txt', 'r')
    for i, line in enumerate(f):
        # Skip header.
        if i == 0:
            continue
        ensg_id, ensmusg_id = line.split()
        if ensg_id not in ensg_to_ensmusg_dct:
            ensg_to_ensmusg_dct[ensg_id] = []
        ensg_to_ensmusg_dct[ensg_id] += [ensmusg_id]
    f.close()
    return ensg_to_ensmusg_dct

def get_external_to_ensg_dct(fname):
    '''
    Get the mappings from Entrez Gene identifiers to Ensembl.
    '''
    external_to_ensg_dct = {}
    f = open('./data/tcga_data/%s.txt' % fname, 'r')
    for i, line in enumerate(f):
        if i == 0:
            continue
        line = line.split()
        if len(line) != 2:
            continue
        external_id = line[1]
        if external_id not in external_to_ensg_dct:
            external_to_ensg_dct[external_id] = []
        external_to_ensg_dct[external_id] += [line[0]]
    f.close()
    return external_to_ensg_dct

def read_annotation_file(gene_type, network_genes):
    '''
    Gets the disease annotations from the disgenet website.
    '''
    if gene_type == 'mouse':
        ensg_to_ensmusg_dct = get_ensg_to_ensmusg_dct()

    # DBGAP uses ENSG IDs already.
    if label_type == 'nci':
        external_to_ensg_dct = get_external_to_ensg_dct('ensg_to_hgnc')
    elif label_type in ['gwas', 'kegg', 'ctd', 'nci']:
        external_to_ensg_dct = get_external_to_ensg_dct('ensg_to_entrez')

    label_to_gene_dct = {}

    # TODO: Change the database of associations.
    if label_type == 'gwas':
        f = open('./data/curated_gene_disease_associations.tsv', 'r')
        # f = open('./data/all_gene_disease_associations.tsv', 'r')
        # f = open('./data/all_snps_sentences_pubmeds_position.tsv', 'r')
    elif label_type == 'kegg':
        f = open('./data/CTD_genes_pathways.tsv', 'r')
    elif label_type == 'ctd':
        f = open('./data/CTD_genes_diseases.tsv', 'r')
    elif label_type == 'nci':
        f = open('./data/nci_pathway_hgnc.txt', 'r')
    elif label_type == 'dbgap':
        f = open('./data/dbgap.txt', 'r')
    
    for line in f:
        # There are a few header lines.
        if label_type == 'gwas' and 'umls' not in line:
            continue
        elif label_type == 'kegg' and 'KEGG:' not in line:
            continue
        elif label_type == 'ctd' and 'MESH:' not in line:
            continue
        elif label_type == 'dbgap' and ('LYMPHOCYTES' in line or 'DNA' in line):
            continue
        line = line.strip().split('\t')

        # TODO: This line is for curated associations.
        if label_type == 'gwas':
            external_id, label, score = line[1], line[5], float(line[2])
            # This line is for all associations.
            # external_id, label, score = line[0], line[4], float(line[5])
            # This line is for SNP associations.
            # external_id, label, score = line[2], line[5], float(line[8])
            # This line is for all/curated associations.
            assert len(line) == 9
            # This line is for SNP associations.
            # assert external_id.isdigit() and len(line) == 15
            # TODO: Tune the score < THRESHOLD boolean.
            if score < 0.001:
                continue
        elif label_type == 'kegg':
            gene_symbol, external_id, label, path_id = line
        elif label_type == 'ctd':
            external_id, label = line[1], line[2]
        elif label_type == 'nci':
            label, external_id = line
        elif label_type == 'dbgap':
            label, external_id = line[0], line[1]

        if label_type not in ['nci', 'dbgap']:
            assert external_id.isdigit()

        if label_type != 'dbgap' and external_id not in external_to_ensg_dct:
            continue

        if label_type == 'dbgap':
            ensembl_gene_id_list = [external_id]
        else:
            ensembl_gene_id_list = external_to_ensg_dct[external_id]
        if gene_type == 'mouse':
            # Convert human to mouse homolog list.
            ensmusg_list = []
            for ensg_id in ensembl_gene_id_list:
                if ensg_id in ensg_to_ensmusg_dct:
                    ensmusg_list += ensg_to_ensmusg_dct[ensg_id]
            ensembl_gene_id_list = ensmusg_list
        # Failed to translate to mouse genes.
        if ensembl_gene_id_list == []:
            continue

        # Replace spaces and tabs with underscores.
        if label_type != 'dbgap':
            label = '_'.join(label.split())

        for ensembl_gene_id in ensembl_gene_id_list:
            if ensembl_gene_id not in network_genes:
                continue
            if label not in label_to_gene_dct:
                label_to_gene_dct[label] = []
            if ensembl_gene_id not in label_to_gene_dct[label]:
                label_to_gene_dct[label] += [ensembl_gene_id]
    f.close()
    return label_to_gene_dct

def dump_label_dct(file_name, label_dct):
    with open(file_name, 'w') as fp:
        json.dump(label_dct, fp)
    fp.close()

def main():
    if len(sys.argv) != 3:
        print ('Usage:python %s mouse/tcga/all/mf_go_go go/dbgap/gwas/kegg/ctd'
            '/nci') % sys.argv[0]
        exit()
    global label_type
    data_type, label_type = sys.argv[1:]
    assert data_type in ['mouse', 'tcga', 'all', 'mf_go_go']
    assert label_type in ['go', 'dbgap', 'gwas', 'kegg', 'ctd', 'nci']

    # Construct just the MF GO-GO edge dictionary.
    if data_type == 'mf_go_go':
        go_go_edge_dct = get_go_go_edge_dct()
        dump_label_dct('./data/mf_go_go_dct.json', go_go_edge_dct)
    else:
        gene_type_list = [data_type]
        if data_type == 'all':
            # If TCGA, make a set of GO dictionaries for each cancer type.
            gene_type_list = file_operations.get_tcga_list()

        for gene_type in gene_type_list:
            network_genes = file_operations.get_network_genes(gene_type)
            folder = './data/%s_data' % gene_type
            if label_type == 'go':
                bp_dct, mf_dct = get_go_dictionaries(gene_type, network_genes)
                dump_label_dct('%s/bp_dct.json' % folder, bp_dct)
                dump_label_dct('%s/mf_dct.json' % folder, mf_dct)
            elif label_type in ['gwas', 'kegg', 'ctd', 'nci', 'dbgap']:
                label_dct = read_annotation_file(gene_type, network_genes)
                dump_label_dct('%s/%s_dct.json' % (folder, label_type),
                    label_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))