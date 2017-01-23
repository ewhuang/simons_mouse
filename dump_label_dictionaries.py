### Author: Edward Huang

import file_operations
import json
import sys
import time

### Dumps two dictionaries, a BP and an MF dictionary, as json files.
### Run time: 50 seconds for mouse, 15 minutes for TCGA. MF GO-GO dictionary
### takes 5 seconds. Takes less time for DBGAP dictionaries.

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

def read_dbgap_file(folder_name):
    '''
    Gets the DBGAP dictionary. Maps a dbgap ID to a list of genes.
    Key: DBGAP ID -> str
    Value: list of ENSMUSG IDs -> list(str)
    '''
    if folder_name == 'mouse':
        ensg_to_ensmusg_dct = get_ensg_to_ensmusg_dct()
    high_std_genes = file_operations.get_high_std_genes(folder_name)

    dbgap_to_gene_dct = {}
    f = open('./data/dbgap.txt', 'r')
    for line in f:
        dbgap_id, ensg_id, bloat_1, bloat_2 = line.split()

        # ENSG values are single genes.
        ensembl_gene_id_list = [ensg_id]
        if folder_name == 'mouse':
            # Convert human to mouse homolog list.
            ensembl_gene_id_list = []
            if ensg_id in ensg_to_ensmusg_dct:
                ensembl_gene_id_list = ensg_to_ensmusg_dct[ensg_id]
        # Failed to convert to mouse genes.
        if ensembl_gene_id_list == []:
            continue

        for ensembl_gene_id in ensembl_gene_id_list:
            if ensembl_gene_id not in high_std_genes:
                continue
            if dbgap_id not in dbgap_to_gene_dct:
                dbgap_to_gene_dct[dbgap_id] = []
            if ensembl_gene_id not in dbgap_to_gene_dct[dbgap_id]:
                dbgap_to_gene_dct[dbgap_id] += [ensembl_gene_id]
    f.close()
    return dbgap_to_gene_dct

def read_gwas_file(folder_name):
    '''
    Gets the disease annotations from the disgenet website.
    '''
    def get_entrez_to_ensg_dct():
        '''
        Get the mappings from Entrez Gene identifiers to Ensembl.
        '''
        entrez_to_ensg_dct = {}
        f = open('./data/tcga_data/entrez_to_ensg.txt', 'r')
        for i, line in enumerate(f):
            if i == 0:
                continue
            line = line.split()
            if len(line) != 2:
                continue
            entrez_id = line[1]
            if entrez_id not in entrez_to_ensg_dct:
                entrez_to_ensg_dct[entrez_id] = []
            entrez_to_ensg_dct[entrez_id] += [line[0]]
        f.close()
        return entrez_to_ensg_dct

    if folder_name == 'mouse':
        ensg_to_ensmusg_dct = get_ensg_to_ensmusg_dct()

    entrez_to_ensg_dct = get_entrez_to_ensg_dct()
    high_std_genes = file_operations.get_high_std_genes(folder_name)

    gwas_to_gene_dct = {}
    
    # TODO: Change the database of associations.
    f = open('./data/curated_gene_disease_associations.tsv', 'r')
    # f = open('./data/all_gene_disease_associations.tsv', 'r')
    # f = open('./data/all_snps_sentences_pubmeds_position.tsv', 'r')
    
    for line in f:
        # There are a few header lines.
        if 'umls' not in line:
            continue
        line = line.strip().split('\t')

        # TODO: This line is for curated associations.
        entrez_id, disease_label, score = line[1], line[5], float(line[2])
        # This line is for all associations.
        # entrez_id, disease_label, score = line[0], line[4], float(line[5])
        # This line is for SNP associations.
        # entrez_id, disease_label, score = line[2], line[5], float(line[8])
        
        # This line is for all/curated associations.
        assert entrez_id.isdigit() and len(line) == 9
        # This line is for SNP associations.
        # assert entrez_id.isdigit() and len(line) == 15
        
        # TODO: Tune the score < THRESHOLD boolean.
        if entrez_id not in entrez_to_ensg_dct or score < 0.001:
            continue

        ensembl_gene_id_list = entrez_to_ensg_dct[entrez_id]
        if folder_name == 'mouse':
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
        disease_label = '_'.join(disease_label.split())

        for ensembl_gene_id in ensembl_gene_id_list:
            if ensembl_gene_id not in high_std_genes:
                continue
            if disease_label not in gwas_to_gene_dct:
                gwas_to_gene_dct[disease_label] = []
            if ensembl_gene_id not in gwas_to_gene_dct[disease_label]:
                gwas_to_gene_dct[disease_label] += [ensembl_gene_id]
    f.close()
    return gwas_to_gene_dct

def dump_label_dct(file_name, label_dct):
    with open(file_name, 'w') as fp:
        json.dump(label_dct, fp)
    fp.close()

def main():
    if len(sys.argv) != 3:
        print 'Usage:python %s mouse/tcga/mf_go_go go/dbgap/gwas' % sys.argv[0]
        exit()
    global gene_type, label_type
    gene_type, label_type = sys.argv[1:]
    assert gene_type in ['mouse', 'tcga', 'mf_go_go']
    assert label_type in ['go', 'dbgap', 'gwas']

    # Construct just the MF GO-GO edge dictionary.
    if gene_type == 'mf_go_go':
        go_go_edge_dct = get_go_go_edge_dct()
        dump_label_dct('./data/mf_go_go_dct.json', go_go_edge_dct)
    else:
        # Otherwise, make the normal GO dictionaries.
        if gene_type == 'mouse':
            folder_list = ['mouse']
        else:
            # If TCGA, make a set of GO dictionaries for each cancer type.
            folder_list = file_operations.get_tcga_disease_list()

        for folder_name in folder_list:
            folder = './data/%s_data/' % folder_name
            if label_type == 'go':
                bp_dct, mf_dct = get_go_dictionaries(folder_name)
                dump_label_dct(folder + 'bp_dct.json', bp_dct)
                dump_label_dct(folder + 'mf_dct.json', mf_dct)
            elif label_type == 'dbgap':
                dbgap_dct = read_dbgap_file(folder_name)
                dump_label_dct(folder + 'dbgap_dct.json', dbgap_dct)
            elif label_type == 'gwas':
                gwas_dct = read_gwas_file(folder_name)
                dump_label_dct(folder + 'gwas_dct.json', gwas_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))