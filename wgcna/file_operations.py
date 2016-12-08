def get_tcga_disease_list():
    disease_list = []
    f = open('../data/tcga_data/tcga_diseases.txt', 'r')
    for line in f:
        disease_list += [line.strip()]
    f.close()
    return disease_list

# compute_dbgap_enrichment_wgcna.py
def read_dbgap_file():
    '''
    Gets the DBGAP dictionary. Maps a dbgap ID to a list of genes.
    Key: DBGAP ID -> str
    Value: list of ENSMUSG IDs -> list(str)
    '''
    def get_ensg_to_ensmusg_dct():
        '''
        Gets a dictionary mapping ENSG ID's to their mouse homologs.
        Key: ENSG ID -> str
        Value: list of ENSMUSG IDs -> list(str)
        '''
        ensg_to_ensmusg_dct = {}
        f = open('../data/mouse_data/mart_export.txt', 'r')
        for i, line in enumerate(f):
            # Skip header.
            if i == 0:
                continue
            ensg_id, ensmusg_id = line.split()
            if ensg_id in ensg_to_ensmusg_dct:
                ensg_to_ensmusg_dct[ensg_id] += [ensmusg_id]
            else:
                ensg_to_ensmusg_dct[ensg_id] = [ensmusg_id]
        f.close()
        return ensg_to_ensmusg_dct

    ensg_to_ensmusg_dct = get_ensg_to_ensmusg_dct()

    dbgap_to_ensmusg_dct = {}

    f = open('../data/dbgap.txt', 'r')
    for i, line in enumerate(f):
        dbgap_id, ensg_id, bloat_1, bloat_2 = line.split()
        # Convert human to mouse homolog list.
        if ensg_id not in ensg_to_ensmusg_dct:
            continue
        ensmusg_id_list = ensg_to_ensmusg_dct[ensg_id]

        if dbgap_id in dbgap_to_ensmusg_dct:
            dbgap_to_ensmusg_dct[dbgap_id] += ensmusg_id_list
        else:
            dbgap_to_ensmusg_dct[dbgap_id] = ensmusg_id_list

    f.close()
    return dbgap_to_ensmusg_dct

def read_ensg_dbgap_file():
    '''
    Reads the dbgap dictionary for TCGA data.
    '''
    dbgap_to_ensg_dct = {}
    f = open('../data/dbgap.txt', 'r')
    for i, line in enumerate(f):
        dbgap_id, ensg_id, bloat_1, bloat_2 = line.split()
        if dbgap_id in dbgap_to_ensg_dct:
            dbgap_to_ensg_dct[dbgap_id] += [ensg_id]
        else:
            dbgap_to_ensg_dct[dbgap_id] = [ensg_id]

    f.close()
    return dbgap_to_ensg_dct
