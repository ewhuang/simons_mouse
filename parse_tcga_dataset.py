### Author: Edward Huang

import file_operations
import os
import time

### This script splits up the TCGA expression dataset into the sub-cancer
### categories. Go to
### http://knowcloud.cse.illinois.edu/index.php/s/DSuJffKJR2hIhAa
### and download the clinical_data file and the expr_converted.txt file. Rename
### expression file to tcga_expr.txt
### Run time: 70 seconds.

def read_clinical_data():
    '''
    Gets the clinical data, and maps the ID for each sample to its corresponding
    disease. Also returns the count for each disease.
    Returns (dct_1, dct_2)
    dct_1 key: primary_disease -> str
    dct_1 value: sample ID's that have the primary disease -> list(str)
    dct_2 key: primary_disease -> str
    dct_2 value: number of samples for the primary disease -> int
    '''
    disease_to_id_dct, disease_count_dct = {}, {}
    f = open('./data/tcga_data/clinical_data', 'r')
    for i, line in enumerate(f):
        # Skip header.
        if i == 0:
            continue
        line = line.split('\t')
        # Skip samples that are not a primary tumor for a patient.
        if 'Primary Tumor' not in line:
            continue
        # The primary disease has spaces in the name. Convert to underscores.
        sample_id, primary_disease = line[0], line[19].replace(' ', '_')

        # Update the dictionaries.
        if primary_disease in disease_count_dct:
            disease_to_id_dct[primary_disease] += [sample_id]
            disease_count_dct[primary_disease] += 1
        else:
            disease_to_id_dct[primary_disease] = [sample_id]
            disease_count_dct[primary_disease] = 1
    f.close()
    # Remove diseases that do not have at least 100 samples.
    for disease in disease_count_dct:
        if disease_count_dct[disease] < 100:
            del disease_to_id_dct[disease]
    return disease_to_id_dct

def create_disease_directories(disease_to_id_dct):
    '''
    Create a directory for every disease with at least 100 samples, if one
    doesn't already exist.
    '''
    for disease in disease_to_id_dct:
        go_dir = './data/%s_data/networks_go' % disease
        no_go_dir = './data/%s_data/networks_no_go' % disease
        if not os.path.exists(go_dir):
            os.makedirs(go_dir)
        if not os.path.exists(no_go_dir):
            os.makedirs(no_go_dir)

def split_expression_data(disease_to_id_dct):
    '''
    This function splits the TCGA gene expression data into multiple matrices,
    one for each disease with at least 100 samples.
    '''
    # Make a dictionary of lists. Initialize each key, a disease, to an
    # empty list. 
    disease_expr_dct, disease_sample_index_dct = {}, {}
    for disease in disease_to_id_dct:
        disease_expr_dct[disease] = []

    # Read the TCGA data.
    gene_list = []
    f = open('./data/tcga_data/tcga_expr.txt', 'r')
    for i, line in enumerate(f):
        line = line.split()
        if i == 0:
            sample_list = line
            # Find the sample indices we should pull for each disease.
            for disease in disease_to_id_dct:
                index_list = [sample_list.index(sample
                    ) for sample in disease_to_id_dct[disease] if (sample in
                    sample_list)]
                if len(index_list) == 0:
                    del disease_expr_dct[disease]
                    continue
                disease_sample_index_dct[disease] = index_list
            continue
        gene, expression_line = line[0], line[1:]
        assert len(expression_line) == len(sample_list)
        gene_list += [gene]
        for disease in disease_expr_dct:
            disease_sample_expr_lst = [expression_line[idx] for idx in (
                disease_sample_index_dct[disease])]
            disease_expr_dct[disease] += [disease_sample_expr_lst]
    f.close()
    # Make sure every disease matrix is of the same length.
    for disease in disease_expr_dct:
        assert len(disease_expr_dct[disease]) == len(gene_list)
        disease_folder = './data/%s_data/' % disease
        if not os.path.exists(disease_folder):
            os.makedirs(disease_folder)

    tcga_folder = './data/tcga_data/'
    if not os.path.exists(tcga_folder):
        os.makedirs(tcga_folder)
    # Make a new expression file for each disease.
    for disease in disease_expr_dct:
        disease_exp_matrix = disease_expr_dct[disease]
        out = open('./data/%s_data/expr.txt' % disease, 'w')
        out.write('\t%s\n' % '\t'.join([sample_list[i] for i in
            disease_sample_index_dct[disease]]))
        for i, row in enumerate(disease_exp_matrix):
            out.write('%s\t%s\n' % (gene_list[i], '\t'.join(map(str, row))))
        out.close()

    # Write out the list of valid diseases.
    out = open('%stcga_diseases.txt' % tcga_folder, 'w')
    for disease in disease_expr_dct:
        out.write(disease + '\n')
    out.close()

    return disease_expr_dct

def main():
    disease_to_id_dct = read_clinical_data()
    disease_expr_dct = split_expression_data(disease_to_id_dct)
    create_disease_directories(disease_expr_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))