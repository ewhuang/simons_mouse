### Author: Edward Huang

import os
import time

### Splits up the TCGA expression dataset into its cancer sub-categories. Go to
### http://knowcloud.cse.illinois.edu/index.php/s/DSuJffKJR2hIhAa
### Download the clinical_data file and the expr_converted.txt file. Rename
### expr_converted.txt to expr.tsv. Move both files to ./data/tcga_data/
### Run time: 50 seconds.

def read_clinical_data():
    '''
    Gets the clinical data, and maps each disease to the set of corresponding
    sample ID's. Returns a dictionary.
    Key: primary_disease -> str
    Value: sample ID's that have the primary disease -> list(str)
    '''
    disease_to_sample_dct = {}
    f = open('./data/tcga_data/clinical_data', 'r')
    for line in f:
        # Skip header and samples that are not a primary tumor for a patient.
        if 'Primary Tumor' not in line:
            continue

        line = line.split('\t')
        # The primary disease has spaces in the name. Convert to underscores.
        sample_id, primary_disease = line[0], line[19].replace(' ', '_')

        # Update the dictionary.
        if primary_disease not in disease_to_sample_dct:
            disease_to_sample_dct[primary_disease] = []
        disease_to_sample_dct[primary_disease] += [sample_id]        
    f.close()

    return disease_to_sample_dct

def split_expression_data(disease_to_sample_dct):
    '''
    This function splits the TCGA gene expression data into multiple matrices,
    one for each disease with at least 100 samples. Each TCGA type should have
    the same number of rows (genes), but different number of columns (samples).
    Returns a dictionary.
    Key: disease names -> str
    Value: rows of expression matrix -> list(list(str))
    '''
    disease_expr_dct, disease_col_idx_dct = {}, {}

    # Read the TCGA data.
    f = open('./data/tcga_data/expr.tsv', 'r')
    for i, line in enumerate(f):
        line = line.split()

        # Process with the header line.
        if i == 0:
            # Find the column indices that we should take for each disease.
            for disease in disease_to_sample_dct:
                sample_list = disease_to_sample_dct[disease]
                # Get the column indices corresponding to each TCGA disease.
                index_list = [sample_idx for sample_idx, sample in enumerate(
                    line) if (sample in sample_list)]
                if len(index_list) < 100:
                    continue
                disease_col_idx_dct[disease] = index_list
                # Initialize the expression list.
                disease_expr_dct[disease] = [['gene_id'] + [line[j] for j in (
                    index_list)]]
            continue

        # Process all other lines.
        gene, expression_line = line[0], line[1:]
        # Add in the gene expression row for each disease.
        for disease in disease_expr_dct:
            row = [expression_line[idx] for idx in disease_col_idx_dct[disease]]
            # Make the gene the first element in the row.
            disease_expr_dct[disease] += [[gene] + row]
    f.close()
    return disease_expr_dct

def write_expression_files(disease_expr_dct):
    '''
    Write out the gene expression matrices out to file.
    '''
    # Make a new expression file for each disease.
    for disease in disease_expr_dct:
        # Create directory if it does not already exist.
        disease_folder = './data/%s_data' % disease
        if not os.path.exists(disease_folder):
            os.makedirs(disease_folder)

        # Write out the matrix.
        out = open('%s/expr.tsv' % disease_folder, 'w')
        for i, row in enumerate(disease_expr_dct[disease]):
            out.write('%s\n' % ('\t'.join(row)))
        out.close()
    return disease_expr_dct

def generate_directories(disease_expr_dct):
    '''
    Create a directory for each disease with at least 100 samples, if one
    doesn't already exist. Also write out the TCGA disease list.
    '''
    tcga_folder = './data/tcga_data'
    if not os.path.exists(tcga_folder):
        os.makedirs(tcga_folder)

    out = open('%s/tcga_diseases.txt' % tcga_folder, 'w')
    for disease in disease_expr_dct:
        out.write(disease + '\n')
        for method_type in ['go', 'no_go']:
            net_dir = './data/%s_data/networks_%s' % (disease, method_type)
            if not os.path.exists(net_dir):
                os.makedirs(net_dir)
    out.close()

def main():
    disease_to_sample_dct = read_clinical_data()
    disease_expr_dct = split_expression_data(disease_to_sample_dct)
    write_expression_files(disease_expr_dct)
    generate_directories(disease_expr_dct)

if __name__ == '__main__':
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))