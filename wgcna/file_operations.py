def get_tcga_disease_list():
    disease_list = []
    f = open('../data/tcga_data/tcga_diseases.txt', 'r')
    for line in f:
        disease_list += [line.strip()]
    f.close()
    return disease_list