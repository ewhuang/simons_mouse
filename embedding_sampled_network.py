### Author: Edward Huang

import file_operations

### This script generates the edges from embedding in the same format as
### ./data/raw_network.txt.

data_folder = './Sheng/data/network/integrated_network/'

def parse_matrix_gene_ids():
    matrix_gene_ids = []
    f = open(data_folder + 'Mouse.embedding.id', 'r')
    for i, line in enumerate(f):
        ensmusg_id = line.strip()
        assert 'ENSMUSG' in ensmusg_id
        matrix_gene_ids += [ensmusg_id]
    f.close()
    return matrix_gene_ids

def main():
    sampled_genes = file_operations.get_sampled_genes()
    matrix_gene_ids = parse_matrix_gene_ids()
    f = open(data_folder + 'Mouse.embedding', 'r')
    out = open('./data/embedding_network.txt', 'w')
    for i, line in enumerate(f):
        row_gene = matrix_gene_ids[i]
        if row_gene not in sampled_genes:
            continue
        line = line.split()
        for j in range(i + 1, len(line)):
            column_gene = matrix_gene_ids[j]
            if column_gene not in sampled_genes:
                continue
            weight = line[j]
            out.write('%s\t%s\t%s\n' % (row_gene, column_gene, weight))
    out.close()
    f.close()

if __name__ == '__main__':
    main()