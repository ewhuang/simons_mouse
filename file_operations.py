### Author: Edward Huang

from collections import OrderedDict
import math

### This file contains functions that parse the data files and return the 
### data objects that we work with in our scripts.

# This function returns a dictionary, with keys as the names of GO annotations
# and values as lists of genes annotated by the keys.
def get_go_labels():
    go_dct = {}
    f = open('./data/go_edges.txt', 'r')
    for line in f:
        gene, go_label = line.split()
        if go_label not in go_dct:
            go_dct[go_label] = [gene]
        else:
            go_dct[go_label] += [gene]
    f.close()
    return go_dct

def create_clean_go_file(run_num):
    f = open('./results/clusters_go_%s.txt' % run_num, 'r')
    out = open('./results/clusters_go_clean_%s.txt' % run_num, 'w')
    for i, line in enumerate(f):
        if i == 0:
            out.write(line)
            continue
        newline = line.split()
        if 'ENSMUSG' not in newline[3]:
            continue
        out.write(line)
    out.close()
    f.close()