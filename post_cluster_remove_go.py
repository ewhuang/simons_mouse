### Author: Edward Huang


### This file removes the GO nodes from the network_go.txt files
### so that we can evaluate them without the GO nodes inside the
### the clusters.

if __name__ == '__main__':
    f = open('./results/clusters_go_1_001.txt', 'r')
    out = open('./results/clusters_go_clean_1_001.txt', 'w')
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