### Author: Edward Huang

### This script compares the GO terms between the ones taken from the database
### and the ones used in Sheng's prediction method.

if __name__ == '__main__':
    # Find the GO terms to GO names.
    # Retrieved data from http://geneontology.org/page/download-annotations
    go_id_to_name_dct = {}
    f = open('./go_hierarchy/go_to_name.txt', 'r')
    while True:
        line = f.readline()
        if line == '':
            break
        if line.strip() == '[Term]':
            go_id = f.readline().split()[1]
            go_name = '_'.join(f.readline()[len('name: '):].split())
            go_id_to_name_dct[go_id.lower()] = go_name
            f.readline()
            next = f.readline()
            while 'alt_id' in next:
                go_id_to_name_dct[next.split()[1].lower()] = go_name
                next = f.readline()
    f.close()

    # Keys are the GO ids, values are the indices in the edge weight matrix.
    go_index_dct = {}
    f = open('./go_hierarchy/noisogoHash.txt', 'r')
    for line in f:
        go_id, index = line.split()
        # Subtract 1 to change to list indices.
        go_index_dct[int(index) - 1] = go_id_to_name_dct[go_id]
    f.close()

    predicted_go = set([])
    f = open('./data/network_go_7.txt', 'r')
    for i, line in enumerate(f):
        if i < 2:
            continue
        node_a, node_b, weight = line.split()
        if 'ENSMUSG' not in node_a:
            predicted_go.add(go_index_dct[int(node_a)])
        if 'ENSMUSG' not in node_b:
            predicted_go.add(go_index_dct[int(node_b)])
    f.close()

    database_go = set([])
    f = open('./data/network_go_1.txt', 'r')
    for i, line in enumerate(f):
        if i < 2:
            continue
        node_a, node_b, weight = line.split()
        if 'ENSMUSG' not in node_a:
            database_go.add(node_a)
        if 'ENSMUSG' not in node_b:
            database_go.add(node_b)
    f.close()

    print "Database size: %d" % len(database_go)
    print "Predicted size: %d" % len(predicted_go)
    print "Intersection size: %d" % len(database_go.intersection(predicted_go))