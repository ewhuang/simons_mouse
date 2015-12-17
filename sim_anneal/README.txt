## InOut Ratio Clustering Model

1. Cost function: \sum_cluster_i (In-Density(cluster_i) / (In-Density(cluster_i) + Out-Density(cluster_i))) * size(cluster_i)
    a. In-Density = In-cluster edges / In-cluster nodes
    b. Out-Density = Outgoing-cluster edges / Out-cluster nodes
    c. size = number of nodes

2. This program has written code (currently commented) to support multithreading for many trials. But this option runs slower than it would be in a single thread, one possible reason could be "maxed out your memory bandwidth". Another issue is it has not been successfully compiled on HAL

3. This program allows a "trashcan" cluster to hold "non-clusterd" nodes, where no delta-cost = 0 if nodes are in this cluster. The implementation detail is as followed:
    a. increase cluster num += 1
    b. whenever nodes go into cluster0, no cost is added from cluster0, and delta cost = 0

4. Runs in synthetic data reveals that cluster0 could be empty, and therefore the clustering result should be similar to the one that no "trashcan" cluster is allowed.

5. A quick-and-loose analytical result shows the "trashcan" cluster will not affect the best cost:
    a. A network: 500 nodes, 10 clusters
    b. y does not change with the increase of x (in expected value)
    c. cost : y = (245 - 4.9x) / ( (245 - 4.9x) / (50-x) + (1125-22.5x) / 450) + 9 * ( 245 / ( (245 / 50) + (1125 - 2.5x) / (450 - x) ) )