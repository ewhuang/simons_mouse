Author: Edward Huang
__________________________________
First, create the edges from the raw data. Compute pearson coefficients
between gene expression values to find correlated genes. We can specify a
parameter, pearson_threshold, which determines the cutoff coefficient for
an edge to exist. Output file is raw_network.txt. 
>>> python gene_edge_weights.py
Output format:
gene_a  gene_b  edge_weight
The edges are not repeated.

_________________________________
After generating the gene-gene edges, we create files that will be sent
to the simulated annealing code. Creates 4 files overall, a network and real
network each for a network with and without GO labels. Two parameters to change.
First is percentage of subgraph to randomly sample. Second is the lambda weight
to assign to all gene-GO edges. Output files network_no_go_RUNNUM.txt, where RUNNUM
indicates the run number. The characteristics of each run number can be found in
run_log.txt. real_network_no_go_RUNNUM.txt, network_go_RUNNUM.txt, and
real_network_go_RUNNUM.txt.
>>> python create_clustering_input.py
Output format for network_go.txt/network_no_go_RUNNUM.txt:
0
NUM_NODES
gene_a  gene_b  edge_weight
gene_b  gene_a  edge_weight
Output format for real_network_no_go_RUNNUM.txt/real_network_go_RUNNUM.txt:
Real network
0   gene_a  gene_b  edge_weight
0   gene_b  gene_a  edge_weight

___________________________________
Compile clustering code inside sim_anneal folder.
>>> g++ -O3 -o bin/cs-grn -Wno-deprecated -std=c++0x *.cpp

orth.txt just needs to contain at least one gene in the network.
Execute clustering code on the created networks.
>>> ./sim_anneal/bin/cs-grn 20 1 0 ./data/orth.txt 1 ./data/network_go_RUNNUM.txt -t 1 2> log > ./results/clusters_go_RUNNUM.txt
>>> ./sim_anneal/bin/cs-grn 20 1 0 ./data/orth.txt 1 ./data/network_no_go_RUNNUM.txt -t 1 2> log > ./results/clusters_no_go_RUNNUM.txt

Remove GO nodes from cluster files before evaluating.
>>> python post_cluster_remove_go.py
Creates ./results/clusters_go_clean_RUNNUM.txt

Evaluate clusters.
>>> perl evaluate_clustering.pl ./results/clusters_no_go_RUNNUM.txt ./data/real_network_no_go_RUNNUM.txt > ./results/cluster_eval_no_go_RUNNUM.txt

>>> perl evaluate_clustering.pl ./results/clusters_go_clean_RUNNUM.txt ./data/real_network_go_RUNNUM.txt > ./results/cluster_eval_go_RUNNUM.txt

__________________________________
Compute GO enrichment of each of the clusterings.
>>> python compute_go_enrichment.py
Produces a histogram of the two clusterings' p-values. These p-values are all
of the p-values from the top 5 highest correlated GO terms, computed by
Fisher's exact test.
Additionally, outputs a file, called ./results/go_top_go_RUNNUM.txt, which shows
the p-values for each of the clusters to see which clusters have high p-values.
About half of the clusters with GO terms have enrichment values much better than
those of the clusters without GO terms.

___________________________________
Analyze the properties of the clusterings.
>>> python analyze_clusters.py RUNNUM
Outputs a file, ./results/clus_info_no_go_RUNNUM.txt
First two lines shows number of genes in input network, number of gene-gene
edges, and number of gene-GO edges. Then, for each cluster, shows the number of
genes in the cluster, number of GO terms, number of gene-gene edges, and number
of gene-GO edges.

Perform the wilcoxon rank-sum test on the cluster in-densities and out-
densities.
>>> python wilcoxon_clusters.py RUNNUM
Prints out the score and p-value of the test.

This script looks at the output files from analyze_clusters.py. First we look at 
clus_inf_no_go_RUNNUM.txt, and then find the median of the in-density for all of
the clusters. We take a threshold of that median, and then look at the lusters 
with GO that meet that threshold. We then count the number of GO terms in these 
clusters, and print that information out.
>>> python objective_function.py RUNNUM
No output file, but if the numbers are bigger than 1, then we can say that the
clusters are reasonable.