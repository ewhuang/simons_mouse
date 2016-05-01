Author: Edward Huang
Simons Foundation Mouse Project

___________________________CREATING THE GENE NETWORK____________________________
1. Plot the standard deviation distribution of the genes, and write to file.

$ python standard_deviation_hist.py

2. Makes three json files corresponding to the 3 GO domains. Keys are GO ID's,
values are lists of ENSMUSG ID's.

$ python dump_go_dictionary_files.py.

3. Compute pearson coefficients between gene expression values to find
correlated genes. We can specify a parameter, pearson_threshold, which
determines the cutoff coefficient for an edge to exist. Output file is
high_std_network.txt. 

$ python gene_edge_weights.py

Output format:
gene_a  gene_b  edge_weight
The edges are not repeated.

_________________ADDING GO NODES AND FORMATTING FOR CLUSTERING__________________
4. Create 4 files overall, a network and real network each for a network with
and without GO labels.
Output files network_go_RUNNUM_FOLD.txt, where RUNNUM indicates the run
number, and FOLD indicates the fold number, as we separate the GO terms into
the three categories: biological process, molecular function, and cellular
component.
Other files:
real_network_go_RUNNUM_FOLD.txt.
network_no_go_RUNNUM_FOLD.txt, and
real_network_no_go_RUNNUM.txt

$ python create_clustering_input.py RUNNUM

Output format for network_go.txt/network_no_go_RUNNUM.txt:
0
NUM_NODES
gene_a  gene_b  edge_weight
gene_b  gene_a  edge_weight
Output format for real_network_no_go_RUNNUM.txt/real_network_go_RUNNUM.txt:
Real network
0   gene_a  gene_b  edge_weight
0   gene_b  gene_a  edge_weight

___________________________________CLUSTERING___________________________________

3. Compile clustering code inside sim_anneal folder.
If static error for EdgeWeightThreshold, add static in front of its declaration
in cs-grn.

$ g++ -O3 -o bin/cs-grn -Wno-deprecated -std=c++0x *.cpp

orth.txt just needs to contain at least one gene in the network.
Execute clustering code on the created networks.

4. Run the simulated annealing clustering code.

$ python simulated_annealing.py RUNNUM go/no_go

Only run the clustering on networks without GO only once, as it will be the
same network for any given percentage of the raw network, since we use a random
seed.

7. Runs the Perl script evaluate_clustering.pl to evaluate cluster densities.

$ python evaluate_clustering.py RUNNUM

Outputs cluster evaluation information in ./results/cluster_eval_go/no_go_RUNNUM
Make sure to copy over clusters_no_go.txt if we didn't cluster the network
without GO.

________________________________CLUSTER ANALYSIS________________________________
Compute GO enrichment of each of the clusterings.

8. Compute GO enrichments for each clustering.

$ python compute_go_enrichment.py RUNNUM

9. Analyze the properties of the clusterings.

$ python cluster_info_summary.py RUNNUM

10. Perform the wilcoxon rank-sum test on the clusters

$ python wilcoxon_clusters.py RUNNUM

Prints out the score and p-value of the test.

11. Objective function analysis.

This script looks at the output files from analyze_clusters.py. First we look at
clus_inf_no_go_RUNNUM.txt, and then find the median of the in-density for all of
the clusters. We take a threshold of that median, and then look at the lusters 
with GO that meet that threshold. We then count the number of GO terms in these 
clusters, and print that information out.

$ python objective_function.py RUNNUM

No output file, but if the numbers are bigger than 1, then we can say that the
clusters are reasonable.

12. Plotting in-density versus top GO enrichment.

$ python plot_indensity_vs_enrich.py RUN_NUM


_____________________WORKING WITH GO EDGE WEIGHT PREDICTION_____________________
Find MGI id to ENSMUSG mappings: http://www.informatics.jax.org/
Find GO id to name mappings: http://geneontology.org/page/download-annotations

1.
$ python gene_list.py

Creates an output file, newline separated, where each line is a gene in the
coexpression matrix. No duplicates. Also creates a sampled list, where each
gene has to occur in the sampled matrix.

2.
$ python parse_GO_weight_predictions.py

Creates an output file, predicted_go_edge_weights.txt, where each newline is
a gene and a GO edge, same format as go_edges.txt. However, the GO terms are the
indices in Mouse_final_Score_matrix.txt, not the actual GO name.
Specify the number of edges to keep with the variable NUM_TOP_WEIGHTS.

3.
Simply add an extra keyword, the literal string 'predicted', to 
create_clustering_input.py to adjust the network to add in the predicted GO
edges.

$ python create_clustering_input.py RUNNUM lambda subgraph_decimal "predicted"

CREATING MATRIX FOR SHENG'S MATLAB CODE.
1. $ python top_left.py
Creates top left of the four-block matrix, which includes the gene-gene edge
weights. Takes the sampled edges from our 1% network, and computes the Pearson
correlation coefficient between each pair of genes. The genes are in order of
of ENSMUSG ID from the original co-expression network.

2. $ python top_right_bottom_left.py
Creates the bottom left and top right blocks of the matrix, or the gene-GO
edges. 1 if there is a gene-GO relationship, 0 otherwise. GO ordering is based
on the index from Sheng's data.

3. $ python convert_embedding_matrix_to_edge_file.py
Outputs a new network, ./data/embedding_edges.txt, which contains
edges between genes with weights computed by embedding.

___________________________________WGCNA________________________________________
1.
cd wgcna/
$ python preprocess_WGCNA.py

2.
Run wgcna.R in R. Move output (module_membership_WGCNA.txt) to results file.

3.
$ python clean_WGCNA_module_results.py

4.
$ perl ./sim_anneal/evaluate_clustering.pl ./wgcna/results/WGCNA_clusters_high_std_genes.txt ./data/networks_no_go/real_network_no_go_42.txt > ./wgcna/results/WGCNA_cluster_eval.txt

5.
$ python compute_go_enrichment_wgcna.py

6.
$ python cluster_info_summary_WGCNA.py