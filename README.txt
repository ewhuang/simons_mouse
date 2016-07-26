Author: Edward Huang
Simons Foundation Mouse Project

___________________________CREATING THE GENE NETWORK____________________________
1. Plot the standard deviation distribution of the genes, and write to file.

$ python standard_deviation_hist.py mouse/tcga

2. Makes three json files corresponding to the 3 GO domains. Keys are GO ID's,
values are lists of ENSMUSG ID's. GO terms must not have spaces in the names (
convert to underscores) or else simulated_annealing will make a node for each
word in the term.

$ python dump_go_dictionary_files.py mouse/tcga

3. Compute pearson coefficients between gene expression values to find
correlated genes. We can specify a parameter, pearson_threshold, which
determines the cutoff coefficient for an edge to exist. Output file is
high_std_network.txt.

$ python gene_edge_weights.py mouse/tcga

Output format:
gene_a  gene_b  edge_weight
The edges are not repeated.

_________________ADDING GO NODES AND FORMATTING FOR CLUSTERING__________________
4. Create 4 files overall, a network and real network each for a network with
and without GO labels.
Output files network_go_run_num_FOLD.txt, where run_num indicates the run
number, and FOLD indicates the fold number, as we separate the GO terms into
the three categories: biological process, molecular function, and cellular
component.
Other files:
real_network_go_run_num_FOLD.txt.
network_no_go_run_num_FOLD.txt, and
real_network_no_go_run_num.txt

$ python create_clustering_input.py data_type run_num -b <bootstrap-optional>

Output format for network_go.txt/network_no_go_run_num.txt:
0
NUM_NODES
gene_a  gene_b  edge_weight
gene_b  gene_a  edge_weight
Output format for real_network_no_go_run_num.txt/real_network_go_run_num.txt:
Real network
0   gene_a  gene_b  edge_weight
0   gene_b  gene_a  edge_weight

___________________________________CLUSTERING___________________________________

5. Compile clustering code inside sim_anneal folder.
If static error for EdgeWeightThreshold, add static in front of its declaration
in cs-grn.h.

$ g++ -O3 -o bin/cs-grn -Wno-deprecated -std=c++0x *.cpp

orth.txt just needs to contain at least one gene in the network.
Execute clustering code on the created networks.

4. Run the simulated annealing clustering code.

$ python simulated_annealing.py data_type objective_function run_num go/no_go
            go_num <if go>

Only run the clustering on networks without GO only once, as it will be the
same network for any given percentage of the raw network, since we use a random
seed.

7. Runs the Perl script evaluate_clustering.pl to evaluate cluster densities.

$ python evaluate_clustering.py data_type objective_function run_num

________________________________CLUSTER ANALYSIS________________________________
Compute GO enrichment of each of the clusterings.

8. Compute GO enrichments for each clustering.

$ python compute_go_enrichment.py data_type objective_function run_num

9. Analyze the properties of the clusterings.

$ python cluster_info_summary.py data_type objective_function run_num

10. Perform the wilcoxon rank-sum test on the clusters

$ python wilcoxon_clusters.py run_num

Prints out the score and p-value of the test.

11. Objective function analysis.

This script looks at the output files from analyze_clusters.py. First we look at
clus_inf_no_go_run_num.txt, and then find the median of the in-density for all of
the clusters. We take a threshold of that median, and then look at the lusters 
with GO that meet that threshold. We then count the number of GO terms in these 
clusters, and print that information out.

$ python objective_function.py run_num

No output file, but if the numbers are bigger than 1, then we can say that the
clusters are reasonable.

12. Plotting in-density versus top GO enrichment.

$ python plot_indensity_vs_enrich.py data_type run_num


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

$ python create_clustering_input.py run_num lambda subgraph_decimal "predicted"

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
$ python preprocess_WGCNA.py genes_only/pca/mean/median

2.
Move results from preprocessing to working directory of R.
Run wgcna.R in 64-bit R. Move output (module_membership_WGCNA.txt) to results
file. Change lines 14 and 71 to suit each domain. Domains are bp, cc, and mf.
Takes about 55 minutes.

3.
$ python clean_WGCNA_module_results.py genes_only/pca/mean/median

4.
Pick the real network with which to evaluate (i.e., real_network_no_go_42.txt)
python evaluate_clustering_wgcna genes_only/pca/mean/median

5.
$ python compute_go_enrichment_wgcna.py genes_only/pca/mean/median

6.
$ python cluster_info_summary_WGCNA.py genes_only/pca/mean/median

7.
$ python plot_indensity_vs_enrich_WGCNA.py genes_only/pca/mean/median bp/cc/mf