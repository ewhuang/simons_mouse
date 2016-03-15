Author: Edward Huang

___________________________CREATING THE RAW NETWORK_____________________________
1. First, create the edges from the raw data. Compute pearson coefficients
between gene expression values to find correlated genes. We can specify a
parameter, pearson_threshold, which determines the cutoff coefficient for
an edge to exist. Output file is raw_network.txt. The second output file,
full_network.txt, contains every single edge weight, including ones that are
low.

$ python gene_edge_weights.py

Output format:
gene_a  gene_b  edge_weight
The edges are not repeated.

1.5 Create the go_edges.txt file from the files Sheng sent in Feb. 29 e-mail.

$ python make_gene_go_file.py.

_________________ADDING GO NODES AND FORMATTING FOR CLUSTERING__________________
2. After generating the gene-gene edges, we create files that will be sent
to the simulated annealing code. Creates 4 files overall, a network and real
network each for a network with and without GO labels. Two parameters to change.
First is percentage of subgraph to randomly sample. Second is the lambda weight
to assign to all gene-GO edges. Output files network_no_go_RUNNUM.txt, where
RUNNUM indicates the run number. The characteristics of each run number can be
found in run_log.txt. real_network_no_go_RUNNUM.txt, network_go_RUNNUM.txt, and
real_network_go_RUNNUM.txt.

$ python create_clustering_input.py RUNNUM lambda subgraph_frac

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

$ python simulated_annealing.py go/no_go temperature num_clusters RUNNUM

Only run the clustering on networks without GO only once, as it will be the
same network for any given percentage of the raw network, since we use a random
seed.

5. Remove GO nodes from cluster files with GO before evaluating.

$ python post_cluster_remove_go.py RUNNUM

Output file: ./results/clusters_go_clean_RUNNUM.txt


6. Runs the Perl script evaluate_clustering.pl to evaluate cluster densities.

$ python evaluate_clustering.py RUNNUM

Outputs cluster evaluation information in ./results/cluster_eval_go/no_go_RUNNUM
Make sure to copy over clusters_no_go.txt if we didn't cluster the network
without GO.

________________________________CLUSTER ANALYSIS________________________________
Compute GO enrichment of each of the clusterings.

7. Compute GO enrichments for each clustering.

$ python compute_go_enrichment.py RUNNUM "predicted"

Needs clusters_go_clean_RUNNUM.txt and clusters_no_go_RUNNUM.txt

Produces a histogram of the two clusterings' p-values. These p-values are all
of the p-values from the top 5 highest correlated GO terms, computed by
Fisher's exact test.
Additionally, outputs a file, called ./results/go_top_go_RUNNUM.txt, which shows
the p-values for each of the clusters to see which clusters have high p-values.
About half of the clusters with GO terms have enrichment values much better than
those of the clusters without GO terms.
Add in the literal string "predicted" without quotes to the end if dealing with
predicted GO edge weights.

8. Analyze the properties of the clusterings.

$ python cluster_info_summary.py RUNNUM

Needs cluster_eval_go/no_go_RUNNUM.txt, output of compute_go_enrichment.py, both
networks, and both raw clusters and networks.

Outputs a file, ./results/clus_info_no_go_RUNNUM.txt
First two lines shows number of genes in input network, number of gene-gene
edges, and number of gene-GO edges. Then, for each cluster, shows the number of
genes in the cluster, number of GO terms, number of gene-gene edges, and number
of gene-GO edges.

9. Perform the wilcoxon rank-sum test on the clusters

$ python wilcoxon_clusters.py RUNNUM

Prints out the score and p-value of the test.

10. Objective function analysis.

This script looks at the output files from analyze_clusters.py. First we look at
clus_inf_no_go_RUNNUM.txt, and then find the median of the in-density for all of
the clusters. We take a threshold of that median, and then look at the lusters 
with GO that meet that threshold. We then count the number of GO terms in these 
clusters, and print that information out.

$ python objective_function.py RUNNUM

No output file, but if the numbers are bigger than 1, then we can say that the
clusters are reasonable.

11. Plotting in-density versus top GO enrichment.

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
indices in Moues_final_Score_matrix.txt, not the actual GO name.
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

___________________________________WGCNA________________________________________
Simply run wgcna.R to cluster on the raw data.
We can run
$ python preprocess_WGCNA.py
to create a file, mm_mrsb_log2_expression_sampled.tsv, which contains only the
genes contained in the randomly sampled network. We change this network with
line 14 in the script. To use this output file to cluster, we change line 14/15
inside wgcna.R

Run the cleaning script to prepare the raw outputs from R for evaluation with
our current python scripts. The script also removes any genes that do not
appear in the sampled network. This can be changed with the block from lines 51
to 60.
$ python clean_WGCNA_module_results.py
Outputs to ./results/WGCNA_results/

To evaluate, we must use a real network from some old network we created.

$ perl evaluate_clustering.pl ./results/WGCNA_results/WGCNA_clusters_all_genes.txt ./data/real_network_no_go_20.txt > ./results/WGCNA_results/WGCNA_cluster_eval.txt

To plot comparisons of p-values between network without GO and WGCNA:

$ python compute_go_enrichment_wgcna.py

