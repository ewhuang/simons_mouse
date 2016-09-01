Author: Edward Huang
Simons Foundation Mouse Project

For TCGA data, before anything else, run
$ python parse_tcga_dataset.py
This splits the TCGA dataset into multiple networks, each corresponding to a
specific type of cancer.

______________________________PARSING GO TERMS__________________________________
1. Download mouse gene association zip file from
http://geneontology.org/page/download-annotations

Remove the first few description lines from the file.

2. Extract genes from the file as a list of MGI ID's.
$ python extract_genes_from_gene_associations.py mouse/tcga

3. Upload the gene list to
http://www.informatics.jax.org/batch
Upload the file, select MGI Gene/Marker ID for input type, check only the
Ensembl ID box, and click search. In the resulting window, click Export: 
Text File. Rename the resulting files mgi_to_ensembl.txt

Also download geneontology.org/ontology/go-basic.obo
Move it to ./data/

4. As it is right now, only creates the MF GO-GO dictionary.
$ python dump_mgi_go_dictionaries.py mouse/tcga

5. Find overlapping BP and MF terms.
$ python find_go_overlaps.py mouse

___________________________CREATING THE GENE NETWORK____________________________
1. Plot the standard deviation distribution of the genes, and write to file.

$ python standard_deviation_hist.py mouse/tcga

If input is 'tcga', runs for all TCGA cancer types.

2. Makes three json files corresponding to the 3 GO domains. Keys are GO ID's,
values are lists of ENSMUSG ID's. GO terms must not have spaces in the names (
convert to underscores) or else simulated_annealing will make a node for each
word in the term. Run this for each TCGA cancer.

$ python dump_go_dictionary_files.py mouse/tcga_cancer_index

3.
Find GO terms from BP and MF that overlap with each other.
$ python find_go_overlaps.py mouse/tcga

4. Create the dictionary mapping GO terms to neighboring GO terms.
$ python make_go_go_dictionary.py mouse

5. Compute pearson coefficients between gene expression values to find
correlated genes. We can specify a parameter, pearson_threshold, which
determines the cutoff coefficient for an edge to exist. Output file is
high_std_network.txt.

$ python gene_edge_weights.py mouse/tcga_cancer_index

Output format:
gene_a  gene_b  edge_weight
The edges are not repeated.

_______________________________FULL PIPELINE____________________________________
We can run everything at once and ignore steps 4-12. Must run the WGCNA pipeline
first in order to compare WlogV to it.
python full_pipeline.py mouse/tcga_cancer_index objective_function run_num

_________________ADDING GO NODES AND FORMATTING FOR CLUSTERING__________________
5. Create 4 files overall, a network and real network each for a network with
and without GO labels.
Output files network_go_run_num_FOLD.txt, where run_num indicates the run
number, and FOLD indicates the fold number, as we separate the GO terms into
the three categories: biological process, molecular function, and cellular
component.
Other files:
real_network_go_run_num_FOLD.txt.
network_no_go_run_num_FOLD.txt, and
real_network_no_go_run_num.txt

$ python create_clustering_input.py mouse/tcga_cancers run_num
                                        -b<bootstrap-optional>

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

6. Compile clustering code inside sim_anneal folder.
If static error for EdgeWeightThreshold, add static in front of its declaration
in cs-grn.h.

cd makedir
rm *
cmake ..
make

orth.txt just needs to contain at least one gene in the network.
Execute clustering code on the created networks.

7. Run the simulated annealing clustering code.

$ python simulated_annealing.py data_type objective_function run_num go/no_go
            go_num <if go>

Only run the clustering on networks without GO only once, as it will be the
same network for any given percentage of the raw network, since we use a random
seed.

8. Runs the Perl script evaluate_clustering.pl to evaluate cluster densities.
Illegal division by zero usually means a file doesn't exist.

$ python evaluate_clustering.py data_type objective_function run_num

________________________________CLUSTER ANALYSIS________________________________
Compute GO enrichment of each of the clusterings.

9. Compute GO enrichments for each clustering.

$ python compute_go_enrichment.py data_type objective_function run_num

10. Analyze the properties of the clusterings.

$ python cluster_info_summary.py data_type objective_function run_num

11. Plotting in-density versus top GO enrichment.

$ python plot_best_clusters.py data_type run_num


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

__________________________________PROSNET_______________________________________
Must first run create_clustering_input.py, and then let Sheng run.

1. Run k-means on low-dimensional vector representations of genes and GO terms.
$ python prosnet_kmeans.py mouse/tcga_cancer_index run_num

2. Evaluate PROSNET.
$ python evaluate_clustering.py mouse/tcga_cancer_index objective_function
                                    run_num

3. Compute GO enrichment
$ python compute_go_enrichment.py mouse/tcga_cancer_index objective_function
                                    run_num

4. Summarize results.
$ python cluster_info_summary.py mouse/tcga_cancer_index objective_function
                                    run_num

___________________________________WGCNA________________________________________
1.
data_type can either be 'genes_only' or an integer denoting a TCGA disease.

cd wgcna/
$ python preprocess_WGCNA.py data_type genes_only/pca/mean/median

2.
Move results from preprocessing to working directory of R.
Run wgcna.R in 64-bit R. Move output (module_membership_WGCNA.txt) to results
file. Change lines 14 and 71 to suit each domain. Domains are bp, cc, and mf.
Takes about 55 minutes.
For TCGA, run wgcna_tcga.R. Change line 14 to suit the type of cancer.

____$ python full_pipeline_wgcna.py data_type genes_only/pca/... network_num
____This runs steps 3-6. network_num should be 20 for mouse, and 1 for TCGA.
____This is because we optimized for mouse, and keep the same parameters for
____TCGA.

3.
$ python clean_wgcna_module_results.py data_type genes_only/pca/mean/median

4.
Pick the real network with which to evaluate (i.e., real_network_no_go_42.txt)
python evaluate_clustering_wgcna.py data_type genes_only/pca/mean/median network_number

5.
$ python compute_go_enrichment_wgcna.py data_type genes_only/pca/mean/median

6.
$ python cluster_info_summary_wgcna.py data_type genes_only/pca/mean/median network_number

7.
$ python plot_indensity_vs_enrich_WGCNA.py genes_only/pca/mean/median bp/cc/mf

old stuff
___
Perform the wilcoxon rank-sum test on the clusters

$ python wilcoxon_clusters.py run_num

Prints out the score and p-value of the test.

___
Objective function analysis.

This script looks at the output files from analyze_clusters.py. First we look at
clus_inf_no_go_run_num.txt, and then find the median of the in-density for all of
the clusters. We take a threshold of that median, and then look at the lusters 
with GO that meet that threshold. We then count the number of GO terms in these 
clusters, and print that information out.

$ python objective_function.py run_num

No output file, but if the numbers are bigger than 1, then we can say that the
clusters are reasonable.