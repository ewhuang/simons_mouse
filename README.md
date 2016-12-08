# Simons Foundation Mouse Project
Author: Edward Huang

## TCGA Preprocessing

1.  For TCGA data, before anything else, split the TCGA dataset into multiple
    networks, each corresponding to a specific type of cancer.
    
    ```bash
    $ python parse_tcga_dataset.py
    ```

## Downloading GO annotations

1.  Download GO annotations.
    Go to [biomart](http://www.ensembl.org/biomart/martview/). Choose database -> "Ensembl Genes 86" -> Mus musculus genes/Homo sapiens genes

    "Attributes" -> "Gene" -> uncheck "Ensembl Transcript ID" -> "External" -> check "GO Term Name" and "GO domain"

    Hit "Results" at the top and export. Filenames will be mart_export.txt.
    Change mouse annotations to ensmusg_to_go.txt and move to ./data/mouse_data.
    Change human annotations to ensg_to_go.txt and move to ./data/tcga_data.

## Creating the gene network

1.  Plot the standard deviation distribution of the genes, and write to file.

    ```bash
    $ python standard_deviation_hist.py mouse/tcga
    ```

    If input is 'tcga', runs for all TCGA cancer types.

2.  Create GO dictionary JSON files for each gene type. Must run with argument
    mf_go_go to create the the MF GO-GO dictionary.

    ```bash
    $ python dump_go_dictionaries.py mouse/tcga/mf_go_go
    ```

3.  Find overlapping BP and MF terms.

    ```bash
    $ python find_go_overlaps.py mouse/tcga
    ```

4.  Compute pearson coefficients between gene expression values to find
    correlated genes. Output file is high_std_network.txt.

    ```bash
    $ python gene_edge_weights.py mouse/tcga
    ```

## Clustering pipeline

Must run the WGCNA pipeline first in order to compare WlogV to it. This script compiles everything below in this section.

```bash
$ python full_pipeline.py mouse/tcga_index objective_function run_num
```

### Adding GO nodes and formatting for clustering

1.  Create 4 files, a network and real network each for a network with and
    without GO labels.
    Output files network_go_RUNNUM.txt, where RUNNUM indicates the run
    number. For networks with GO labels, we add in the full set of MF terms.
    Other files:
        real_network_go_RUNNUM.txt.
        network_no_go_RUNNUM.txt, and
        real_network_no_go_RUNNUM.txt

    ```bash
    $ python create_clustering_input.py mouse/tcga_cancers run_num -b<bootstrap-optional>
    ```

### Clustering

2.  Compile clustering code inside sim_anneal folder.
    If static error for EdgeWeightThreshold, add static in front of its
    declaration in cs-grn.h.

    ```bash
    $ cd makedir
    $ rm *
    $ cmake ..
    $ make
    ```
    orth.txt just needs to contain at least one gene in the network.

3.  Run the simulated annealing clustering cod to cluster on the networks.

    ```bash
    $ python simulated_annealing.py data_type objective_function run_num go/no_go go_num <if go>
    ```

    Only run the clustering on networks without GO only once, as it will be the
    same network for any given percentage of the raw network, since we use a
    random seed.

4.  Runs the Perl script evaluate_clustering.pl to evaluate cluster densities.
    Illegal division by zero usually means a file doesn't exist.

    ```bash
    $ python evaluate_clustering.py data_type objective_function run_num
    ```

### Cluster analysis

5.  Compute GO enrichments for each clustering.

    ```bash
    $ python compute_go_enrichment.py data_type objective_function run_num
    ```

6.  Compute DBGAP enrichments for each clustering.

    ```bash
    $ python compute_dbgap_enrichment.py data_type objective_function run_num
    ```

7.  Check if genes in clusters labeled by the most enriched BP term in that
    cluster roughly have the same in/(in + out) as genes not labeled by the
    term.

    ```bash
    $ python cheating_evaluation.py data_type objective_function run_num
    ```

8.  Analyze the properties of the clusterings.

    ```bash
    $ python cluster_info_summary.py data_type objective_function run_num
    ```

9.  Plotting in-density versus top GO enrichment.
    
    ```bash
    $ python plot_best_clusters.py data_type run_num
    ```

10. Plotting box plots. One for enrichment, one for in/in + out. 
    Plots for run_num's 1-10.

    ```bash
    $ python box_plot_density_and_enrichment.py data_type clustering_method
    ```

## WGCNA Baseline
1.  data_type can either be 'genes_only' or an integer denoting a TCGA disease.

    ```bash
    cd wgcna/
    $ python preprocess_WGCNA.py data_type genes_only/pca/mean/median
    ```

2.  Move results from preprocessing to working directory of R. Run wgcna.R in
    64-bit R. Move output (module_membership_WGCNA.txt) to results file. Change 
    lines 14 and 71 to suit each domain. Domains are bp, cc, and mf. For TCGA,
    run wgcna_tcga.R. Change line 14 to suit the type of cancer.

    ```bash
    $ python full_pipeline_wgcna.py data_type genes_only/pca/... network_num
    ```

    This runs steps 3-6. network_num should be 1 for mouse, and 1 for TCGA.
    This is because we optimized for mouse, and keep the same parameters for
    TCGA.

3.
    ```bash
    $ python clean_wgcna_module_results.py data_type genes_only/pca/mean/median
    ```

4.  Pick the real network with which to evaluate (i.e., real_network_no_go_42.txt)

    ```bash
    $python evaluate_clustering_wgcna.py data_type genes_only/pca/mean/median network_number
    ```

5.
    ```bash
    $ python compute_go_enrichment_wgcna.py data_type genes_only/pca/mean/median
    ```

6.
    ```bash
    $ python cluster_info_summary_wgcna.py data_type genes_only/pca/mean/median network_number
    ```

7.
    ```bash
    $ python plot_indensity_vs_enrich_WGCNA.py genes_only/pca/mean/median bp/cc/mf
    ```


## Edge weight prediction (deprecated)
Find MGI id to ENSMUSG mappings: http://www.informatics.jax.org/. Find GO id to name mappings: http://geneontology.org/page/download-annotations

1.  Creates an output file, newline separated, where each line is a gene in the
    coexpression matrix. No duplicates. Also creates a sampled list, where each
    gene has to occur in the sampled matrix.

    ```bash
    $ python gene_list.py
    ```

2.  Creates an output file, predicted_go_edge_weights.txt, where each newline is
    a gene and a GO edge, same format as go_edges.txt. However, the GO terms are
    the indices in Mouse_final_Score_matrix.txt, not the actual GO name. Specify
    the number of edges to keep with the variable NUM_TOP_WEIGHTS.

    ```bash
    $ python parse_GO_weight_predictions.py
    ```

3. Simply add an extra keyword, the literal string 'predicted', to create_clustering_input.py to adjust the network to add in the predicted GO edges.

    ```bash
    $ python create_clustering_input.py run_num lambda subgraph_decimal "predicted"
    ```

### Creating matrix for Sheng's Matlab code.

1.  Creates top left of the four-block matrix, which includes the gene-gene edge
    weights. Takes the sampled edges from our 1% network, and computes the
    Pearson correlation coefficient between each pair of genes. The genes are in
    order of of ENSMUSG ID from the original co-expression network.

    ```bash
    $ python top_left.py
    ```

2.  Creates the bottom left and top right blocks of the matrix, or the gene-GO
    edges. 1 if there is a gene-GO relationship, 0 otherwise. GO ordering is
    based on the index from Sheng's data.

    ```bash
    $ python top_right_bottom_left.py
    ```

3.  Outputs a new network, ./data/embedding_edges.txt, which contains
    edges between genes with weights computed by embedding.

    ```bash
    $ python convert_embedding_matrix_to_edge_file.py
    ```

## ProSNet

Must first run create_clustering_input.py, and then let Sheng run.

1.  Run prosnet on the networks created by create_clustering_input.py

    ```bash
    $ python low_dimensional_nodes_prosnet.py mouse/tcga_cancer run_num
    ```

2.  Run k-means on low-dimensional vector representations of genes and GO terms.

    ```bash
    $ python prosnet_kmeans.py mouse/tcga_cancer_index run_num
    ```

3. Evaluate ProSNet.

    ```bash
    $ python evaluate_clustering.py mouse/tcga_cancer_index objective_function run_num
    ```

4. Compute GO enrichment

    ```bash
    $ python compute_go_enrichment.py mouse/tcga_cancer_index objective_function run_num
    ```

5. Summarize results.

    ```bash
    $ python cluster_info_summary.py mouse/tcga_cancer_index objective_function run_num
    ```

## Cluster-One and MCL

1.  Runs the full pipeline after Sheng sends clusters in ./Sheng/Module/
    Does not plot.

    ```bash
    $ python evaluate_co_and_mcl.py mouse/tcga_cancer_index run_num
    ```