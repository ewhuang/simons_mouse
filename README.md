# Simons Foundation Mouse Project
Author: Edward Huang

## TCGA Preprocessing

1.  Split the TCGA dataset into multiple networks, each corresponding to a
    specific type of cancer. This should be run before anything else.
    
    ```bash
    $ python split_tcga_dataset.py
    ```

## Downloading annotations

1.  Download GO annotations.
    Go to [biomart](http://www.ensembl.org/biomart/martview/). Choose database -> "Ensembl Genes 87" -> Mus musculus genes/Homo sapiens genes

    "Attributes" -> "Gene" -> uncheck "Ensembl Transcript ID" -> "External" -> check "GO Term Name" and "GO domain"

    Hit "Results" at the top and export to file as TSV. Filenames will be
    mart_export.txt. Remember to tick "unique results only".
    Change mouse annotations to ensmusg_to_go.txt and move to ./data/mouse_data. 
    Change human annotations to ensg_to_go.txt and move to ./data/tcga_data.

2.  Download DBGAP annotations.
    Go to http://veda.cs.uiuc.edu/TCGA_classify/msigdb/gene_sets/dbgap_all/,
    download dbgap.edge, and rename to dbgap.txt. Move to ./data/.
    To get the translation files, obtain mart_export.txt by going to
    ensembl.org/biomart:
    Dataset -> Homo sapiens genes (GRCh38.p7)
    Filters -> Multi Species Comparisons -> Orthologous Mouse Genes: Only
    Attributes -> Ensembl Gene ID, uncheck transcript ID
    Add another dataset, [Ensembl genes 87] Mouse genes, then results.
    Export as TSV, tick "unique results only". Move to ./data/mouse_data/,
    rename as ensg_to_ensmusg.txt.

3.  Download DisGeNET (GWAS) annotations.
    Go to http://www.disgenet.org/web/DisGeNET/menu/downloads and download
    the curated gene-disease associations. Move to ./data/.
    Go to biomart (as in GO and DBGAP).
    Dataset -> Homo sapiens genes
    Attributes -> Gene ID, uncheck transcript, EntrezGene ID
    Export, rename to entrez_to_ensg.txt Move to ./data/tcga_data/. Only get
    human genes because for mouse, we translate from Entrez->ENSG->ENSMUSG. The
    direct Entrez-ENSMUSG database is quite sparse.

4.  Download CTD gene-pathway associations.
    Go to http://ctdbase.org/downloads/, and download the tsv files for different gene associations. Gene-pathway associations are under CTD_genes_pathways.tsv. Move it to ./data/.

## Creating the gene network

1.  Plot the standard deviation distribution of the genes, and write to file.
    Only keeps genes with standard deviation > 0.1.

    ```bash
    $ python standard_deviation_hist.py mouse/tcga
    ```

    If command line argument is 'tcga', runs for all TCGA cancer types.

2.  Create GO dictionary JSON files for each gene type. Must run with argument
    mf_go_go to create the the MF GO-GO dictionary.

    ```bash
    $ python dump_label_dictionaries.py mouse/tcga/mf_go_go go/dbgap/gwas
    ```
    Last argument doesn't matter if argument is mf_go_go

<!-- 3.  Find overlapping BP and MF terms.

    ```bash
    $ python find_go_overlaps.py mouse/tcga
    ``` -->

4.  Compute Pearson coefficients between gene expression values to find
    correlated genes. Output file is high_std_network.txt.

    ```bash
    $ python gene_edge_weights.py mouse/tcga
    ```

## WGCNA Pre-processing
1.  Must have previously run split_tcga_dataset.py and standard_deviation_hist.py.

    ```bash
    $ python preprocess_wgcna.py mouse/tcga
    ```

2.  Move results from preprocessing to working directory of R. Run wgcna.R in
    64-bit R (you can just copy paste the contents into the R shell). Move
    output (%s_module_membership.txt) to ./data/wgcna_data. Takes roughly 45
    minutes per dataset.

## Clustering pipeline

This script compiles everything below in this section.
Note: must run create_clustering_input.py prior to running evaluate_clustering
for WGCNA. This is so WGCNA has the "true" network to evaluate in/out ratio.
Must have run ./wgcna/wgcna.R prior to running simulated_annealing.py. This is
so we know how many clusters to use.
Must have run everything for WGCNA prior to plotting. This is so we have
something to plot for WGCNA.

```bash
$ python full_pipeline.py mouse/tcga_index wlogv/wgcna run_num
```

### Adding GO nodes and formatting for clustering

1.  Create 4 files, a network and real network each for a network with and
    without GO labels.
    Output files network_go_RUNNUM.txt, where RUNNUM indicates the run
    number. For networks with GO labels, we add in the full set of MF terms.
    Other files:
        real_network_go_RUNNUM.txt.
        network_no_go_RUNNUM.txt
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
    $ python compute_label_enrichments.py data_type objective_function run_num go/dbgap/gwas
    ```

6.  Check if genes in clusters labeled by the most enriched BP term in that
    cluster roughly have the same in/(in + out) as genes not labeled by the
    term.

    ```bash
    $ python cheating_evaluation.py data_type objective_function run_num
    ```

7.  Analyze the properties of the clusterings.

    ```bash
    $ python cluster_info_summary.py data_type objective_function run_num
    ```

8.  Plotting in-density versus top GO enrichment.

    ```bash
    $ python plot_best_clusters.py data_type run_num go/go_auc/dbgap
    ```

9.  Plotting box plots. One for enrichment, one for in/in + out. 
    Plots for run_num's 1-10.

    ```bash
    $ python box_plot_density_and_enrichment.py data_type clustering_method
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