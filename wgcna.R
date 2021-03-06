### Author: Michael Saul

stopifnot(require("WGCNA"))
stopifnot(require("reshape"))
stopifnot(require("ggplot2"))
stopifnot(require("tools"))

# Setting working directory
setwd("C:/Users/ewhuang3/Documents/simons_mouse")

options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 16)

# Either 'mouse' or a TCGA cancer.
# for (data_type in c('mouse', 'lung_squamous_cell_carcinoma',
#     'colon_adenocarcinoma', 'brain_lower_grade_glioma',
#     'liver_hepatocellular_carcinoma', 'head_&_neck_squamous_cell_carcinoma',
#     'thyroid_carcinoma', 'glioblastoma_multiforme',
#     'uterine_corpus_endometrioid_carcinoma', 'pheochromocytoma_&_paraganglioma',
#     'ovarian_serous_cystadenocarcinoma', 'kidney_papillary_cell_carcinoma',
#     'prostate_adenocarcinoma', 'cervical_&_endocervical_cancer',
#     'pancreatic_adenocarcinoma', 'sarcoma', 'bladder_urothelial_carcinoma',
#     'kidney_clear_cell_carcinoma', 'breast_invasive_carcinoma',
#     'lung_adenocarcinoma')) {
for (data_type in c('tcga')) {

    hb_mrsb_cpm = read.table(paste(data_type, "_expr.tsv", sep=""), sep="\t",
        header=T)
    gene_id = hb_mrsb_cpm$gene_id
    hb_mrsb_cpm$gene_id = NULL
    hb_mrsb_cpm_t <- t(hb_mrsb_cpm)
    hb_mrsb_cpm_gsg = goodSamplesGenes(hb_mrsb_cpm_t, verbose = 3)

    hb_mrsb_wgcna_in = hb_mrsb_cpm_t[, hb_mrsb_cpm_gsg$goodGenes]
    gene_id = gene_id[hb_mrsb_cpm_gsg$goodGenes]
    hb_mrsb_wgcna_in = hb_mrsb_wgcna_in[, goodSamplesGenes(hb_mrsb_wgcna_in,
        verbose = 3)$goodGenes]

    stopifnot(goodSamplesGenes(hb_mrsb_wgcna_in, verbose = 3)$allOK)

    # powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
    # hb_mrsb_sft = pickSoftThreshold(hb_mrsb_wgcna_in, powerVector = powers,
    #     networkType = "signed", verbose = 5)

    # plot(x = hb_mrsb_sft$fitIndices$Power, hb_mrsb_sft$fitIndices$SFT.R.sq,
    #     xlab = "Soft Threshold (power)", ylab = "R-Squared", type = "l",
    #     col = "dark gray", main = "Scale Independence")
    # text(hb_mrsb_sft$fitIndices$Power, hb_mrsb_sft$fitIndices$SFT.R.sq,
    #     labels = powers, col = "#3399mf")
    # abline(h = 0.85, col = "#FF3333")

    # plot(x = hb_mrsb_sft$fitIndices$Power, y = hb_mrsb_sft$fitIndices$mean.k,
    #     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
    #     main = "Mean Connectivity", col = "dark gray", type = "l")
    # text(hb_mrsb_sft$fitIndices$Power, hb_mrsb_sft$fitIndices$mean.k.,
    #     labels = powers, col = "#3399mf")

    hb_mrsb_modules = blockwiseModules(hb_mrsb_wgcna_in, power = 20,
        networkType = "signed", minModuleSize = 30, corType = "pearson",
        maxBlockSize = 30000, numericLabels = TRUE, saveTOMs = FALSE,
        verbose = 3)

    # Eigenvalues
    hb_mrsb_colors = labels2colors(hb_mrsb_modules$colors)
    hb_mrsb_eigengenes = hb_mrsb_modules$MEs

    hb_mrsb_module_membership = as.data.frame(matrix(ncol = 4,
        nrow = ncol(hb_mrsb_wgcna_in)))
    row.names(hb_mrsb_module_membership) = gene_id # colnames(hb_mrsb_wgcna_in)
    colnames(hb_mrsb_module_membership) = c("ID", "module", "color",
        "module_membership")
    colnames(hb_mrsb_wgcna_in) = gene_id
    hb_mrsb_module_membership$ID = gene_id #colnames(hb_mrsb_wgcna_in)
    hb_mrsb_module_membership$module = hb_mrsb_modules$colors
    hb_mrsb_module_membership$color = hb_mrsb_colors
    for (i in 1:nrow(hb_mrsb_module_membership)) {
        # current_id = hb_mrsb_module_membership[i, "ID"]
        # Michael's correction.
        current_id = as.character(hb_mrsb_module_membership[i, "ID"])
        current_module = hb_mrsb_module_membership[i, "module"]
        current_eigengene = paste("ME", current_module, sep = "")
        hb_mrsb_module_membership[i, "module_membership"] = cor(
            hb_mrsb_eigengenes[,current_eigengene], hb_mrsb_wgcna_in[,
            current_id])
    }

    write.table(hb_mrsb_module_membership, file=paste(data_type,
        "_module_membership.txt", sep=""), sep="\t", row.names=FALSE)
}
