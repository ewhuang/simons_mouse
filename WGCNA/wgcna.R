### Author: Michael Saul

stopifnot(require("WGCNA"))
stopifnot(require("reshape"))
stopifnot(require("ggplot2"))
stopifnot(require("tools"))

# Setting working directory
setwd("C:/Users/ewhuang3/Documents/simons_mouse")

options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 16)

# hb_mrsb_cpm_t = read.table("mm_mrsb_log2_expression_sampled.tsv", sep="\t")
hb_mrsb_cpm <- read.table("mm_mrsb_log2_expression.tsv", sep="\t", header=T)
gene_id = hb_mrsb_cpm$gene_id
hb_mrsb_cpm$gene_id = NULL
hb_mrsb_cpm_t <- t(hb_mrsb_cpm)
hb_mrsb_cpm_gsg = goodSamplesGenes(hb_mrsb_cpm_t, verbose = 3)

hb_mrsb_wgcna_in = hb_mrsb_cpm_t[, hb_mrsb_cpm_gsg$goodGenes]
gene_id = gene_id[hb_mrsb_cpm_gsg$goodGenes]
hb_mrsb_wgcna_in = hb_mrsb_wgcna_in[, goodSamplesGenes(hb_mrsb_wgcna_in, verbose = 3)$goodGenes]

stopifnot(goodSamplesGenes(hb_mrsb_wgcna_in, verbose = 3)$allOK)

# powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
# hb_mrsb_sft = pickSoftThreshold(hb_mrsb_wgcna_in, powerVector = powers, networkType = "signed", 
#     verbose = 5)

# plot(x = hb_mrsb_sft$fitIndices$Power, hb_mrsb_sft$fitIndices$SFT.R.sq, xlab = "Soft Threshold (power)", 
#     ylab = "R-Squared", type = "l", col = "dark gray", main = "Scale Independence")
# text(hb_mrsb_sft$fitIndices$Power, hb_mrsb_sft$fitIndices$SFT.R.sq, labels = powers, 
#     col = "#3399CC")
# abline(h = 0.85, col = "#FF3333")


# plot(x = hb_mrsb_sft$fitIndices$Power, y = hb_mrsb_sft$fitIndices$mean.k, xlab = "Soft Threshold (power)", 
#     ylab = "Mean Connectivity", main = "Mean Connectivity", col = "dark gray", 
#     type = "l")
# text(hb_mrsb_sft$fitIndices$Power, hb_mrsb_sft$fitIndices$mean.k., labels = powers, 
#     col = "#3399CC")

hb_mrsb_modules = blockwiseModules(hb_mrsb_wgcna_in, power = 8, networkType = "signed", 
    minModuleSize = 30, corType = "pearson", maxBlockSize = 30000, numericLabels = TRUE, 
    saveTOMs = FALSE, verbose = 3)

# Eigenvalues
hb_mrsb_colors = labels2colors(hb_mrsb_modules$colors)
hb_mrsb_eigengenes = hb_mrsb_modules$MEs

hb_mrsb_module_membership = as.data.frame(matrix(ncol = 4, nrow = ncol(hb_mrsb_wgcna_in)))
row.names(hb_mrsb_module_membership) = gene_id # colnames(hb_mrsb_wgcna_in)
colnames(hb_mrsb_module_membership) = c("ID", "module", "color", "module_membership")
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
    hb_mrsb_module_membership[i, "module_membership"] = cor(hb_mrsb_eigengenes[, 
        current_eigengene], hb_mrsb_wgcna_in[, current_id])
}

write.table(hb_mrsb_cpm_gsg$goodGenes, file="good_gene_booleans_WGCNA.txt", sep="\t", row.names=FALSE, col.names=FALSE)

write.table(hb_mrsb_module_membership, file="module_membership_WGCNA.txt", sep="\t", row.names=FALSE)
