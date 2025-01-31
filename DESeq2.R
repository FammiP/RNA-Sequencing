########################################
# RNA-Seq Analysis Using DESeq2 & GO Enrichment
########################################

# Description:
# This script performs differential expression analysis and gene ontology (GO) enrichment
# for RNA-seq data using DESeq2 and clusterProfiler. It includes exploratory data visualization 
# (PCA plots, volcano plots) and GO term visualizations (dot plots, bar plots, and tree plots).
#
# Inputs:
# - "fCounts_any_end_aligned.txt": Raw count matrix where rows are genes and columns are samples
# - "metadata.csv": Metadata describing experimental conditions for each sample
#
# Outputs:
# - Volcano plots of differentially expressed genes
# - GO enrichment plots (dot plots, bar plots, tree plots)
#
# Usage:
# srun --time=01:00:00 --mem=4G --cpus-per-task=1 --pty /bin/bash
# Rscript DESeq2.R

########################################
# 1. Load necessary libraries and set up environment
########################################

# Set working directory (adjust path as needed)
setwd("/Users/fammiparokkaran/Desktop/RNA-seq")

# Install and load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "biomaRt", "org.Mm.eg.db", "EnhancedVolcano", "clusterProfiler", "GO.db"))

# Load libraries
library(DESeq2)
library(biomaRt)
library(org.Mm.eg.db)
library(EnhancedVolcano)
library(clusterProfiler)
library(enrichplot)
library(GO.db)
library(stringr)
library(ggplot2)

########################################
# 2. Import and preprocess input files
########################################

# Import raw count matrix and metadata
count_raw <- read.table("fCounts_any_end_aligned.txt", header = TRUE, row.names = 1, sep = "\t")

# Remove first 5 unnecessary columns (e.g., annotation information)
counts <- count_raw[, -c(1:5)]

# Extract sample names using SRR IDs
colnames(counts) <- str_extract(colnames(counts), "SRR\\d+")

# Import metadata and create experimental condition columns
metadata <- read.csv("metadata.csv", row.names = 1)
metadata$Tissue <- ifelse(grepl("Lung", metadata$Group), "Lung", "Blood")
metadata$Condition <- ifelse(grepl("Case", metadata$Group), "Case", "Control")

# Match metadata order to the count matrix
metadata <- metadata[match(colnames(counts), rownames(metadata)), ]

# Check for matching samples between counts and metadata
if (!all(colnames(counts) == rownames(metadata))) {
  stop("Sample names in count matrix do not match the metadata. Please check the inputs.")
}

########################################
# 3. Differential Expression Analysis with DESeq2
########################################

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata, design = ~ Tissue + Condition)

# Run DESeq2 analysis
dds <- DESeq(dds)

########################################
# 4. Exploratory Data Analysis (PCA plot)
########################################

vsd <- vst(dds, blind = TRUE)
plotPCA(vsd, intgroup = c("Group")) + ggtitle("PCA of Samples") + theme(plot.title = element_text(hjust = 0.5))

########################################
# 5. Extract DE genes for Blood and Lung Case-Control Comparisons
########################################

# Blood Case vs Blood Control
results_blood <- results(dds, contrast = c("Condition", "Case", "Control"))
DE_blood <- results_blood[which(results_blood$padj < 0.05), ]

# Lung Case vs Lung Control
results_lung <- results(dds, contrast = c("Condition", "Case", "Control"))
DE_lung <- results_lung[which(results_lung$padj < 0.05), ]

########################################
# 6. Annotate DE genes with gene names using biomaRt
########################################

# Connect to Ensembl using biomaRt
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", version = 110)

# Annotate Blood DE genes
gene_info_b <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = rownames(DE_blood),
                     mart = ensembl)
DE_blood_annot <- merge(data.frame(ensembl_gene_id = rownames(DE_blood), DE_blood), gene_info_b, by = "ensembl_gene_id")

# Annotate Lung DE genes
gene_info_l <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "ensembl_gene_id",
                     values = rownames(DE_lung),
                     mart = ensembl)
DE_lung_annot <- merge(data.frame(ensembl_gene_id = rownames(DE_lung), DE_lung), gene_info_l, by = "ensembl_gene_id")

########################################
# 7. Volcano Plots of DE Genes
########################################

EnhancedVolcano(DE_blood_annot, lab = DE_blood_annot$external_gene_name, x = 'log2FoldChange', y = 'pvalue',
                title = "Blood Case vs Blood Control")

EnhancedVolcano(DE_lung_annot, lab = DE_lung_annot$external_gene_name, x = 'log2FoldChange', y = 'pvalue',
                title = "Lung Case vs Lung Control")

########################################
# 8. GO Enrichment Analysis and Visualization
########################################

# Filter valid ENSEMBL IDs for Blood DE genes
valid_blood_genes <- DE_blood_annot$ensembl_gene_id[DE_blood_annot$ensembl_gene_id %in% gene_info_b$ensembl_gene_id]

# GO enrichment for blood
go_blood <- enrichGO(gene = valid_blood_genes, universe = rownames(counts), OrgDb = org.Mm.eg.db, ont = "BP", keyType = "ENSEMBL")

# Blood GO visualization
dotplot(go_blood) + ggtitle("Blood (Case vs Control): GeneRatio")
barplot(go_blood, showCategory = 10) + ggtitle("Blood (Case vs Control): P-value")

# Filter valid ENSEMBL IDs for Lung DE genes
valid_lung_genes <- DE_lung_annot$ensembl_gene_id[DE_lung_annot$ensembl_gene_id %in% gene_info_l$ensembl_gene_id]

# GO enrichment for lung
go_lung <- enrichGO(gene = valid_lung_genes, universe = rownames(counts), OrgDb = org.Mm.eg.db, ont = "BP", keyType = "ENSEMBL")

# Lung GO visualization
dotplot(go_lung) + ggtitle("Lung (Case vs Control): GeneRatio")
barplot(go_lung, showCategory = 10) + ggtitle("Lung (Case vs Control): P-value")

########################################
# 9. Tree Plots of GO Terms (Alternative to goplot)
########################################

# Calculate similarity for treeplot
go_blood <- pairwise_termsim(go_blood)
go_lung <- pairwise_termsim(go_lung)

treeplot(go_blood, showCategory = 10) + ggtitle("GO Tree - Blood Case vs Blood Control")
treeplot(go_lung, showCategory = 10) + ggtitle("GO Tree - Lung Case vs Lung Control")
