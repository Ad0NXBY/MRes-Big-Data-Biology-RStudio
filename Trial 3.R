#Installation of necessary packages
install.packages("ggplot2") # for plotting
install.packages("tidyverse")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("biomaRt")

#Load necessary libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(EnhancedVolcano)
library(biomaRt)

#Define file path
files <- list.files("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts",
                    pattern = "*.txt", full.names = TRUE)

#Read each file into a list of dataframes
counts_list <- lapply(files, function(file) {
  df <- read.table(file, header = TRUE)
  df <- df %>% select(Geneid, ends_with(".bam"))
  colnames(df) <- c("Geneid", basename(file))
  return(df)
})

#Merge all dataframes by "Geneid"
merged_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), counts_list)

# Map Ensembl gene IDs to gene names using biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = merged_counts$Geneid,
                      mart = mart)

# Merge gene names with the count data
merged_counts_with_names <- merge(merged_counts, gene_mapping, by.x = "Geneid", by.y = "ensembl_gene_id", all.x = TRUE)

# Add gene names as row names for count data
rownames(merged_counts_with_names) <- merged_counts_with_names$external_gene_name

# Keep only the count columns, excluding Geneid and gene name
count_data <- merged_counts_with_names[, -c(1, ncol(merged_counts_with_names))]

# Write the merged dataframe to a new file if needed
write.table(merged_counts_with_names, file = "merged_featurecounts_with_gene_names.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Create metadata dataframe from experimental design
# plenti(wildtype) is A1-A3, KO22 is B1-B3, KO23 is C1-C3
condition <- factor(c(rep("plenti", 3), rep("KO22", 3), rep("KO23", 3)))

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = data.frame(condition),
                              design = ~condition)

# DESeq2 analysis
dds <- DESeq(dds)

# Transform the data for PCA
vsd <- vst(dds, blind = FALSE) #variance stabilizing transformation

# Generate a PCA plot
plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA Plot of RNA-Seq Data")

# DESeq2 results for comparisons
# plenti vs KO22
results_plenti_vs_KO22 <- results(dds, contrast = c("condition", "KO22", "plenti"))
summary(results_plenti_vs_KO22)

# plenti vs KO23
results_plenti_vs_KO23 <- results(dds, contrast = c("condition", "KO23", "plenti"))
summary(results_plenti_vs_KO23)

# Generate Volcano Plots
# Volcano plot of plenti vs KO22
EnhancedVolcano(results_plenti_vs_KO22,
                lab = rownames(results_plenti_vs_KO22),
                x = "log2FoldChange",
                y = "pvalue",
                title = "plenti vs KO22",
                pCutoff = 0.05,
                FCcutoff = 2.0,
                labSize = 3,
                selectLab = rownames(results_plenti_vs_KO22)[which(abs(results_plenti_vs_KO22$log2FoldChange) > 2 & results_plenti_vs_KO22$pvalue < 0.01)])

# Volcano plot of plenti vs KO23
EnhancedVolcano(results_plenti_vs_KO23,
                lab = rownames(results_plenti_vs_KO23),
                x = "log2FoldChange",
                y = "pvalue",
                title = "plenti vs KO23",
                pCutoff = 0.05,
                FCcutoff = 1.5)
