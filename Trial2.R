# Installation of necessary packages (if not already installed)
#Install packages as needed (some of these may already be installed)
install.packages("ggplot2")
install.packages("tidyverse")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Bioconductor installations
BiocManager::install(c("GenomicFeatures", "DESeq2", "EnhancedVolcano", "biomaRt"), force = TRUE)

# Load necessary libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(EnhancedVolcano)



# Define file path and read feature count files
files <- list.files("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts",
                    pattern = "*.txt", full.names = TRUE)

# Read each file into a list of dataframes
counts_list <- lapply(files, function(file) {
  df <- read.table(file, header = TRUE)
  df <- df %>% select(Geneid, ends_with(".bam"))
  colnames(df) <- c("Geneid", basename(file))
  return(df)
})

# Merge all dataframes by "Geneid"
merged_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), counts_list)


# Fetch gene names using biomaRt
library(biomaRt)
# Use Ensembl to map Gene IDs (e.g., ENSG000001) to gene names
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = merged_counts$Geneid,
                      mart = mart)

# Merge gene names with the count data
merged_counts <- merge(merged_counts, gene_mapping, by.x = "Geneid", by.y = "ensembl_gene_id", all.x = TRUE)

# Replace rownames with gene names (if available), otherwise keep Ensembl IDs
rownames(merged_counts) <- ifelse(!is.na(merged_counts$external_gene_name), 
                                  merged_counts$external_gene_name, 
                                  merged_counts$Geneid)

# Drop unnecessary columns (like Geneid and external_gene_name)
count_data <- merged_counts[, -c(1, ncol(merged_counts))]

# Create metadata dataframe from experimental design
# plenti (wildtype) is A1-A3, KO22 is B1-B3, KO23 is C1-C3
condition <- factor(c(rep("plenti", 3), rep("KO22", 3), rep("KO23", 3)))

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = data.frame(condition),
                              design = ~condition)

# DESeq2 analysis
dds <- DESeq(dds)

# Transform the data for PCA
vsd <- vst(dds, blind = FALSE) # variance stabilizing transformation

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
                lab = rownames(merged_counts),
                x = "log2FoldChange",
                y = "pvalue",
                title = "plenti vs KO22",
                pCutoff = 0.05,
                FCcutoff = 2.0,
                labSize = 3,
                selectLab = rownames(results_plenti_vs_KO22)[which(abs(results_plenti_vs_KO22$log2FoldChange) > 2 & results_plenti_vs_KO22$pvalue < 0.01)])

# Volcano plot of plenti vs KO23
EnhancedVolcano(results_plenti_vs_KO23,
                lab = rownames(merged_counts),
                x = "log2FoldChange",
                y = "pvalue",
                title = "plenti vs KO23",
                pCutoff = 0.05,
                FCcutoff = 1.5)

