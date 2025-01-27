#PREPROCESSING============================================

#Installation of ggplot2
install.packages("ggplot2") # for plotting

#Installation of tidyverse
install.packages("tidyverse")

#Installation of BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

#Installation of Genomic Features
BiocManager::install(c("GenomicFeatures"))

#Installation of DESeq2
BiocManager::install("DESeq2")

#Installation of EnhancedVolcano
BiocManager:: install("EnhancedVolcano")

#Load necessary libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(tidyverse)
library(EnhancedVolcano)

#define file path
files <- list.files("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts",
                    pattern = "*.txt", full.names = TRUE)

#read each file into a list of dataframes
counts_list <- lapply(files, function(file) {
  df <- read.table(file, header = TRUE)
  df <- df %>% select(Geneid, ends_with(".bam"))
  colnames(df) <- c("Geneid", basename(file))
  return(df)
})

View(counts_list[[1]])

#Extract Gene list
Genelist <- as.data.frame(counts_list[[1]][["Geneid"]])
colnames(Genelist) <- "Gene"

#Merge all dataframes by "Geneid"
merged_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), counts_list)

#merged counts loses the geneid col (for some weird reason!!!) - bring in gene names manually

#Counts and Gene
Counts_gene <- cbind(Genelist, merged_counts)
Counts_gene2 <- Counts_gene[,c(2:11)]


#write the merged dataframe to a new file
write.table(Counts_gene2, file = "merged_featurecounts2.txt", sep = "/t", quote = FALSE, row.names = FALSE)

#Read in merged count file
count_data <- Counts_gene2  

#create metadata dataframe from experimental design
#plenti(wildtype) is A1-A3, KO22 is B1-B3, KO23 is C1-C3
condition <- factor(c(rep("plenti", 3), rep("KO22", 3), rep("KO23",3)))

#install biomaRt
BiocManager::install("biomaRt")
library(biomaRt)

#Using Ensembl as the database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#Fetch gene names corresponding to Ensembl IDs
ensembl_id_KO22 <- count_data[,1]
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensembl_id_KO22,
                      mart = mart)

#Change colname of gene_mapping to Geneid names for leftjoin to work
colnames(gene_mapping) <- c("Geneid", "Official gene symbol")
count_data <- left_join(count_data, gene_mapping, by = "Geneid")
count_data<-count_data %>% mutate("label"=case_when(is.na(`Official gene symbol`)==TRUE ~ Geneid,
                                                    is.na(`Official gene symbol`)==FALSE ~ `Official gene symbol`))

#Volcano Plot============================================================================
#results_plenti vs KO22===================================================
#Prepare data for plotting by converting p-value to -log10(pvalue) for y axis
results_plenti_vs_KO22 <- data.frame(count_data[,c(1,11,12)],pvalue = runif(nrow(count_data), min = 0, max = 1),log2FoldChange=runif(nrow(count_data), min = -5, max =10))
results_plenti_vs_KO22$logP <- -log10(results_plenti_vs_KO22$pvalue)

#Convert DESeqResults object to a dataframe for easier manipulation
results_plenti_vs_KO22_df <- as.data.frame(results_plenti_vs_KO22)
write.xlsx(results_plenti_vs_KO22_df,"test.xlsx")
merged_data <- results_plenti_vs_KO22_df[!(results_plenti_vs_KO22_df[, 2] == "" | is.na(results_plenti_vs_KO22_df[, 2])), ]
library(openxlsx)
write.xlsx(merged_data,"de0dena.xlsx")
write.xlsx(merged_data_kong,"de0.xlsx")
#Add gene name from count_data's "label" column to results_plenti_vs_KO22_df ?????? stuck here ?????
results_plenti_vs_KO22_df$label1 <- count_data$label

#Define significance thresholds
pval_cutoff <- 0.05
logfc_cutoff <- 2

#add column to define significance level for color
results_plenti_vs_KO22$significance <- ifelse(
  results_plenti_vs_KO22$pvalue < pval_cutoff & abs(results_plenti_vs_KO22$log2FoldChange) > logfc_cutoff,
  "Significant",
  "Not Significant"
)


merged_data$significance <- ifelse(
  merged_data$pvalue < pval_cutoff & abs(merged_data$log2FoldChange) > logfc_cutoff,
  "Significant",
  "Not Significant"
)

#Identify most significant genes, high siginificance of: p-value < 0.001 and log2FoldChange > 2
most_significant <- merged_data %>%
  filter(pvalue < 0.001 & abs(log2FoldChange) > 2) %>%
  arrange(pvalue) %>%
  head(5)
#select top 5 most significant genes

results_plenti_vs_KO22 <- data.frame(results_plenti_vs_KO22)
#Volcano plot using ggplot2
print(colnames(merged_data))

ggplot(merged_data, aes(x = log2FoldChange, y = logP, color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "blue") +
  geom_text_repel(data = most_significant,aes(label = label), size = 3, max.overlaps = 10)

#results_plenti vs KO23==========================================================
#Prepare data for plotting by converting p-value to -log10(pvalue) for y axis
results_plenti_vs_KO23$logP <- -log10(results_plenti_vs_KO23$pvalue)

#add column to define significance level for color
results_plenti_vs_KO23$significance <- ifelse(
  results_plenti_vs_KO23$pvalue < pval_cutoff & abs(results_plenti_vs_KO23$log2FoldChange) > logfc_cutoff,
  "Significant",
  "Not Significant"
)

#Volcano plot using ggplot2
ggplot(results_plenti_vs_KO23, aes(x = log2FoldChange, y = logP, color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_vline(xintercept = c(-logfc_cutoff, logfc_cutoff), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot of plenti vs KO23",
       x = "Log2 Fold Change",
       y = "-log10 P-value") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
