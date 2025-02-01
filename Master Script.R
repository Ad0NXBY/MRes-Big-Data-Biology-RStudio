#Installing the necessary packages
install.packages("ggplot2")# ggplot2: A system for declaratively creating graphics, based on The Grammar of Graphics.
install.packages("tidyverse")# tidyverse: A collection of R packages designed for data science.
if (!require("BiocManager", quietly = TRUE))# BiocManager: A package to manage the installation of Bioconductor packages.
  install.packages("BiocManager")
BiocManager::install("GenomicFeatures")# GenomicFeatures: Tools for making and manipulating transcript centric annotations.
BiocManager::install("DESeq2")# DESeq2: Differential gene expression analysis based on the negative binomial distribution.
BiocManager::install("EnhancedVolcano")# EnhancedVolcano: Publication-ready volcano plots.
BiocManager::install("biomaRt")# biomaRt: Interface to BioMart databases.

# Load necessary libraries
library(DESeq2)           # For differential gene expression analysis
library(dplyr)            # For data manipulation
library(ggplot2)          # For plotting
library(tidyverse)        # For data science tasks
library(biomaRt)          # For accessing BioMart databases

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

#DESEQ2======================================================================

#Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = Counts_gene[,3:11],
                              colData = data.frame(condition),
                              design = ~condition)

#DESeq2 analysis
dds <- DESeq(dds)
                  
#Transform the data for PCA
vsd <- vst(dds, blind = FALSE) #variance stabilizing transformation

#PCA PLOT ====================================================================

#Generate a PCA plot
plotPCA(vsd, intgroup = "condition") +
  ggtitle("PCA Plot of RNA-Seq Data")

#DESeq2 results for comparisons
#plenti vs KO22
results_plenti_vs_KO22 <- results(dds, contrast = c("condition", "KO22", "plenti"))
summary(results_plenti_vs_KO22) #summary of results
#plenti vs KO23
results_plenti_vs_KO23 <- results(dds, contrast = c("condition", "KO23", "plenti"))
summary(results_plenti_vs_KO23) #summary of results


#Volcano Plot==============================================================
##Volcano plot of Plent vs KO22---------------------------
###Ensure gene names are set as rownames
rownames(results_plenti_vs_KO22) <- count_data$label
volcano_data_22 <- as.data.frame(results_plenti_vs_KO22)
volcano_data_22$Gene <- rownames(volcano_data_22)

###Define significance thresholds
volcano_data_22$Significance <- "Not Significant"
volcano_data_22$Significance[volcano_data_22$padj < 0.05 & volcano_data_22$log2FoldChange > log2(2)] <- "Upregulated"
volcano_data_22$Significance[volcano_data_22$padj < 0.05 & volcano_data_22$log2FoldChange < -log2(2)] <- "Downregulated"

###Separate top 15 upregulated and downregulated genes
top_upregulated_22 <- volcano_data_22 %>%
  filter(Significance == "Upregulated") %>%
  arrange(padj) %>%
  head(15)

top_downregulated_22 <- volcano_data_22 %>%
  filter(Significance == "Downregulated") %>%
  arrange(padj) %>%
  head(15)

###Combine top 30 genes
top_genes_22 <- rbind(top_upregulated_22, top_downregulated_22)

###Create the volcano plot
ggplot(volcano_data_22, aes(x = log2FoldChange, y = -log10(padj), color = Significance, label = Gene)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) + 
  geom_vline(xintercept = c(-log2(2), log2(2)), linetype = "dashed", color = "black") + #FC threshold
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # p-value threshold
  geom_text(data = top_upregulated_22, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "red", check_overlap = TRUE) + #Red for upregulated
  geom_text(data = top_downregulated_22, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "blue", check_overlap = TRUE) + #Blue for downregulated
  xlim(-10, 10) +  # Limit the Log2 Fold Change range as it will go up to 20
  labs(title = "Volcano Plot (Plenti vs KO22) - Top 15 Up & Downregulated Genes", 
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()

##Volcano plot of Plent vs KO23---------------------------
###Ensure gene names are set as rownames
rownames(results_plenti_vs_KO23) <- count_data$label
volcano_data_23 <- as.data.frame(results_plenti_vs_KO23)
volcano_data_23$Gene <- rownames(volcano_data_23)

###Define significance thresholds
volcano_data_23$Significance <- "Not Significant"
volcano_data_23$Significance[volcano_data_23$padj < 0.05 & volcano_data_23$log2FoldChange > log2(2)] <- "Upregulated"
volcano_data_23$Significance[volcano_data_23$padj < 0.05 & volcano_data_23$log2FoldChange < -log2(2)] <- "Downregulated"

###Top 15 up/downregulated
top_upregulated_23 <- volcano_data_23 %>%
  filter(Significance == "Upregulated") %>%
  arrange(padj) %>%
  head(15)

top_downregulated_23 <- volcano_data_23 %>%
  filter(Significance == "Downregulated") %>%
  arrange(padj) %>%
  head(15)

###Combine top 30 genes
top_genes_23 <- rbind(top_upregulated_23, top_downregulated_23)

###Volcano plot for KO23
ggplot(volcano_data_23, aes(x = log2FoldChange, y = -log10(padj), color = Significance, label = Gene)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) + 
  geom_vline(xintercept = c(-log2(2), log2(2)), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text(data = top_upregulated_23, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "red", check_overlap = TRUE) +
  geom_text(data = top_downregulated_23, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "blue", check_overlap = TRUE) +
  labs(title = "Volcano Plot (Plenti vs KO23) - Top 15 Up & Downregulated Genes", 
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()
