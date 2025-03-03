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

#install biomaRt
BiocManager::install("biomaRt")

#install vscDebugger for VSCode
install.packages("vscDebugger", repos = "https://manuelhentschel.r-universe.dev")

#Load necessary libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(EnhancedVolcano)
library(biomaRt)

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

#PCA PLOT ===========================

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

#Generate Volcano Plots
#Volcano plot of plenti vs KO22
rownames(results_plenti_vs_KO22) <- count_data$label
EnhancedVolcano(results_plenti_vs_KO22,
                lab = rownames(results_plenti_vs_KO22),
                x = "log2FoldChange",
                y = "pvalue",
                title = "plenti vs KO22",
                pCutoff = 0.05,
                FCcutoff = 2.0,
                labSize = 3,
                selectLab = rownames(results_plenti_vs_KO22)[which(abs(results_plenti_vs_KO22$log2FoldChange) > 2 & results_plenti_vs_KO22$pvalue < 0.01)])
                
#Volcano plot of plenti vs KO23
rownames(results_plenti_vs_KO23) <- count_data$label
EnhancedVolcano(results_plenti_vs_KO23,
                lab = rownames(results_plenti_vs_KO23),
                x = "log2FoldChange",
                y = "pvalue",
                title = "plenti vs KO23",
                pCutoff = 0.05,
                FCcutoff = 1.5,
                labSize = 3,
                selectLab = rownames(results_plenti_vs_KO23)[which(abs(results_plenti_vs_KO23$log2FoldChange) > 2 & results_plenti_vs_KO23$pvalue < 0.01)])


