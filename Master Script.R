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
library(EnhancedVolcano)  # For creating volcano plots
library(biomaRt)          # For accessing BioMart databases

#Define file Path
files <- list.files(files <- list.files("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts", # nolint
                        pattern = "*.txt", full.names = TRUE))

#read each file into a list of dataframes
counts_list <- lapply(files, function(file) {
  df <- read.table(file, header = TRUE)
  df <- df %>% select(Geneid, ends_with(".bam"))
  colnames(df) <- c("Geneid", basename(file))
  return(df)
})

View(counts_list[[1]]) ## do i really need this?

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
