#PREPROCESSING==================================================================
#Installing the necessary packages
install.packages("ggplot2")# ggplot2: A system for declaratively creating graphics, based on The Grammar of Graphics.
install.packages("tidyverse")# tidyverse: A collection of R packages designed for data science.
install.packages("factoextra")
install.packages("ggpubr")
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
library(factoextra)
library(ggpubr)

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

#Filter genes with zero counts across all samples
nonzero_counts <- merged_counts[rowSums(merged_counts[, -1]) > 0, ]  # Exclude Geneid column
merged_counts <- nonzero_counts


#merged counts loses the geneid col (for some weird reason!!!) - bring in gene names manually

#Counts and Gene
Counts_gene <- cbind(Genelist, merged_counts)

# Assign GeneID (column 2) as row names
rownames(Counts_gene) <- Counts_gene[, 2] 

Counts_gene2 <- Counts_gene[,c(2:11)]


#write the merged dataframe to a new file
write.table(Counts_gene2, file = "merged_featurecounts2.txt", sep = "/t", quote = FALSE, row.names = FALSE)

#Read in merged count file
count_data <- Counts_gene2  

#create metadata dataframe from experimental design
#plenti(wildtype) is A1-A3, KO22 is B1-B3, KO23 is C1-C3
condition <- factor(c(rep("plenti", 3), rep("KO22", 3), rep("KO23",3)))



#Using Ensembl as the database
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") #sometimes the server is temporarily unavailble, go do something else

#Fetch gene names corresponding to Ensembl IDs
ensembl_id_KO22_KO23 <- count_data[,1]
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensembl_id_KO22_KO23,
                      mart = mart)

#Change colname of gene_mapping to Geneid names for leftjoin to work
colnames(gene_mapping) <- c("Geneid", "Official gene symbol")
count_data <- left_join(count_data, gene_mapping, by = "Geneid")
count_data<-count_data %>% mutate("label"=case_when(is.na(`Official gene symbol`)==TRUE ~ Geneid,
                                                    is.na(`Official gene symbol`)==FALSE ~ `Official gene symbol`))

#DESEQ2=========================================================================

#Create DESeq2 dataset object
# Create DESeq2 dataset object (with row names)
dds <- DESeqDataSetFromMatrix(
  countData = Counts_gene[, 3:11],  # Columns 3-11: count data
  colData = data.frame(condition),  # Metadata
  design = ~ condition               # Experimental design
)
#DESeq2 analysis
dds <- DESeq(dds)

#Transform the data for PCA
vsd <- vst(dds, blind = FALSE) #variance stabilizing transformation

#SCREE & PCA PLOT ==============================================================
##Scree plot -------------------------------------------------------------------
DS1.svd <- assay(vsd) |> 
  t() |> 
  prcomp(scale = FALSE) # PCA using prcomp()
summary(DS1.svd)

pScree <- fviz_eig(DS1.svd, addlabels = TRUE) + 
  theme_pubr(base_size = 9)

# PCA Plot with colors
pPCA <- fviz_pca_ind(DS1.svd, 
                     label = "all",  # Ensure all labels are visible
                     habillage = condition,  # Color by condition
                     repel = TRUE, # Prevent label overlap
                     mean.point = FALSE) +  #remove centroid marker
  labs(title = "PCA Plot",
       x = "PC1",
       y = "PC2")

# Arrange plots
pScreePCA <- ggarrange(pScree, pPCA,
                       labels = c("A", "B"),
                       ncol = 2, nrow = 1)

print(pScreePCA)
ggsave("PCA_ScreePlot.png", plot = pScreePCA, width = 20, height = 10)
#DESeq2 results for comparisons
#plenti vs KO22
results_plenti_vs_KO22 <- results(dds, contrast = c("condition", "KO22", "plenti")) %>%
  as.data.frame() %>%
  na.omit()  # Removes NA in padj

#plenti vs KO23
results_plenti_vs_KO23 <- results(dds, contrast = c("condition", "KO23", "plenti")) %>%
  as.data.frame() %>%
  na.omit()  # Removes NA in padj
# Save results
write.csv(as.data.frame(results_plenti_vs_KO22), "DEG_plenti_vs_KO22.csv")
write.csv(as.data.frame(results_plenti_vs_KO23), "DEG_plenti_vs_KO23.csv")

#Prepare GSEA ranked lists (.rnk files)
library(org.Hs.eg.db)

# For KO22 vs plenti
ranked_list_KO22 <- results_plenti_vs_KO22 %>%
  rownames_to_column("Geneid") %>%
  dplyr::mutate(
    EntrezID = mapIds(
      org.Hs.eg.db,
      keys = Geneid,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  ) %>%
  dplyr::filter(!is.na(EntrezID)) %>%
  dplyr::select(EntrezID, log2FoldChange) %>%
  dplyr::arrange(desc(log2FoldChange))

write.table(
  ranked_list_KO22,
  file = "GSEA_ranked_list_KO22.rnk",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

# For KO23 vs plenti (repeat similarly)
ranked_list_KO23 <- results_plenti_vs_KO23 %>%
  rownames_to_column("Geneid") %>%
  dplyr::mutate(
    EntrezID = mapIds(
      org.Hs.eg.db,
      keys = Geneid,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  ) %>%
  dplyr::filter(!is.na(EntrezID)) %>%
  dplyr::select(EntrezID, log2FoldChange) %>%
  dplyr::arrange(desc(log2FoldChange))

write.table(
  ranked_list_KO23,
  file = "GSEA_ranked_list_KO23.rnk",
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)


#VOLCANO PLOT===================================================================
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
Volcano_Plot_Plenti_v_KO22 <- ggplot(volcano_data_22, aes(x = log2FoldChange, y = -log10(padj), color = Significance, label = Gene)) +
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
ggsave("Volcano_Plot_Plenti_vs_KO22.png", plot = Volcano_Plot_Plenti_v_KO22, width = 10, height = 6)
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
Volcano_Plot_Plenti_v_KO23 <- ggplot(volcano_data_23, aes(x = log2FoldChange, y = -log10(padj), color = Significance, label = Gene)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) + 
  geom_vline(xintercept = c(-log2(2), log2(2)), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text(data = top_upregulated_23, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "red", check_overlap = TRUE) +
  geom_text(data = top_downregulated_23, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "blue", check_overlap = TRUE) +
  labs(title = "Volcano Plot (Plenti vs KO23) - Top 15 Up & Downregulated Genes", 
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()
ggsave("Volcano_Plot_Plenti_vs_KO23.png", plot = Volcano_Plot_Plenti_v_KO23, width = 10, height = 6)

#GENE ONTOLOGY==================================================================
# Load necessary libraries for GO analysis for KO22
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(DOSE)

# Define the function to convert fractions to decimals
fraction_to_decimal <- function(fraction_string) {
  sapply(strsplit(fraction_string, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

# Perform GO enrichment analysis for upregulated genes in different ontologies in KO22
upregulated_genes_KO22 <- rownames(results_plenti_vs_KO22)[results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange > log2(2)]

enrich.go.up.bp_KO22 <- enrichGO(gene = upregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.up.cc_KO22 <- enrichGO(gene = upregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.up.mf_KO22 <- enrichGO(gene = upregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

result.enrich.go.up.bp_KO22 <- enrich.go.up.bp_KO22@result
result.enrich.go.up.cc_KO22 <- enrich.go.up.cc_KO22@result
result.enrich.go.up.mf_KO22 <- enrich.go.up.mf_KO22@result

# Perform GO enrichment analysis for downregulated genes in different ontologies in KO22
downregulated_genes_KO22 <- rownames(results_plenti_vs_KO22)[results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange < -log2(2)]

enrich.go.down.bp_KO22 <- enrichGO(gene = downregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.down.cc_KO22 <- enrichGO(gene = downregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.down.mf_KO22 <- enrichGO(gene = downregulated_genes_KO22, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

result.enrich.go.down.bp_KO22 <- enrich.go.down.bp_KO22@result
result.enrich.go.down.cc_KO22 <- enrich.go.down.cc_KO22@result
result.enrich.go.down.mf_KO22 <- enrich.go.down.mf_KO22@result

# Convert GeneRatio from fraction to decimal for KO22
result.enrich.go.up.bp_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.up.bp_KO22$GeneRatio)
result.enrich.go.up.cc_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.up.cc_KO22$GeneRatio)
result.enrich.go.up.mf_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.up.mf_KO22$GeneRatio)

result.enrich.go.down.bp_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.down.bp_KO22$GeneRatio)
result.enrich.go.down.cc_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.down.cc_KO22$GeneRatio)
result.enrich.go.down.mf_KO22$GeneRatio <- fraction_to_decimal(result.enrich.go.down.mf_KO22$GeneRatio)

# Save GO enrichment results for KO22
write.csv(result.enrich.go.up.bp_KO22, "GO_enrichment_upregulated_BP_KO22.csv", row.names = FALSE)
write.csv(result.enrich.go.up.cc_KO22, "GO_enrichment_upregulated_CC_KO22.csv", row.names = FALSE)
write.csv(result.enrich.go.up.mf_KO22, "GO_enrichment_upregulated_MF_KO22.csv", row.names = FALSE)

write.csv(result.enrich.go.down.bp_KO22, "GO_enrichment_downregulated_BP_KO22.csv", row.names = FALSE)
write.csv(result.enrich.go.down.cc_KO22, "GO_enrichment_downregulated_CC_KO22.csv", row.names = FALSE)
write.csv(result.enrich.go.down.mf_KO22, "GO_enrichment_downregulated_MF_KO22.csv", row.names = FALSE)

# Filter and visualize top GO terms for KO22
top_15_up_bp_KO22 <- head(result.enrich.go.up.bp_KO22, 15)
top_15_up_cc_KO22 <- head(result.enrich.go.up.cc_KO22, 15)
top_15_up_mf_KO22 <- head(result.enrich.go.up.mf_KO22, 15)

top_15_down_bp_KO22 <- head(result.enrich.go.down.bp_KO22, 15)
top_15_down_cc_KO22 <- head(result.enrich.go.down.cc_KO22, 15)
top_15_down_mf_KO22 <- head(result.enrich.go.down.mf_KO22, 15)

# Dotplot for upregulated GO terms in BP for KO22
up_GO_dotplot_bp_KO22 <- ggplot(top_15_up_bp_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (BP) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_BP_KO22.png", plot = up_GO_dotplot_bp_KO22, width = 10, height = 6)

# Dotplot for upregulated GO terms in CC for KO22
up_GO_dotplot_cc_KO22 <- ggplot(top_15_up_cc_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (CC) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_CC_KO22.png", plot = up_GO_dotplot_cc_KO22, width = 10, height = 6)

# Dotplot for upregulated GO terms in MF for KO22
up_GO_dotplot_mf_KO22 <- ggplot(top_15_up_mf_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (MF) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_MF_KO22.png", plot = up_GO_dotplot_mf_KO22, width = 10, height = 6)

# Dotplot for downregulated GO terms in BP for KO22
down_GO_dotplot_bp_KO22 <- ggplot(top_15_down_bp_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (BP) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_BP_KO22.png", plot = down_GO_dotplot_bp_KO22, width = 10, height = 6)

# Dotplot for downregulated GO terms in CC for KO22
down_GO_dotplot_cc_KO22 <- ggplot(top_15_down_cc_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (CC) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_CC_KO22.png", plot = down_GO_dotplot_cc_KO22, width = 10, height = 6)

# Dotplot for downregulated GO terms in MF for KO22
down_GO_dotplot_mf_KO22 <- ggplot(top_15_down_mf_KO22, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (MF) KO22", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_MF_KO22.png", plot = down_GO_dotplot_mf_KO22, width = 10, height = 6)




# Perform GO enrichment analysis for upregulated genes in different ontologies in KO23
upregulated_genes_KO23 <- rownames(results_plenti_vs_KO23)[results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange > log2(2)]

enrich.go.up.bp_KO23 <- enrichGO(gene = upregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.up.cc_KO23 <- enrichGO(gene = upregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.up.mf_KO23 <- enrichGO(gene = upregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

result.enrich.go.up.bp_KO23 <- enrich.go.up.bp_KO23@result
result.enrich.go.up.cc_KO23 <- enrich.go.up.cc_KO23@result
result.enrich.go.up.mf_KO23 <- enrich.go.up.mf_KO23@result

# Perform GO enrichment analysis for downregulated genes in different ontologies in KO23
downregulated_genes_KO23 <- rownames(results_plenti_vs_KO23)[results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange < -log2(2)]

enrich.go.down.bp_KO23 <- enrichGO(gene = downregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.down.cc_KO23 <- enrichGO(gene = downregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)
enrich.go.down.mf_KO23 <- enrichGO(gene = downregulated_genes_KO23, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = TRUE)

result.enrich.go.down.bp_KO23 <- enrich.go.down.bp_KO23@result
result.enrich.go.down.cc_KO23 <- enrich.go.down.cc_KO23@result
result.enrich.go.down.mf_KO23 <- enrich.go.down.mf_KO23@result

# Convert GeneRatio from fraction to decimal for KO23
result.enrich.go.up.bp_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.up.bp_KO23$GeneRatio)
result.enrich.go.up.cc_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.up.cc_KO23$GeneRatio)
result.enrich.go.up.mf_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.up.mf_KO23$GeneRatio)

result.enrich.go.down.bp_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.down.bp_KO23$GeneRatio)
result.enrich.go.down.cc_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.down.cc_KO23$GeneRatio)
result.enrich.go.down.mf_KO23$GeneRatio <- fraction_to_decimal(result.enrich.go.down.mf_KO23$GeneRatio)

# Save GO enrichment results for KO23
write.csv(result.enrich.go.up.bp_KO23, "GO_enrichment_upregulated_BP_KO23.csv", row.names = FALSE)
write.csv(result.enrich.go.up.cc_KO23, "GO_enrichment_upregulated_CC_KO23.csv", row.names = FALSE)
write.csv(result.enrich.go.up.mf_KO23, "GO_enrichment_upregulated_MF_KO23.csv", row.names = FALSE)

write.csv(result.enrich.go.down.bp_KO23, "GO_enrichment_downregulated_BP_KO23.csv", row.names = FALSE)
write.csv(result.enrich.go.down.cc_KO23, "GO_enrichment_downregulated_CC_KO23.csv", row.names = FALSE)
write.csv(result.enrich.go.down.mf_KO23, "GO_enrichment_downregulated_MF_KO23.csv", row.names = FALSE)

# Filter and visualize top GO terms for KO23
top_15_up_bp_KO23 <- head(result.enrich.go.up.bp_KO23, 15)
top_15_up_cc_KO23 <- head(result.enrich.go.up.cc_KO23, 15)
top_15_up_mf_KO23 <- head(result.enrich.go.up.mf_KO23, 15)

top_15_down_bp_KO23 <- head(result.enrich.go.down.bp_KO23, 15)
top_15_down_cc_KO23 <- head(result.enrich.go.down.cc_KO23, 15)
top_15_down_mf_KO23 <- head(result.enrich.go.down.mf_KO23, 15)

# Dotplot for upregulated GO terms in BP for KO23
up_GO_dotplot_bp_KO23 <- ggplot(top_15_up_bp_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (BP) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_BP_KO23.png", plot = up_GO_dotplot_bp_KO23, width = 10, height = 6)

# Dotplot for upregulated GO terms in CC for KO23
up_GO_dotplot_cc_KO23 <- ggplot(top_15_up_cc_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (CC) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_CC_KO23.png", plot = up_GO_dotplot_cc_KO23, width = 10, height = 6)

# Dotplot for upregulated GO terms in MF for KO23
up_GO_dotplot_mf_KO23 <- ggplot(top_15_up_mf_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Upregulated GO Terms (MF) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_upregulated_dotplot_MF_KO23.png", plot = up_GO_dotplot_mf_KO23, width = 10, height = 6)

# Dotplot for downregulated GO terms in BP for KO23
down_GO_dotplot_bp_KO23 <- ggplot(top_15_down_bp_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (BP) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_BP_KO23.png", plot = down_GO_dotplot_bp_KO23, width = 10, height = 6)

# Dotplot for downregulated GO terms in CC for KO23
down_GO_dotplot_cc_KO23 <- ggplot(top_15_down_cc_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (CC) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_CC_KO23.png", plot = down_GO_dotplot_cc_KO23, width = 10, height = 6)

# Dotplot for downregulated GO terms in MF for KO23
down_GO_dotplot_mf_KO23 <- ggplot(top_15_down_mf_KO23, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  labs(title = "Top 15 Downregulated GO Terms (MF) KO23", x = "Gene Ratio", y = "GO Term") +
  theme_minimal()
ggsave("GO_downregulated_dotplot_MF_KO23.png", plot = down_GO_dotplot_mf_KO23, width = 10, height = 6)

#KEGG Enrichment Analysis (there are two versions of this, same result?)========
# Function to Convert Gene Symbols to Entrez IDs
convert_to_entrez <- function(gene_list) {
  entrez_ids <- mapIds(org.Hs.eg.db, 
                       keys = gene_list, 
                       column = "ENTREZID", 
                       keytype = "SYMBOL", 
                       multiVals = "first")
  
  # Remove NA values
  return(entrez_ids[!is.na(entrez_ids)])
}

# Function to Perform KEGG Enrichment
perform_kegg_enrichment <- function(gene_list) {
  entrez_ids <- convert_to_entrez(gene_list)  # Reuse the convert_to_entrez() function
  enrichKEGG(gene = entrez_ids,
             organism = "hsa",  # Human KEGG pathways
             pvalueCutoff = 0.05,
             keyType = "kegg")
}

# Function to Plot KEGG Results with ggplot2
plot_kegg_ggplot <- function(kegg_res, title, filename, top_n = 20) {
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res@result)) > 0) {
    kegg_df <- as.data.frame(kegg_res@result)
    kegg_df <- head(kegg_df[order(kegg_df$p.adjust), ], top_n)
    
    # Calculate GeneRatio as a numeric value
    kegg_df$GeneRatio <- sapply(kegg_df$GeneRatio, function(x) eval(parse(text = x)))
    
    p <- ggplot(kegg_df, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
      geom_point() +
      scale_color_gradient(low = "red", high = "blue") +
      labs(title = title, x = "Gene Ratio", y = "KEGG Pathway", color = "Adjusted p-value", size = "Gene Count") +
      theme_minimal() +
      theme(text = element_text(size = 12))
    
    # Save plot
    ggsave(filename, plot = p, width = 10, height = 6)
  } else {
    print(paste("No significant KEGG pathways found for", title))
  }
}

# Convert Gene Symbols to Entrez IDs
up_entrez_KO22 <- convert_to_entrez(upregulated_genes_KO22)
down_entrez_KO22 <- convert_to_entrez(downregulated_genes_KO22)

# Perform KEGG Enrichment for Upregulated and Downregulated Genes
kegg_up_KO22 <- enrichKEGG(gene = up_entrez_KO22, organism = "hsa", pvalueCutoff = 0.05)
kegg_down_KO22 <- enrichKEGG(gene = down_entrez_KO22, organism = "hsa", pvalueCutoff = 0.05)

# Plot and Save KEGG Results
plot_kegg_ggplot(kegg_up_KO22, "KEGG Pathway for KO22 - Upregulated Genes", "KEGG_Upregulated_KO22.png")
plot_kegg_ggplot(kegg_down_KO22, "KEGG Pathway for KO22 - Downregulated Genes", "KEGG_Downregulated_KO22.png")

# Save Results
write.csv(as.data.frame(kegg_up_KO22@result), "KEGG_Upregulated_KO22.csv")
write.csv(as.data.frame(kegg_down_KO22@result), "KEGG_Downregulated_KO22.csv")

# Convert Gene Symbols to Entrez IDs
up_entrez_KO23 <- convert_to_entrez(upregulated_genes_KO23)
down_entrez_KO23 <- convert_to_entrez(downregulated_genes_KO23)

# Perform KEGG Enrichment for Upregulated and Downregulated Genes
kegg_up_KO23 <- enrichKEGG(gene = up_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)
kegg_down_KO23 <- enrichKEGG(gene = down_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)

# Plot and Save KEGG Results
plot_kegg_ggplot(kegg_up_KO23, "KEGG Pathway for KO23 - Upregulated Genes", "KEGG_Upregulated_KO23.png")
plot_kegg_ggplot(kegg_down_KO23, "KEGG Pathway for KO23 - Downregulated Genes", "KEGG_Downregulated_KO23.png")

# Save Results
write.csv(as.data.frame(kegg_up_KO23@result), "KEGG_Upregulated_KO23.csv")
write.csv(as.data.frame(kegg_down_KO23@result), "KEGG_Downregulated_KO23.csv")

