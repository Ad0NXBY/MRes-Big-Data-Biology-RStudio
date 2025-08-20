#PREPROCESSING==================================================================
##Installing the necessary packages---------------------------------------------
install.packages("ggplot2")        # For advanced plotting using the Grammar of Graphics.
install.packages("tidyverse")      # For data manipulation and analysis (dplyr, tidyr, etc.).
install.packages("factoextra")     # For visualizing PCA results and multivariate data.
install.packages("ggpubr")         # For publication-ready plots.
install.packages("data.table")     # For fast file reading with fread().
install.packages("limma")          # For filtering lowly expressed genes (filterByExpr()).
install.packages("ggrepel")        # For improved text labeling on plots.
install.packages("VennDiagram")    # For creating Venn diagrams.
install.packages("scales")         # For scaling functions in plots.

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("GenomicFeatures")   # For handling and annotating genomic features.
BiocManager::install("DESeq2")            # For differential expression analysis.
BiocManager::install("EnhancedVolcano")   # For generating enhanced volcano plots.
BiocManager::install("biomaRt")           # For interfacing with BioMart (e.g., Ensembl).
BiocManager::install("edgeR")             # For filtering low expressed genes using filterByExpr().
BiocManager::install("clusterProfiler")   # For GO and KEGG enrichment analyses.
BiocManager::install("enrichplot")        # For visualizing enrichment results.
BiocManager::install("DOSE")              # For disease ontology and enrichment analysis.
BiocManager::install("ReactomePA")        # For Reactome pathway analysis.


##Load Necessary Libraries------------------------------------------------------
# General data handling and plotting packages:
library(ggplot2)        
library(tidyverse)      
library(factoextra)   
library(ggpubr)         
library(data.table)   
library(limma)        
library(ggrepel)        
library(VennDiagram)   
library(scales)         

# Differential expression and genomic analysis packages:
library(DESeq2)         
library(EnhancedVolcano) 
library(biomaRt)        
library(edgeR)          

# Enrichment and pathway analysis packages:
library(clusterProfiler) 
library(enrichplot)     
library(DOSE)            
library(ReactomePA)    

##Create output directories-----------------------------------------------------
plot_dir <- "Analysis plots with combined code with apply"
data_dir <- "List_for_analysis with combined code with apply"

dir.create(plot_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)

##Read files into R-------------------------------------------------------------
Control1 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_A1_featurecounts.txt")
Control2 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_A2_featurecounts.txt")
Control3 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_A3_featurecounts.txt")
FIH22KO1 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_B1_featurecounts.txt")
FIH22KO2 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_B2_featurecounts.txt")
FIH22KO3 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_B3_featurecounts.txt")
FIH23KO1 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_C1_featurecounts.txt")
FIH23KO2 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_C2_featurecounts.txt")
FIH23KO3 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_C3_featurecounts.txt")

##Set column names--------------------------------------------------------------
colnames(Control1)[7] <- "Control1_Count"
colnames(Control2)[7] <- "Control2_Count"
colnames(Control3)[7] <- "Control3_Count"
colnames(FIH22KO1)[7] <- "FIH22KO1_Count"
colnames(FIH22KO2)[7] <- "FIH22KO2_Count"
colnames(FIH22KO3)[7] <- "FIH22KO3_Count"
colnames(FIH23KO1)[7] <- "FIH23KO1_Count"
colnames(FIH23KO2)[7] <- "FIH23KO2_Count"
colnames(FIH23KO3)[7] <- "FIH23KO3_Count"
colnames(Control1)[1] <- "GeneID"
colnames(Control2)[1] <- "GeneID"
colnames(Control3)[1] <- "GeneID"
colnames(FIH22KO1)[1] <- "GeneID"
colnames(FIH22KO2)[1] <- "GeneID"
colnames(FIH22KO3)[1] <- "GeneID"
colnames(FIH23KO1)[1] <- "GeneID"
colnames(FIH23KO2)[1] <- "GeneID"
colnames(FIH23KO3)[1] <- "GeneID"

##Merge data--------------------------------------------------------------------
combined_data <- Control1[, .(GeneID, Control1_Count)]
combined_data <- merge(combined_data, Control2[, .(GeneID, Control2_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, Control3[, .(GeneID, Control3_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH22KO1[, .(GeneID, FIH22KO1_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH22KO2[, .(GeneID, FIH22KO2_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH22KO3[, .(GeneID, FIH22KO3_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH23KO1[, .(GeneID, FIH23KO1_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH23KO2[, .(GeneID, FIH23KO2_Count)], by = "GeneID", all = TRUE)
combined_data <- merge(combined_data, FIH23KO3[, .(GeneID, FIH23KO3_Count)], by = "GeneID", all = TRUE)
combined_data <- as.data.frame(combined_data)
rownames(combined_data) <- combined_data$GeneID

##Filter genes with zero counts across all samples------------------------------(USING APPLY, REMOVES ROWS WITH EVEN 1!!!! ZERO COUNT)
combined_data <- combined_data[!apply(combined_data == 0, 1, any), ]

##Using Ensembl as the database-------------------------------------------------
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensemblID_to_genesymbol <- combined_data[, 1]

gene_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensemblID_to_genesymbol,
                      mart = mart)

colnames(gene_mapping) <- c("GeneID", "GeneSymbol")

combined_data <- as.data.frame(combined_data)
combined_data$GeneID <- rownames(combined_data)

#Merge with gene mapping
combined_data <- left_join(combined_data, gene_mapping, by = "GeneID")

#Replace GeneID with GeneSymbol if available
combined_data$FinalName <- ifelse(is.na(combined_data$GeneSymbol), combined_data$GeneID, combined_data$GeneSymbol)

#Keep only the first occurrence of each gene symbol
combined_data <- combined_data %>%
  distinct(FinalName, .keep_all = TRUE)

#Set row names and remove unnecessary columns
rownames(combined_data) <- combined_data$FinalName
combined_data <- combined_data %>% select(-GeneID, -GeneSymbol, -FinalName)



##Convert the data to a DGEList object for edgeR to work------------------------
counts <- as.matrix(combined_data)
condition <- factor(c(rep("plenti", 3), rep("KO22", 3), rep("KO23",3)))
dge <- DGEList(counts = counts, group = condition)

#Filter out low-expressed genes
keep <- filterByExpr(dge, group = condition, min.count = 15, min.total.count = 70, min.prop = 0.2)
dge_filtered <- dge[keep,]
filtered_data <- as.data.frame(dge_filtered$counts)

##Conver the data to the format for DESeq2--------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(filtered_data),
  colData = data.frame(condition),
  design = ~condition
)

#Run DESeq2 normalization and differentail expression analysis
dds <- DESeq(dds)


##Testing for normal distribution-----------------------------------------------
normalized_counts <- assay(rlog(dds))
apply(normalized_counts, 2, function(column_data) {
  ks.test(column_data, "pnorm", mean = mean(column_data), sd = sd(column_data))
})

##Test for Poisson distribution-------------------------------------------------
#Fit a Poisson distribution for each sample and calculate the likelihood
poisson_fit <- function(column_data) {
  lambda <- mean(column_data)  # The λ for Poisson distribution is the mean of the data
  log_likelihood <- sum(dpois(column_data, lambda, log = TRUE))  # Calculate log-likelihood of Poisson distribution
  return(log_likelihood)
}

#Perform the fit for each sample
log_likelihood_results <- apply(normalized_counts, 2, poisson_fit)

#View results
log_likelihood_results

#Save results
write.csv(normalized_counts, file.path(data_dir, "normalized_counts.csv"))

#DEG============================================================================
experimental_data <- read.csv("C:/Users/Brandon/Documents/MRes Big Data Biology/R studio/List_for_analysis with combined code//normalized_counts.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(experimental_data) <- experimental_data$X
experimental_data <- experimental_data[,-1]

# Filter low-expressed genes using limma's filterByExpr()
keep_genes <- filterByExpr(experimental_data, 
                           group = condition, 
                           min.count = 15, 
                           min.total.count = 70, 
                           min.prop = 0.2)

# Subset the data (correct way)
experimental_data_filtered <- experimental_data[keep_genes, ]

condition <- factor(c(rep("plenti", 3), rep("KO22", 3), rep("KO23", 3)))


# 创建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(experimental_data),
  colData = data.frame(condition = condition),
  design = ~ condition
)

dds <- DESeq(dds)

#Creating .csv DEG list for FIHKO22
results_plenti_vs_KO22 <- results(dds, contrast = c("condition", "KO22", "plenti"))
results_plenti_vs_KO22_df <- as.data.frame(results_plenti_vs_KO22)
results_plenti_vs_KO22_df$gene <- rownames(results_plenti_vs_KO22_df)

# Save full DEG results
write.csv(results_plenti_vs_KO22_df, 
          file.path(data_dir, "DEG_KO22_vs_Plenti_filterByExpr.csv"), 
          row.names = TRUE)

# Filter for significant DEGs for GSEA (padj < 0.05) ==========================
GSEA_FIH22KO <- results_plenti_vs_KO22_df %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(desc(log2FoldChange))

# Save the filtered DEG list
write.csv(GSEA_FIH22KO,
          file.path(data_dir, "GSEA_plenti_vs_KO22.csv"),
          row.names = TRUE)

# Prepare ranked list for GSEA ===============================================
GSEA_ranked_list_KO22 <- GSEA_FIH22KO %>%
  dplyr::select(gene, log2FoldChange)

# Write .rnk file (no headers or row names as required by GSEA)
write.table(
  GSEA_ranked_list_KO22,
  file.path(data_dir, "GSEA_ranked_list_KO22_filter_by_padj_lessthan0.05.rnk"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)



# 提取 FIH23KO vs Control 的差异表达结果
results_plenti_vs_KO23 <- results(dds, contrast = c("condition", "KO23", "plenti"))
results_plenti_vs_KO23_df <- as.data.frame(results_plenti_vs_KO23)
results_plenti_vs_KO23_df$gene <- rownames(results_plenti_vs_KO23_df)

# Save full DEG results
write.csv(results_plenti_vs_KO23_df, 
          file.path(data_dir, "DEG_KO23_vs_Plenti_filterByExpr.csv"), 
          row.names = TRUE)

# Filter for significant DEGs for GSEA (padj < 0.05) ==========================
GSEA_FIH23KO <- results_plenti_vs_KO23_df %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(desc(log2FoldChange))

# Save the filtered DEG list
write.csv(GSEA_FIH23KO,
          file.path(data_dir, "GSEA_plenti_vs_KO23.csv"),
          row.names = TRUE)

# Prepare ranked list for GSEA ===============================================
GSEA_ranked_list_KO23 <- GSEA_FIH23KO %>%
  dplyr::select(gene, log2FoldChange)

# Write .rnk file (no headers or row names as required by GSEA)
write.table(
  GSEA_ranked_list_KO23,
  file.path(data_dir, "GSEA_ranked_list_KO23_filter_by_padj_lessthan0.05.rnk"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

#Scree & PCA Plot===========================================================================
#Transform the data for PCA
vsd <- vst(dds, blind = FALSE) #variance stabilizing transformation

condition <- factor(c(rep("plenti", 3), rep("KO22", 3), rep("KO23",3)))
##Scree plot -------------------------------------------------------------------
DS1.svd <- assay(vsd) |> 
  t() |> 
  prcomp(scale = FALSE) # PCA using prcomp()
summary(DS1.svd)

pScree <- fviz_eig(DS1.svd, addlabels = TRUE) + 
  theme_pubr(base_size = 9)

##PCA Plot with colors----------------------------------------------------------
pPCA <- fviz_pca_ind(DS1.svd, 
                     label = "all",  # Ensure all labels are visible
                     habillage = condition,  # Color by condition
                     repel = TRUE, # Prevent label overlap
                     mean.point = FALSE) +  #remove centroid marker
  labs(title = "PCA Plot",
       x = "PC1",
       y = "PC2")

#Arrange plots
pScreePCA <- ggarrange(pScree, pPCA,
                       labels = c("A", "B"),
                       ncol = 2, nrow = 1)

print(pScreePCA)
ggsave(file.path(plot_dir, "PCA_ScreePlot.png"), plot = pScreePCA, width = 20, height = 10)




#VENN PLOT=================================================================================
up_genes_KO22 <- rownames(results_plenti_vs_KO22[!is.na(results_plenti_vs_KO22$padj) & results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange > 0, ])
down_genes_KO22 <- rownames(results_plenti_vs_KO22[!is.na(results_plenti_vs_KO22$padj) & results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange < 0, ])

up_genes_KO23 <- rownames(results_plenti_vs_KO23[!is.na(results_plenti_vs_KO23$padj) & results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange > 0, ])
down_genes_KO23 <- rownames(results_plenti_vs_KO23[!is.na(results_plenti_vs_KO23$padj) & results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange < 0, ])

gene_lists <- list(
  "FIH22KO_up" = up_genes_KO22,
  "FIH22KO_down" = down_genes_KO22,
  "FIH23KO_up" = up_genes_KO23,
  "FIH23KO_down" = down_genes_KO23
)

venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = c("KO22 Upregulated", "KO22 Downregulated", "KO23 Upregulated", "KO23 Downregulated"),
  filename = NULL,
  imagetype = "png",
  height = 750,
  width = 750,
  resolution = 300,
  compression = "lzw",
  lwd = 2,
  col = "transparent",
  fill = c("#F8766D", "#00BA38", "#619CFF", "purple4"),
  alpha = 0.50,
  label.col = "darkblue",
  cex = 2,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkred", "darkgreen", "darkblue", "purple4"),
  cat.cex = 1.3,
  cat.fontfamily = "serif",
  cat.default.pos = "outer",
  cat.pos = c(-40, 40, -30, 30),
  cat.dist = c(0.3, 0.3, 0.2, 0.2),
  cat.fontface = "bold",
  rotation.degree = 0,
  margin = 0.2
)
grid.draw(venn.plot)


# 绘制维恩图
ggsave(file.path(plot_dir, "vennplot.jpg"), plot = venn.plot, width = 12, height = 8, dpi = 800)


#VOLCANO PLOT===================================================================
##Volcano plot of Plent vs KO22-------------------------------------------------
#Ensure gene names are set as rownames
volcano_data_22 <- as.data.frame(results_plenti_vs_KO22)
volcano_data_22$Gene <- rownames(volcano_data_22)

#Define significance thresholds
volcano_data_22$Significance <- "Not Significant"
volcano_data_22$Significance[volcano_data_22$padj < 0.05 & volcano_data_22$log2FoldChange > log2(2)] <- "Upregulated"
volcano_data_22$Significance[volcano_data_22$padj < 0.05 & volcano_data_22$log2FoldChange < -log2(2)] <- "Downregulated"

#Separate top 25 upregulated and downregulated genes
top_upregulated_22 <- volcano_data_22 %>%
  filter(Significance == "Upregulated") %>%
  arrange(padj) %>%
  head(25)

top_downregulated_22 <- volcano_data_22 %>%
  filter(Significance == "Downregulated") %>%
  arrange(padj) %>%
  head(25)

#Combine top 30 genes
top_genes_22 <- rbind(top_upregulated_22, top_downregulated_22)

#Create the volcano plot
Volcano_Plot_Plenti_v_KO22 <- ggplot(volcano_data_22, aes(x = log2FoldChange, y = -log10(padj), color = Significance, label = Gene)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) + 
  geom_vline(xintercept = c(-log2(2), log2(2)), linetype = "dashed", color = "black") + #FC threshold
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # p-value threshold
  geom_text(data = top_upregulated_22, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "black", check_overlap = TRUE) + #Red for upregulated
  geom_text(data = top_downregulated_22, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "black", check_overlap = TRUE) + #Blue for downregulated
  xlim(-10, 10) +  # Limit the Log2 Fold Change range as it will go up to 20
  labs(title = "Volcano Plot (Plenti vs KO22) - Top 25 Up & Downregulated Genes", 
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()
ggsave(file.path(plot_dir, "Volcano_Plot_Plenti_vs_KO22.png"), plot = Volcano_Plot_Plenti_v_KO22, width = 10, height = 6)
##Volcano plot of Plent vs KO23-------------------------------------------------
#Ensure gene names are set as rownames
volcano_data_23 <- as.data.frame(results_plenti_vs_KO23)
volcano_data_23$Gene <- rownames(volcano_data_23)

#Define significance thresholds
volcano_data_23$Significance <- "Not Significant"
volcano_data_23$Significance[volcano_data_23$padj < 0.05 & volcano_data_23$log2FoldChange > log2(2)] <- "Upregulated"
volcano_data_23$Significance[volcano_data_23$padj < 0.05 & volcano_data_23$log2FoldChange < -log2(2)] <- "Downregulated"

#Top 25 up/downregulated
top_upregulated_23 <- volcano_data_23 %>%
  filter(Significance == "Upregulated") %>%
  arrange(padj) %>%
  head(25)

top_downregulated_23 <- volcano_data_23 %>%
  filter(Significance == "Downregulated") %>%
  arrange(padj) %>%
  head(25)

#Combine top 30 genes
top_genes_23 <- rbind(top_upregulated_23, top_downregulated_23)

#Volcano plot for KO23
Volcano_Plot_Plenti_v_KO23 <- ggplot(volcano_data_23, aes(x = log2FoldChange, y = -log10(padj), color = Significance, label = Gene)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) + 
  geom_vline(xintercept = c(-log2(2), log2(2)), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_text(data = top_upregulated_23, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "black", check_overlap = TRUE) +
  geom_text(data = top_downregulated_23, aes(label = Gene), vjust = 1, hjust = 1, size = 3, color = "black", check_overlap = TRUE) +
  labs(title = "Volcano Plot (Plenti vs KO23) - Top 5 Up & Downregulated Genes", 
       x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
  theme_minimal()
ggsave(file.path(plot_dir,"Volcano_Plot_Plenti_vs_KO23.png"), plot = Volcano_Plot_Plenti_v_KO23, width = 10, height = 6)


#GENE ONTOLOGY=======================================================================================
# FIH22KO 确保筛选时排除掉 NA 值--------------------------------------------------------------------
up_genes_FIH22KO <- rownames(results_plenti_vs_KO22[!is.na(results_plenti_vs_KO22$padj) & results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange > 0, ])
down_genes_FIH22KO <- rownames(results_plenti_vs_KO22[!is.na(results_plenti_vs_KO22$padj) & results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange < 0, ])
background_genes <- rownames(experimental_data)

# 上调基因的GO分析
ego_up_KO22 <- enrichGO(gene         = up_genes_FIH22KO,
                        OrgDb        = org.Hs.eg.db,
                        keyType      = "SYMBOL",
                        ont          = "ALL",  
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)

ego_up_df_KO22 <- as.data.frame(ego_up_KO22)

# 下调基因的GO分析
ego_down_KO22 <- enrichGO(gene         = down_genes_FIH22KO,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",
                          ont          = "ALL",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

ego_down_df_KO22 <- as.data.frame(ego_down_KO22)

# 上调基因表达的GO分析结果处理

sorted_up_bp_KO22 <- ego_up_df_KO22 %>%
  filter(ONTOLOGY == "BP") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_up_cc_KO22 <- ego_up_df_KO22 %>%
  filter(ONTOLOGY == "CC") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_up_mf_KO22 <- ego_up_df_KO22 %>%
  filter(ONTOLOGY == "MF") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_ego_up_df_KO22 <- bind_rows(sorted_up_bp_KO22, sorted_up_cc_KO22, sorted_up_mf_KO22)

# 下调基因表达的GO分析结果处理，与上调处理相似

sorted_down_bp_KO22 <- ego_down_df_KO22 %>%
  filter(ONTOLOGY == "BP") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_down_cc_KO22 <- ego_down_df_KO22 %>%
  filter(ONTOLOGY == "CC") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_down_mf_KO22 <- ego_down_df_KO22 %>%
  filter(ONTOLOGY == "MF") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_ego_down_df_KO22 <- bind_rows(sorted_down_bp_KO22, sorted_down_cc_KO22, sorted_down_mf_KO22)


# 上调基因的点图
Gene_ontology_KO22_UP <- ggplot(sorted_ego_up_df_KO22, aes(x = Count, y = Description)) +
  geom_point(aes(color = p.adjust), size = 8) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_wrap(~ ONTOLOGY, scales = "free_y", nrow = 3, strip.position = "right") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey80"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    plot.title = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),     
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),   
    legend.title = element_text(size = 14)
  ) +
  labs(
    title = "GO Terms Up-regulated Genes KO22",
    x = "",
    y = "",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Up-regulated Genes KO22.jpg"), plot = Gene_ontology_KO22_UP, width = 12, height = 8, dpi = 800)

# 下调基因的点图
Gene_ontology_KO22_DOWN <- ggplot(sorted_ego_down_df_KO22, aes(x = Count, y = Description)) +
  geom_point(aes(color = p.adjust), size = 8) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_wrap(~ ONTOLOGY, scales = "free_y", nrow = 3, strip.position = "right") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey80"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),     
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),   
    legend.title = element_text(size = 14)
  ) +
  labs(
    title = "GO Terms for Down-regulated Genes KO22",
    x = "",
    y = "",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Down-regulated Genes KO22.jpg"), plot = Gene_ontology_KO22_DOWN, width = 12, height = 8, dpi = 800)


# FIH23KO 确保筛选时排除掉 NA 值--------------------------------------------------------------------
up_genes_FIH23KO <- rownames(results_plenti_vs_KO23[!is.na(results_plenti_vs_KO23$padj) & results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange > 0, ])
down_genes_FIH23KO <- rownames(results_plenti_vs_KO23[!is.na(results_plenti_vs_KO23$padj) & results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange < 0, ])
background_genes <- rownames(experimental_data)

# 上调基因的GO分析
ego_up_KO23 <- enrichGO(gene         = up_genes_FIH23KO,
                        OrgDb        = org.Hs.eg.db,
                        keyType      = "SYMBOL",
                        ont          = "ALL",  
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)

ego_up_df_KO23 <- as.data.frame(ego_up_KO23)

# 下调基因的GO分析
ego_down_KO23 <- enrichGO(gene         = down_genes_FIH23KO,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",
                          ont          = "ALL",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

ego_down_df_KO23 <- as.data.frame(ego_down_KO23)

# 上调基因表达的GO分析结果处理

sorted_up_bp_KO23 <- ego_up_df_KO23 %>%
  filter(ONTOLOGY == "BP") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_up_cc_KO23 <- ego_up_df_KO23 %>%
  filter(ONTOLOGY == "CC") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_up_mf_KO23 <- ego_up_df_KO23 %>%
  filter(ONTOLOGY == "MF") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_ego_up_df_KO23 <- bind_rows(sorted_up_bp_KO23, sorted_up_cc_KO23, sorted_up_mf_KO23)

# 下调基因表达的GO分析结果处理，与上调处理相似

sorted_down_bp_KO23 <- ego_down_df_KO23 %>%
  filter(ONTOLOGY == "BP") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_down_cc_KO23 <- ego_down_df_KO23 %>%
  filter(ONTOLOGY == "CC") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_down_mf_KO23 <- ego_down_df_KO23 %>%
  filter(ONTOLOGY == "MF") %>%
  arrange(p.adjust, desc(GeneRatio), desc(Count)) %>%
  head(10)

sorted_ego_down_df_KO23 <- bind_rows(sorted_down_bp_KO23, sorted_down_cc_KO23, sorted_down_mf_KO23)


# 上调基因的点图
Gene_ontology_KO23_UP <- ggplot(sorted_ego_up_df_KO23, aes(x = Count, y = Description)) +
  geom_point(aes(color = p.adjust), size = 8) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_wrap(~ ONTOLOGY, scales = "free_y", nrow = 3, strip.position = "right") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey80"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", fill=NA, linewidth =1),
    plot.title = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),     
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),   
    legend.title = element_text(size = 14)
  ) +
  labs(
    title = "GO Terms Up-regulated Genes KO23",
    x = "",
    y = "",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Up-regulated Genes KO23.jpg"), plot = Gene_ontology_KO23_UP, width = 12, height = 8, dpi = 800)

# 下调基因的点图
Gene_ontology_KO23_DOWN <- ggplot(sorted_ego_down_df_KO23, aes(x = Count, y = Description)) +
  geom_point(aes(color = p.adjust), size = 8) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_wrap(~ ONTOLOGY, scales = "free_y", nrow = 3, strip.position = "right") +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey80"),
    strip.placement = "outside",
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.title = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),     
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 12),   
    legend.title = element_text(size = 14)
  ) +
  labs(
    title = "GO Terms for Down-regulated Genes KO23",
    x = "",
    y = "",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Down-regulated Genes KO23.jpg"), plot = Gene_ontology_KO23_DOWN, width = 12, height = 8, dpi = 800)


#KEGG Enrichment Analysis=======================================================
##Function to Convert Gene Symbols to Entrez IDs--------------------------------
convert_to_entrez <- function(gene_list) {
  entrez_ids <- mapIds(org.Hs.eg.db, 
                       keys = gene_list, 
                       column = "ENTREZID", 
                       keytype = "SYMBOL", 
                       multiVals = "first")
  
  # Remove NA values
  return(entrez_ids[!is.na(entrez_ids)])
}

##Function to Perform KEGG Enrichment-------------------------------------------
perform_kegg_enrichment <- function(gene_list) {
  entrez_ids <- convert_to_entrez(gene_list)  # Reuse the convert_to_entrez() function
  enrichKEGG(gene = entrez_ids,
             organism = "hsa",  # Human KEGG pathways
             pvalueCutoff = 0.05,
             keyType = "kegg")}

##Function to Plot KEGG Results with ggplot2------------------------------------
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
    ggsave(file.path(plot_dir, filename), plot = p, width = 10, height = 6)
  } else {
    print(paste("No significant KEGG pathways found for", title))
  }
}

##Creating .csv and plots-------------------------------------------------------
#Convert Gene Symbols to Entrez IDs
up_entrez_KO22 <- convert_to_entrez(upregulated_genes_KO22)
down_entrez_KO22 <- convert_to_entrez(downregulated_genes_KO22)

#Perform KEGG Enrichment for Upregulated and Downregulated Genes
kegg_up_KO22 <- enrichKEGG(gene = up_entrez_KO22, organism = "hsa", pvalueCutoff = 0.05)
kegg_down_KO22 <- enrichKEGG(gene = down_entrez_KO22, organism = "hsa", pvalueCutoff = 0.05)

#Plot and Save KEGG Results
plot_kegg_ggplot(kegg_up_KO22, "KEGG Pathway for KO22 - Upregulated Genes", "KEGG_Upregulated_KO22.png")
plot_kegg_ggplot(kegg_down_KO22, "KEGG Pathway for KO22 - Downregulated Genes", "KEGG_Downregulated_KO22.png")

#Save Results
write.csv(as.data.frame(kegg_up_KO22@result), file.path(data_dir, "KEGG_Upregulated_KO22.csv"))
write.csv(as.data.frame(kegg_down_KO22@result), file.path(data_dir, "KEGG_Downregulated_KO22.csv"))

#Convert Gene Symbols to Entrez IDs
up_entrez_KO23 <- convert_to_entrez(upregulated_genes_KO23)
down_entrez_KO23 <- convert_to_entrez(downregulated_genes_KO23)

#Perform KEGG Enrichment for Upregulated and Downregulated Genes
kegg_up_KO23 <- enrichKEGG(gene = up_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)
kegg_down_KO23 <- enrichKEGG(gene = down_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)

#Plot and Save KEGG Results
plot_kegg_ggplot(kegg_up_KO23, "KEGG Pathway for KO23 - Upregulated Genes", "KEGG_Upregulated_KO23.png")
plot_kegg_ggplot(kegg_down_KO23, "KEGG Pathway for KO23 - Downregulated Genes", "KEGG_Downregulated_KO23.png")

#Save Results
write.csv(as.data.frame(kegg_up_KO23@result), file.path(data_dir,"KEGG_Upregulated_KO23.csv"))
write.csv(as.data.frame(kegg_down_KO23@result), file.path(data_dir,"KEGG_Downregulated_KO23.csv"))



