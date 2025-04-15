#PREPROCESSING==================================================================
##Installing the necessary packages---------------------------------------------
install.packages("ggplot2")       
install.packages("tidyverse")      
install.packages("factoextra")     
install.packages("ggpubr")         
install.packages("data.table")     
install.packages("limma")          
install.packages("ggrepel")        
install.packages("VennDiagram")    
install.packages("scales")         

if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("DESeq2")           
BiocManager::install("biomaRt")           
BiocManager::install("edgeR")             
BiocManager::install("clusterProfiler")   
BiocManager::install("org.Hs.eg.db")      
BiocManager::install("AnnotationDbi")     
BiocManager::install("GSVA")             
BiocManager::install("msigdbr")          
BiocManager::install("GSEABase")          
BiocManager::install("ComplexHeatmap")   


##Create output directories-----------------------------------------------------
plot_dir <- "Analysis plots with combined code"
data_dir <- "List_for_analysis with combined code"

dir.create(plot_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)

##Load necessary libraries------------------------------------------------------
library(data.table)  # For fread()
library(dplyr)       # For data manipulation (select, left_join, distinct, etc.)
library(biomaRt)     # For mapping Ensembl IDs to gene symbols

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

##Filter genes with zero counts across all samples------------------------------ (USING ROW SUMS, REMOVES ROW IF THE ENTIRE ROW = ZERO)
zero_counts <- rowSums(combined_data == 0)
combined_data <- combined_data[zero_counts <= 8, ]

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

# Keep only the first occurrence of each gene symbol
combined_data <- combined_data %>%
  distinct(FinalName, .keep_all = TRUE)

# Set row names and remove unnecessary columns
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
##Load necessary libraries------------------------------------------------------
library(DESeq2)
library(limma)
library(dplyr)

##Reading in the normalized counts file----------------------------------------
experimental_data <- read.csv("C:/Users/Brandon/Documents/MRes Big Data Biology/R studio/List_for_analysis with combined code//normalized_counts.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
rownames(experimental_data) <- experimental_data$X
experimental_data <- experimental_data[,-1]

##Filter low-expressed genes using limma's filterByExpr()----------------------
keep_genes <- filterByExpr(experimental_data, 
                           group = condition, 
                           min.count = 15, 
                           min.total.count = 70, 
                           min.prop = 0.2)

##Subset the data (correct way)-------------------------------------------------
experimental_data_filtered <- experimental_data[keep_genes, ]

condition <- factor(c(rep("plenti", 3), rep("KO22", 3), rep("KO23", 3)))


##Creating dds object-----------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(experimental_data),
  colData = data.frame(condition = condition),
  design = ~ condition
)

dds <- DESeq(dds)

##Creating .csv DEG list for FIHKO22--------------------------------------------
results_plenti_vs_KO22 <- results(dds, contrast = c("condition", "KO22", "plenti"))
results_plenti_vs_KO22_df <- as.data.frame(results_plenti_vs_KO22)
results_plenti_vs_KO22_df$gene <- rownames(results_plenti_vs_KO22_df)

#Save full DEG results
write.csv(results_plenti_vs_KO22_df, 
          file.path(data_dir, "DEG_KO22_vs_Plenti_filterByExpr.csv"), 
          row.names = TRUE)

##Filter for significant DEGs for KO22 GSEA (padj < 0.05)-----------------------
GSEA_FIH22KO <- results_plenti_vs_KO22_df %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(desc(log2FoldChange))

#Save the filtered DEG list
write.csv(GSEA_FIH22KO,
          file.path(data_dir, "GSEA_plenti_vs_KO22.csv"),
          row.names = TRUE)

##Prepare ranked list for GSEA--------------------------------------------------
GSEA_ranked_list_KO22 <- GSEA_FIH22KO %>%
  dplyr::select(gene, log2FoldChange)

#Write .rnk file (no headers or row names as required by GSEA)
write.table(
  GSEA_ranked_list_KO22,
  file.path(data_dir, "GSEA_ranked_list_KO22_filter_by_padj_lessthan0.05.rnk"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)



#提取 FIH23KO vs Control 的差异表达结果
results_plenti_vs_KO23 <- results(dds, contrast = c("condition", "KO23", "plenti"))
results_plenti_vs_KO23_df <- as.data.frame(results_plenti_vs_KO23)
results_plenti_vs_KO23_df$gene <- rownames(results_plenti_vs_KO23_df)

#Save full DEG results
write.csv(results_plenti_vs_KO23_df, 
          file.path(data_dir, "DEG_KO23_vs_Plenti_filterByExpr.csv"), 
          row.names = TRUE)

##Filter for significant DEGs for KO23 GSEA(padj < 0.05)------------------------
GSEA_FIH23KO <- results_plenti_vs_KO23_df %>%
  dplyr::filter(padj < 0.05) %>%
  dplyr::arrange(desc(log2FoldChange))

#Save the filtered DEG list
write.csv(GSEA_FIH23KO,
          file.path(data_dir, "GSEA_plenti_vs_KO23.csv"),
          row.names = TRUE)

##Prepare ranked list for GSEA--------------------------------------------------
GSEA_ranked_list_KO23 <- GSEA_FIH23KO %>%
  dplyr::select(gene, log2FoldChange)

#Write .rnk file (no headers or row names as required by GSEA)
write.table(
  GSEA_ranked_list_KO23,
  file.path(data_dir, "GSEA_ranked_list_KO23_filter_by_padj_lessthan0.05.rnk"),
  sep = "\t",
  col.names = FALSE,
  row.names = FALSE,
  quote = FALSE
)

#Scree & PCA Plot===========================================================================
#Load necessary libraries------------------------------------------------------
library(DESeq2)
library(factoextra)
library(ggpubr)
library(ggplot2)

##Transform the data for PCA----------------------------------------------------
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
##Load necessary libraries------------------------------------------------------
library(VennDiagram)
library(grid)
library(ggplot2)

##Filtering out genes and creating a gene list----------------------------------
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

##Creating Venn plot------------------------------------------------------------
venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = c("KO22 Upregulated", "KO22 Downregulated", "KO23 Upregulated", "KO23 Downregulated"),
  filename = NULL,
  imagetype = "jpg",  
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
  fontfamily = "Arial",  
  fontface = "bold",
  cat.col = c("darkred", "darkgreen", "darkblue", "purple4"),
  cat.cex = 1.3,
  cat.fontfamily = "Arial", 
  cat.default.pos = "outer",
  cat.pos = c(-40, 40, -30, 30),
  cat.dist = c(0.3, 0.3, 0.2, 0.2),
  cat.fontface = "bold",
  rotation.degree = 0,
  margin = 0.2
)

grid.draw(venn.plot)

ggsave(file.path(plot_dir, "vennplot.jpg"), plot = venn.plot, width = 12, height = 8, dpi = 800)


#VOLCANO PLOT===================================================================
##Load necessary libraries------------------------------------------------------
library(ggplot2)
library(ggrepel)
        
##Volcano plot of Plent vs KO22-------------------------------------------------
#Ensure gene names are set as rownames
volcano_data_22 <- as.data.frame(results_plenti_vs_KO22)
volcano_data_22$Gene <- rownames(volcano_data_22)

#Remove NA rows
volcano_data_22 <- na.omit(volcano_data_22)

#Define significance thresholds
volcano_data_22$Significance <- "Not Significant"
volcano_data_22$Significance[volcano_data_22$padj < 0.05 & volcano_data_22$log2FoldChange > log2(2)] <- "Upregulated"
volcano_data_22$Significance[volcano_data_22$padj < 0.05 & volcano_data_22$log2FoldChange < -log2(2)] <- "Downregulated"

#Create volcano plot with geom_text_repel
Volcano_Plot_Plenti_v_KO22 <- ggplot(volcano_data_22, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
  geom_vline(xintercept = c(-log2(2), 0, log2(2)), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(
    data = subset(volcano_data_22, padj < 0.05 & abs(log2FoldChange) > log2(2)),
    aes(label = Gene),
    size = 3,
    color = "black",
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = 25
  ) +
  xlim(-10, 10) +
  labs(title = "Volcano Plot (Plenti vs KO22)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal()

#Save plot
ggsave(file.path(plot_dir, "Volcano_Plot_Plenti_vs_KO22.jpg"), plot = Volcano_Plot_Plenti_v_KO22,
       width = 10, height = 6, dpi = 800)
##Volcano plot of Plent vs KO23-------------------------------------------------
#Ensure gene names are set as rownames
volcano_data_23 <- as.data.frame(results_plenti_vs_KO23)
volcano_data_23$Gene <- rownames(volcano_data_23)

#Remove NA rows
volcano_data_23 <- na.omit(volcano_data_23)

#Define significance thresholds
volcano_data_23$Significance <- "Not Significant"
volcano_data_23$Significance[volcano_data_23$padj < 0.05 & volcano_data_23$log2FoldChange > log2(2)] <- "Upregulated"
volcano_data_23$Significance[volcano_data_23$padj < 0.05 & volcano_data_23$log2FoldChange < -log2(2)] <- "Downregulated"

#Create volcano plot with geom_text_repel
Volcano_Plot_Plenti_v_KO23 <- ggplot(volcano_data_23, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
  geom_vline(xintercept = c(-log2(2), 0, log2(2)), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
  geom_text_repel(
    data = subset(volcano_data_23, padj < 0.05 & abs(log2FoldChange) > log2(2)),
    aes(label = Gene),
    size = 3,
    color = "black",
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines"),
    max.overlaps = 25
  ) +
  xlim(-10, 10) +
  labs(title = "Volcano Plot (Plenti vs KO23)",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme_minimal()

#Save plot
ggsave(file.path(plot_dir, "Volcano_Plot_Plenti_vs_KO23.jpg"), plot = Volcano_Plot_Plenti_v_KO23,
       width = 10, height = 6, dpi = 800)

#GENE ONTOLOGY=================================================================
##Load necessary libraries------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tibble)
library(forcats)
library(ggplot2)

##FIH22KO --------------------------------------------------------------------
up_genes_FIH22KO <- rownames(results_plenti_vs_KO22[!is.na(results_plenti_vs_KO22$padj) & results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange > 0, ])
down_genes_FIH22KO <- rownames(results_plenti_vs_KO22[!is.na(results_plenti_vs_KO22$padj) & results_plenti_vs_KO22$padj < 0.05 & results_plenti_vs_KO22$log2FoldChange < 0, ])
background_genes <- rownames(experimental_data)

###Processing of GO analysis results of upregulated gene expression-------------
ego_up_KO22 <- enrichGO(gene         = up_genes_FIH22KO,
                        OrgDb        = org.Hs.eg.db,
                        keyType      = "SYMBOL",
                        ont          = "ALL",  
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)

ego_up_df_KO22 <- as.data.frame(ego_up_KO22)

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


#Sort within each ONTOLOGY group by Count
sorted_ego_up_df_KO22 <- sorted_ego_up_df_KO22 %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(Count), .by_group = TRUE) %>%
  mutate(Description = factor(Description, levels = unique(Description)))


###Creating upregulated Gene ontology graph-----------------------------------
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
    x = "Gene Count",
    y = "GO Term",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Up-regulated Genes KO22.jpg"), plot = Gene_ontology_KO22_UP, width = 12, height = 8, dpi = 800)


###Processing of GO analysis results of downregulated gene expression-------------
ego_down_KO22 <- enrichGO(gene         = down_genes_FIH22KO,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",
                          ont          = "ALL",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

ego_down_df_KO22 <- as.data.frame(ego_down_KO22)

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

#Sort within each ONTOLOGY group by Count
sorted_ego_down_df_KO22 <- sorted_ego_down_df_KO22 %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(Count), .by_group = TRUE) %>%
  mutate(Description = factor(Description, levels = unique(Description)))


###Creating downregulated Gene ontology graph-----------------------------------
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
    x = "Gene Count",
    y = "GO Term",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Down-regulated Genes KO22.jpg"), plot = Gene_ontology_KO22_DOWN, width = 12, height = 8, dpi = 800)


##FIH23KO------------------------------------------------------------------------
up_genes_FIH23KO <- rownames(results_plenti_vs_KO23[!is.na(results_plenti_vs_KO23$padj) & results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange > 0, ])
down_genes_FIH23KO <- rownames(results_plenti_vs_KO23[!is.na(results_plenti_vs_KO23$padj) & results_plenti_vs_KO23$padj < 0.05 & results_plenti_vs_KO23$log2FoldChange < 0, ])
background_genes <- rownames(experimental_data)

###Processing of GO analysis results of upregulated gene expression-------------
ego_up_KO23 <- enrichGO(gene         = up_genes_FIH23KO,
                        OrgDb        = org.Hs.eg.db,
                        keyType      = "SYMBOL",
                        ont          = "ALL",  
                        pAdjustMethod = "BH",
                        qvalueCutoff = 0.05)

ego_up_df_KO23 <- as.data.frame(ego_up_KO23)

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

#Sort within each ONTOLOGY group by Count
sorted_ego_up_df_KO23 <- sorted_ego_up_df_KO23 %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(Count), .by_group = TRUE) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

###Creating upregulated Gene ontology graph-----------------------------------
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
    x = "Gene Count",
    y = "GO Term",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Up-regulated Genes KO23.jpg"), plot = Gene_ontology_KO23_UP, width = 12, height = 8, dpi = 800)

###Processing of GO analysis results of downregulated gene expression-------------
ego_down_KO23 <- enrichGO(gene         = down_genes_FIH23KO,
                          OrgDb        = org.Hs.eg.db,
                          keyType      = "SYMBOL",
                          ont          = "ALL",
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.05)

ego_down_df_KO23 <- as.data.frame(ego_down_KO23)

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

#Sort within each ONTOLOGY group by Count
sorted_ego_down_df_KO23 <- sorted_ego_down_df_KO23 %>%
  group_by(ONTOLOGY) %>%
  arrange(desc(Count), .by_group = TRUE) %>%
  mutate(Description = factor(Description, levels = unique(Description)))


###Creating downregulated Gene ontology graph-----------------------------------
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
    x = "Gene Count",
    y = "GO Term",
    color = "adj.P.Val"
  )

ggsave(file.path(plot_dir, "GO Terms Down-regulated Genes KO23.jpg"), plot = Gene_ontology_KO23_DOWN, width = 12, height = 8, dpi = 800)


#KEGG Enrichment Analysis=======================================================
##Load necessary libraries------------------------------------------------------
library(clusterProfiler)
library(ggplot2)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)

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
up_entrez_KO22 <- convert_to_entrez(up_genes_FIH22KO)
down_entrez_KO22 <- convert_to_entrez(down_genes_FIH22KO)

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
up_entrez_KO23 <- convert_to_entrez(up_genes_FIH23KO)
down_entrez_KO23 <- convert_to_entrez(down_genes_FIH23KO)

#Perform KEGG Enrichment for Upregulated and Downregulated Genes
kegg_up_KO23 <- enrichKEGG(gene = up_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)
kegg_down_KO23 <- enrichKEGG(gene = down_entrez_KO23, organism = "hsa", pvalueCutoff = 0.05)

#Plot and Save KEGG Results
plot_kegg_ggplot(kegg_up_KO23, "KEGG Pathway for KO23 - Upregulated Genes", "KEGG_Upregulated_KO23.png")
plot_kegg_ggplot(kegg_down_KO23, "KEGG Pathway for KO23 - Downregulated Genes", "KEGG_Downregulated_KO23.png")

#Save Results
write.csv(as.data.frame(kegg_up_KO23@result), file.path(data_dir,"KEGG_Upregulated_KO23.csv"))
write.csv(as.data.frame(kegg_down_KO23@result), file.path(data_dir,"KEGG_Downregulated_KO23.csv"))


#GSVA===========================================================================
##Load necessary libraries------------------------------------------------------
library(GSVA)
library(msigdbr)
library(GSEABase)
library(ComplexHeatmap)

##Get human gene sets from msigdbr----------------------------------------------
msigdbr_species <- msigdbr(species = "Homo sapiens")
hallmark_genesets_df <- msigdbr_species[msigdbr_species$gs_cat == "H", ]


##Split gene sets---------------------------------------------------------------
gset.idx.list <- split(hallmark_genesets_df$gene_symbol, hallmark_genesets_df$gs_name)

##Create GeneSet objects with unique genes per set------------------------------
geneSets <- lapply(names(gset.idx.list), function(name) {
  GeneSet(geneIds = unique(gset.idx.list[[name]]), 
          setName = name, 
          geneIdType = SymbolIdentifier())
})

geneSetCollection <- GeneSetCollection(geneSets)

expr_matrix <- assay(vsd)  


##Create parameter object and run GSVA------------------------------------------
param <- gsvaParam(expr = expr_matrix, 
                   geneSets = geneSetCollection, 
                   kcdf = "Poisson")
gsva_results <- gsva(param)


head(gsva_results)

write.csv(gsva_results, file.path(data_dir, "GSVA_results.csv"), row.names = TRUE)



##Calculate t-test statistics and p-values for Control vs KO22------------------
ttest_results <- t(apply(gsva_results, 1, function(row) {
  test <- t.test(row[c("Control1_Count", "Control2_Count", "Control3_Count")], 
                 row[c("FIH22KO1_Count", "FIH22KO2_Count", "FIH22KO3_Count")])
  c(t_value = test$statistic, p_value = test$p.value)
}))

##Clean up gene set names and create data frame--------------------------------
bar_data_KO22 <- data.frame(
  GeneSet = sub("HALLMARK_", "", rownames(ttest_results)),
  T_Value = ttest_results[, "t_value"],
  P_Value = ttest_results[, "p_value"],
  FillColor = ifelse(ttest_results[, "p_value"] < 0.05,
                     ifelse(ttest_results[, "t_value"] > 0, "red", "blue"),
                     "grey")
)

#Sort by t-value for better visual
bar_data_KO22 <- bar_data_KO22[order(bar_data_KO22$T_Value, decreasing = TRUE),]

##Plot bar graph with p-value-based coloring------------------------------------
GSVA_BarPlot_KO22 <- ggplot(bar_data_KO22, aes(x = reorder(GeneSet, T_Value), y = T_Value, fill = FillColor)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  coord_flip() +
  labs(x = "Gene Set", y = "t value of GSVA", title = "GSVA Analysis for Control vs KO22") +
  theme_minimal() +
  theme(text = element_text(color = "black"),
        legend.position = "none")

#Save the plot to file
ggsave(file.path(plot_dir, "GSVA_BarPlot_Plenti_vs_KO22.png"), plot = GSVA_BarPlot_KO22, width = 10, height = 6, dpi = 800)



##Calculate t-test statistics and p-values for Control vs KO23------------------
ttest_results_KO23 <- t(apply(gsva_results, 1, function(row) {
  test <- t.test(row[c("Control1_Count", "Control2_Count", "Control3_Count")], 
                 row[c("FIH23KO1_Count", "FIH23KO2_Count", "FIH23KO3_Count")])
  c(t_value = test$statistic, p_value = test$p.value)
}))

##Clean up gene set names and create data frame for KO23------------------------
bar_data_KO23 <- data.frame(
  GeneSet = sub("HALLMARK_", "", rownames(ttest_results_KO23)),
  T_Value = ttest_results_KO23[, "t_value"],
  P_Value = ttest_results_KO23[, "p_value"],
  FillColor = ifelse(ttest_results_KO23[, "p_value"] < 0.05,
                     ifelse(ttest_results_KO23[, "t_value"] > 0, "red", "blue"),
                     "grey")
)

#Sort by t-value for better visual in KO23
bar_data_KO23 <- bar_data_KO23[order(bar_data_KO23$T_Value, decreasing = TRUE),]

##Plot bar graph for KO23 with p-value-based coloring---------------------------
GSVA_BarPlot_KO23 <- ggplot(bar_data_KO23, aes(x = reorder(GeneSet, T_Value), y = T_Value, fill = FillColor)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  coord_flip() +
  labs(x = "Gene Set", y = "t value of GSVA", title = "GSVA Analysis for Control vs KO23") +
  theme_minimal() +
  theme(text = element_text(color = "black"),
        legend.position = "none")

#Save the plot to file for KO23
ggsave(file.path(plot_dir, "GSVA_BarPlot_Plenti_vs_KO23.png"), plot = GSVA_BarPlot_KO23, width = 10, height = 6, dpi = 800)



##Calculate t-test statistics and p-values for KO22 vs KO23---------------------
#ttest_results_KO22_vs_KO23 <- t(apply(gsva_results, 1, function(row) {
  #test <- t.test(row[c("FIH22KO1_Count", "FIH22KO2_Count", "FIH22KO3_Count")], 
                 row[c("FIH23KO1_Count", "FIH23KO2_Count", "FIH23KO3_Count")])
 # c(t_value = test$statistic, p_value = test$p.value)
#}))

##Clean up gene set names and create data frame for KO22 vs KO23----------------
#bar_data_KO22_vs_KO23 <- data.frame(
  #GeneSet = sub("HALLMARK_", "", rownames(ttest_results_KO22_vs_KO23)),
  #T_Value = ttest_results_KO22_vs_KO23[, "t_value.t"],  # Corrected column name
 # P_Value = ttest_results_KO22_vs_KO23[, "p_value"],     # No change needed
  #FillColor = ifelse(ttest_results_KO22_vs_KO23[, "p_value"] < 0.05,
                    # ifelse(ttest_results_KO22_vs_KO23[, "t_value.t"] > 0, "red", "blue"),
                     #"grey")
)


#Sort by t-value for better visual in KO22 vs KO23
#bar_data_KO22_vs_KO23 <- bar_data_KO22_vs_KO23[order(bar_data_KO22_vs_KO23$T_Value, decreasing = TRUE),]

##Plot bar graph for KO22 vs KO23 with p-value-based coloring-------------------
#GSVA_BarPlot_KO22_vs_KO23 <- ggplot(bar_data_KO22_vs_KO23, aes(x = reorder(GeneSet, T_Value), y = T_Value, fill = FillColor)) +
  #geom_bar(stat = "identity") +
  #scale_fill_identity() +
  #coord_flip() +
  #labs(x = "Gene Set", y = "t value of GSVA", title = "GSVA Analysis for KO22 vs KO23") +
 # theme_minimal() +
 # theme(text = element_text(color = "black"),
       # legend.position = "none")

#Save the plot to file for KO22 vs KO23
#ggsave(file.path(plot_dir, "GSVA_BarPlot_KO22_vs_KO23.png"), plot = GSVA_BarPlot_KO22_vs_KO23, width = 10, height = 6, dpi = 800)



##Extract relevant GSVA scores for all 3 conditions for heat map----------------
gsva_subset <- gsva_results[, c("Control1_Count", "Control2_Count", "Control3_Count",
                                "FIH22KO1_Count", "FIH22KO2_Count", "FIH22KO3_Count",
                                "FIH23KO1_Count", "FIH23KO2_Count", "FIH23KO3_Count")]

#Optionally, remove HALLMARK_ prefix
rownames(gsva_subset) <- sub("HALLMARK_", "", rownames(gsva_subset))

#Optionally scale rows (mean-center each gene set)
gsva_scaled <- t(scale(t(gsva_subset)))

##Create annotation for columns-------------------------------------------------
annotation_col <- data.frame(
  Condition = c(rep("Control", 3), rep("KO22", 3), rep("KO23", 3))
)
rownames(annotation_col) <- colnames(gsva_scaled)

ha <- HeatmapAnnotation(
  Condition = annotation_col$Condition,
  col = list(Condition = c("Control" = "#4DBBD5", "KO22" = "#E64B35", "KO23" = "#00A087"))
)

##Create the heatmap using ComplexHeatmap---------------------------------------
GSVA_heatmap <- Heatmap(gsva_scaled,  # GSVA results matrix
                        name = "GSVA Score",  # Colorbar label
                        top_annotation = ha,  # Annotations for columns (conditions)
                        cluster_rows = TRUE,  # Cluster rows (gene sets)
                        cluster_columns = FALSE,  # Don't cluster columns to preserve order (optional)
                        show_row_names = TRUE,  # Show row (gene set) names
                        show_column_names = TRUE,  # Show column (sample) names
                        col = colorRampPalette(c("navy", "white", "firebrick3"))(100),  # Custom color palette
                        heatmap_legend_param = list(title = "GSVA Score"),  # Color bar title
                        row_names_gp = gpar(fontsize = 7)# Adjust row name font size
)

jpeg(file.path(plot_dir,"GSVA_heatmap.jpg"), width = 1200, height = 1000, res = 150)
draw(GSVA_heatmap)
dev.off()

