plot_dir <- "Analysis_plots with own data and chenyin code"
data_dir <- "List_for_analysis with own data and chenyin code"

dir.create(plot_dir, showWarnings = FALSE)
dir.create(data_dir, showWarnings = FALSE)



# Load necessary packages
library(data.table)
library(DESeq2)
library(edgeR)

# Read files
Control1 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_A1_featurecounts.txt")
Control2 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_A2_featurecounts.txt")
Control3 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_A3_featurecounts.txt")
FIH22KO1 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_B1_featurecounts.txt")
FIH22KO2 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_B2_featurecounts.txt")
FIH22KO3 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_B3_featurecounts.txt")
FIH23KO1 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_C1_featurecounts.txt")
FIH23KO2 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_C2_featurecounts.txt")
FIH23KO3 <- fread("C:/Users/Brandon/Documents/MRes Big Data Biology/Data analysis/ChenY_RNA_seq/Feature counts/Ch_C3_featurecounts.txt")

# Set column names
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

# Merge data
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
zero_counts <- rowSums(combined_data == 0)
combined_data <- combined_data[zero_counts <= 8, ]


##Using Ensembl as the database-------------------------------------------------
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") #sometimes the server is temporarily unavailble, go do something else

#Fetch gene names corresponding to Ensembl IDs
ensemblID_to_genesymbol <- combined_data[,1]
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") 
gene_mapping <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id",
                      values = ensemblID_to_genesymbol,
                      mart = mart)

##Rename columns for merging----------------------------------------------------
colnames(gene_mapping) <- c("GeneID", "GeneSymbol")

#Convert combined_data to dataframe (if it’s still a matrix)
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
combined_data <- combined_data %>% dplyr::select(-GeneID, -GeneSymbol, -FinalName)

# Convert the data to a DGEList object
counts <- as.matrix(combined_data)
group <- factor(c(rep("Control", 3), rep("FIH22KO", 3), rep("FIH23KO", 3)))
dge <- DGEList(counts = counts, group = group)

# Filter low-expressed genes
keep <- filterByExpr(dge, group = group, min.count = 15, min.total.count = 70, min.prop = 0.2)
dge_filtered <- dge[keep,]

# Update the filtered data
filtered_data <- as.data.frame(dge_filtered$counts)

# Convert the data to the format required by DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(filtered_data),
  colData = DataFrame(condition = group),
  design = ~ condition
)

# Run DESeq2 normalization and differential expression analysis
dds <- DESeq(dds)

# View normalized counts
normalized_counts <- assay(rlog(dds))

# Test for normal distribution

apply(normalized_counts, 2, function(column_data) {
  ks.test(column_data, "pnorm", mean = mean(column_data), sd = sd(column_data))
})

# Test for Poisson distribution

# Fit a Poisson distribution for each sample and calculate the likelihood
poisson_fit <- function(column_data) {
  lambda <- mean(column_data)  # The λ for Poisson distribution is the mean of the data
  log_likelihood <- sum(dpois(column_data, lambda, log = TRUE))  # Calculate log-likelihood of Poisson distribution
  return(log_likelihood)
}

# Perform the fit for each sample
log_likelihood_results <- apply(normalized_counts, 2, poisson_fit)

# View results
log_likelihood_results

# Save results
write.csv(normalized_counts, file.path(data_dir, "normalized_counts.csv"))



#DEG======================================================================================================
# 加载包 & 读取数据-------------------------------------------------------------
library(limma)
library(ggplot2)
library(ggrepel)
library(VennDiagram)
library(dplyr)
library(scales)
library(ReactomePA)

# 读取已经标准化的数据
Ch_data <- read.csv("C:/Users/Brandon/Documents/MRes Big Data Biology/R studio/List_for_analysis with own data and chenyin code/normalized_counts.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
rownames(Ch_data) <- Ch_data$X
Ch_data <- Ch_data[, -1]


# Create group information
group <- factor(c(rep("Control", 3), rep("FIH22KO", 3), rep("FIH23KO", 3)))

# Filter low-expressed genes using limma's filterByExpr()
keep_genes <- filterByExpr(Ch_data, 
                           group = group, 
                           min.count = 15, 
                           min.total.count = 70, 
                           min.prop = 0.2)

# Subset the data (correct way)
Ch_data_filtered <- Ch_data[keep_genes, ]


# DEG --------------------------------------------------------------------------

# 创建分组信息
group <- factor(c(rep("Control", 3), rep("FIH22KO", 3), rep("FIH23KO", 3)))

# 创建 DESeq2 对象
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(Ch_data),
  colData = data.frame(condition = group),
  design = ~ condition
)

dds <- DESeq(dds)

# 提取 FIH22KO vs Control 的差异表达结果
results_FIH22KO <- results(dds, contrast=c("condition", "FIH22KO", "Control"))

# 提取 FIH23KO vs Control 的差异表达结果
results_FIH23KO <- results(dds, contrast=c("condition", "FIH23KO", "Control"))

# For FIH22KO vs Control:
results_FIH22KO_df <- as.data.frame(results_FIH22KO)
results_FIH22KO_df$gene <- rownames(results_FIH22KO_df)  # Add gene names as a column
# Save the full DEG results as a CSV file
write.csv(results_FIH22KO_df, file.path(data_dir, "DEG_FIH22KO_vs_Control.csv"), row.names = FALSE)

# Create a ranked list file for GSEA (using log2FoldChange as the ranking metric)
rnk_FIH22KO <- results_FIH22KO_df[!is.na(results_FIH22KO_df$padj) & results_FIH22KO_df$padj < 0.05, c("gene", "log2FoldChange")]
rnk_FIH22KO <- rnk_FIH22KO[order(rnk_FIH22KO$log2FoldChange, decreasing = TRUE), ]
write.table(rnk_FIH22KO, file.path(data_dir, "DEG_FIH22KO_vs_Control.rnk"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


# For FIH23KO vs Control:
results_FIH23KO_df <- as.data.frame(results_FIH23KO)
results_FIH23KO_df$gene <- rownames(results_FIH23KO_df)  # Add gene names as a column
# Save the full DEG results as a CSV file
write.csv(results_FIH23KO_df, file.path(data_dir, "DEG_FIH23KO_vs_Control.csv"), row.names = FALSE)

# Create a ranked list file for GSEA (using log2FoldChange as the ranking metric)
rnk_FIH23KO <- results_FIH23KO_df[!is.na(results_FIH23KO_df$padj) & results_FIH23KO_df$padj < 0.05, c("gene", "log2FoldChange")]
rnk_FIH23KO <- rnk_FIH23KO[order(rnk_FIH23KO$log2FoldChange, decreasing = TRUE), ]
write.table(rnk_FIH23KO, file.path(data_dir, "DEG_FIH23KO_vs_Control.rnk"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


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

# 提取显著上调和下调基因取交集，绘制韦恩图--------------------------------------


up_genes_FIH22KO <- rownames(results_FIH22KO[!is.na(results_FIH22KO$padj) & results_FIH22KO$padj < 0.05 & results_FIH22KO$log2FoldChange > 0, ])
down_genes_FIH22KO <- rownames(results_FIH22KO[!is.na(results_FIH22KO$padj) & results_FIH22KO$padj < 0.05 & results_FIH22KO$log2FoldChange < 0, ])

up_genes_FIH23KO <- rownames(results_FIH23KO[!is.na(results_FIH23KO$padj) & results_FIH23KO$padj < 0.05 & results_FIH23KO$log2FoldChange > 0, ])
down_genes_FIH23KO <- rownames(results_FIH23KO[!is.na(results_FIH23KO$padj) & results_FIH23KO$padj < 0.05 & results_FIH23KO$log2FoldChange < 0, ])

gene_lists <- list(
  "FIH22KO_up" = up_genes_FIH22KO,
  "FIH22KO_down" = down_genes_FIH22KO,
  "FIH23KO_up" = up_genes_FIH23KO,
  "FIH23KO_down" = down_genes_FIH23KO
)

#VENN PLOT=================================================================================
venn.plot <- venn.diagram(
  x = gene_lists,
  category.names = c("FIH KO 1Up Regulated", "FIH KO 1Down Regulated", "FIH KO 2Up Regulated", "FIH KO 2Down Regulated"),
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


# VOLCANO PLOT==================================================================================

# 准备 FIH22KO 的数据，并排除 NA 值
res_df_KO22 <- results_FIH22KO
colnames(res_df_KO22)[colnames(res_df_KO22) == "logFC"] <- "log2FoldChange"
colnames(res_df_KO22)[colnames(res_df_KO22) == "adj.P.Val"] <- "padj"
res_df_KO22$gene_name <- rownames(res_df_KO22)

# 移除包含 NA 值的行
res_df_KO22 <- na.omit(res_df_KO22)

# 添加标记字段
res_df_KO22$threshold <- with(res_df_KO22, ifelse(log2FoldChange > 1 & padj < 0.05, "up-regulated",
                                                  ifelse(log2FoldChange < -1 & padj < 0.05, "down-regulated", 
                                                         "not significant")))

# 创建火山图
Volcano_plot_KO22 <- ggplot(res_df_KO22, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
  geom_point(alpha=1, size=2) +
  scale_color_manual(values = c("up-regulated" = "red", "down-regulated" = "blue", "not significant" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha=0.5) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", alpha=0.5) + 
  geom_vline(xintercept = -1, linetype = "dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha=0.5) +
  geom_text_repel(
    data = subset(res_df_KO22, padj < 0.05 & abs(log2FoldChange) > 1),
    aes(label = gene_name),
    size = 4,
    colour = "black",
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey80"),
    strip.placement = "outside",
    plot.title = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),     
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),   
    legend.title = element_text(size = 14)
  ) +
labs(title = "FIH22KO vs Control Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value")

# 打印火山图
ggsave(file.path(plot_dir, "volcano_KO22.jpg"), plot = Volcano_plot_KO22, width = 12, height = 8, dpi = 800)


# 准备 FIH23KO 的数据，并排除 NA 值
res_df_KO23 <- results_FIH23KO
colnames(res_df_KO23)[colnames(res_df_KO23) == "logFC"] <- "log2FoldChange"
colnames(res_df_KO23)[colnames(res_df_KO23) == "adj.P.Val"] <- "padj"
res_df_KO23$gene_name <- rownames(res_df_KO23)

# 移除包含 NA 值的行
res_df_KO23 <- na.omit(res_df_KO23)

# 添加标记字段
res_df_KO23$threshold <- with(res_df_KO23, ifelse(log2FoldChange > 1 & padj < 0.05, "up-regulated",
                                                  ifelse(log2FoldChange < -1 & padj < 0.05, "down-regulated", 
                                                         "not significant")))

# 创建火山图
Volcano_plot_KO23 <- ggplot(res_df_KO23, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
  geom_point(alpha=1, size=2) +
  scale_color_manual(values = c("up-regulated" = "red", "down-regulated" = "blue", "not significant" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha=0.5) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", alpha=0.5) + 
  geom_vline(xintercept = -1, linetype = "dashed", color = "black", alpha=0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha=0.5) +
  geom_text_repel(
    data = subset(res_df_KO23, padj < 0.05 & abs(log2FoldChange) > 1),
    aes(label = gene_name),
    size = 4,
    colour = "black",
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.5, "lines")
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    strip.background = element_rect(fill = "grey80"),
    strip.placement = "outside",
    plot.title = element_text(size = 18), 
    axis.title.x = element_text(size = 16), 
    axis.title.y = element_text(size = 16),     
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    legend.text = element_text(size = 14),   
    legend.title = element_text(size = 14)
  ) +
  labs(title = "FIH23KO vs Control Volcano Plot", x = "Log2 Fold Change", y = "-Log10 p-value")

# 打印火山图
ggsave(file.path(plot_dir, "volcano_KO23.jpg"), plot = Volcano_plot_KO23, width = 12, height = 8, dpi = 800)


#GENE ONTOLOGY=======================================================================================
# FIH22KO 确保筛选时排除掉 NA 值--------------------------------------------------------------------
up_genes_FIH22KO <- rownames(results_FIH22KO[!is.na(results_FIH22KO$padj) & results_FIH22KO$padj < 0.05 & results_FIH22KO$log2FoldChange > 0, ])
down_genes_FIH22KO <- rownames(results_FIH22KO[!is.na(results_FIH22KO$padj) & results_FIH22KO$padj < 0.05 & results_FIH22KO$log2FoldChange < 0, ])
background_genes <- rownames(Ch_data)

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

# FIH23KO 确保筛选时排除掉 NA 值-----------------------------------------------
up_genes_FIH23KO <- rownames(results_FIH23KO[!is.na(results_FIH23KO$padj) & results_FIH23KO$padj < 0.05 & results_FIH23KO$log2FoldChange > 0, ])
down_genes_FIH23KO <- rownames(results_FIH23KO[!is.na(results_FIH23KO$padj) & results_FIH23KO$padj < 0.05 & results_FIH23KO$log2FoldChange < 0, ])
background_genes <- rownames(Ch_data)

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
