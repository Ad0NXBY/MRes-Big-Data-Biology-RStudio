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


