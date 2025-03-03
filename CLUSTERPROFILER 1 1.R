#ORA enrichment clusterprofiler -ORA - Gene Ontology, KEGG, REACTOM
####Load required libraries#####

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)


library(dplyr)
library(tidyverse)
library(DESeq2)
library(msigdbr)
library(GSEABase)
library(msigdbr)



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi"))
BiocManager::install("clusterProfiler")
library(AnnotationDbi)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(data.table)
library(DOSE)

require(magrittr)
require(DOSE)

setwd("/Users/macbook/Library/CloudStorage/OneDrive-UniversityofSouthampton/Documents/WorkDir/Yihua_Wang/thesis data/yomna_RNA_Seq/ATII-KO-Seq-raw/deseq2_output/Analysis_trim vs un trim DEGS/Untrimmed/NoTCGA/DeSeq2_ouput/CLUSTERPROFILER")
####SET UP DATA ######

res<- read.csv("/Users/macbook/Library/CloudStorage/OneDrive-UniversityofSouthampton/Documents/WorkDir/Yihua_Wang/thesis data/yomna_RNA_Seq/ATII-KO-Seq-raw/deseq2_output/Analysis_trim vs un trim DEGS/Untrimmed/NoTCGA/DeSeq2_ouput/Notrim_Nofilter_ATII_KO_vs_WT_CON_2gps_DEGs.csv")
#res <- read.csv("/Users/macbook/Desktop/yomna/SoCoBio_DTP/Yihua_Wang/thesis data/yomna_RNA_Seq/Dec23/SUBSETS/WT-VS-KO/comparingvs hypooxia data/common_rows_A549_up.csv")
#res <- read.csv ("/Users/macbook/Desktop/yomna/SoCoBio_DTP/Yihua_Wang/thesis data/yomna_RNA_Seq/Dec23/SUBSETS/WT-VS-KO/comparingvs hypooxia data/unique_deg_A549.csv")
#res <- read.csv("/Users/macbook/Desktop/yomna/SoCoBio_DTP/Yihua_Wang/thesis data/yomna_RNA_Seq/Hypoxia ATII/Hypoxia_AT2_16h 2.csv", header = TRUE, stringsAsFactors = FALSE)

# Select FPKM, q-value and log2FoldChange thresholds
select.FPKM <- res$baseMean > 3
select.qval <- res$padj < 0.05
#head(res$padj[res$padj < 0.05])
select.qval <- res$padj < 0.05
summary(select.qval)
select.qval <- res$padj < 0.05
summary(select.qval)

# Filtering differentially expressed genes

# Get upregulated genes
select.log2FC.up <- res$log2FoldChange > 0
select.vec.up <- select.FPKM & select.log2FC.up & select.qval
up.degs.list1 <- as.character(na.omit(res$gene_name[select.vec.up]))
#up.degs.list1.ordered = sort(up.degs.list1, decreasing = TRUE)
up.degs.list1.df <- as.data.frame(up.degs.list1)

#table(select.log2FC.up)
#ta

#up.degs.list1.df <- as.data.frame(up.degs.list1)
#head(up.degs.list1)
#up.degs.list1.ordered = sort(up.degs.list1, decreasing = TRUE)


# Get downregulated genes
select.log2FC.down <- res$log2FoldChange < 0
select.vec.down <- select.FPKM & select.log2FC.down & select.qval
down.degs.list1 <- as.character(na.omit(res$gene_name[select.vec.down]))
#down.degs.list1.ordered = sort(down.degs.list1, decreasing = TRUE)
down.degs.list1.df <- as.data.frame(down.degs.list1)

#table(select.log2FC.down)
#table(select.vec.down)
#down.degs.list1.df <- as.data.frame(down.degs.list1)
#head(down.degs.list1)
#down.degs.list1.ordered = sort(down.degs.list1, decreasing = TRUE)

####GO enrichment analysis#######

# Upregulated genes
enrich.go.up <- enrichGO(gene = up.degs.list1,
                            OrgDb = org.Hs.eg.db,
                            keyType = "SYMBOL",
                            ont = "ALL",
                        pAdjustMethod = "none",
                            pvalueCutoff = 0.1,
                            qvalueCutoff = 0.25)


# Downregulated genes
enrich.go.down <- enrichGO(gene = down.degs.list1,
                              OrgDb = org.Hs.eg.db,
                              keyType = "SYMBOL",
                              ont = "ALL",
                              pAdjustMethod = "none",
                           pvalueCutoff = 0.1,
                           qvalueCutoff = 0.25
                              )


#extract results
result.enrich.go.up <- enrich.go.up@result
result.enrich.go.down <- enrich.go.down@result

#######save results######
write.table(result.enrich.go.up, "erich.go.KO.up.pathways.txt", sep = "\t", col.names = NA)
write.csv(result.enrich.go.up, "erich.go.KO.up.csv", row.names = FALSE)

write.table(result.enrich.go.down.CC, "erich.go.KO.down.pathways.txt", sep = "\t", col.names = NA)
write.csv(result.enrich.go.down, "erich.go.KO.down.csv", row.names = FALSE)



####GO visualization - set up ########
#Important : If you check your result, you would notice GeneRatio values are in Fractional format. To avoid them being represented as fractional format on the plot, you need to transform them into Decimal formats.

#### Function to convert fractions to decimals 
fraction_to_decimal <- function(fraction_string) {
  parts <- unlist(strsplit(fraction_string, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
}

#apply fraction to decimal function to whatever df you want to visualize
result.enrich.go.up$GeneRatio <- sapply(result.enrich.go.up$GeneRatio, fraction_to_decimal)
result.enrich.go.down$GeneRatio <- sapply(result.enrich.go.down$GeneRatio, fraction_to_decimal)

####### Create a new object that contains only rows where the 'Description' column contains the word "extracellular matrix"#######
result_enrich_filtered <- result.enrich.go.up[grepl("extracellular matrix", result.enrich.go.up$Description, ignore.case = TRUE), ]

#you could filter your resultt based on p.adj
filtered_up <- result.enrich.go.up[result.enrich.go.up$p.adjust < 0.05, ]
filtered_down <- result.enrich.go.down.CC[result.enrich.go.down.CC$p.adjust < 0.05, ]
#you can reoder your filtered df based on Count
filtered_up_ordered <- filtered_up [order(filtered_up$GeneRatio, decreasing = T), ]
filtered_down_ordered <- filtered_down[order(filtered_down$Count, decreasing = T), ] 

# Prepare top 50 enriched pathways for up and down adjust this to whatever number of pathways you are interested to visualize
top_25_up <- head(filtered_up, 50)
top_50_down <- head(filtered_down, 50)

####GO plotting ######

#plotting up or down pathways in separatte plots
### A. visualization using ggplot2 package
#dotplot

down_GO_CC_dotplot <- ggplot(top_25_up, aes(x = GeneRatio, y = reorder(Description, Count), color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "#CC0033", high = "#0072B2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(margin = margin(10, 10, 0, 0), size = 8),  # Adjust the right margin for x-axis labels
    axis.text.y = element_text(margin = margin(10, 10, 0, 0), size = 10))   # Adjust the right margin for y-axis labels

down_GO_CC_dotplot
ggsave("top15_down_MF_dotplot.png", width = 10, height = 6)


#barplot
down_GO_CC_barplot <- ggplot(top_15_down, aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#CC0033", high = "#0072B2") +  # Customize colors
  labs(x = "Count", y = "Description", fill = "p.adjust") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(margin = margin(10, 10, 0, 0), size = 8),  # Adjust the right margin for x-axis labels
    axis.text.y = element_text(margin = margin(10, 10, 0, 0), size = 10))   # Adjust the right margin for y-axis labels

down_GO_CC_barplot
ggsave("top15_down_MF_barplot.png", width = 10, height = 6)


#
####

################### You can combine up and down pathways in one object and visualize them
top_15_up <- head(filtered_up, 15)
top_15_down <- head(result.enrich.go.down.CC, 15)

# Add Direction column
top_15_up$Direction <- "Up"
top_15_down$Direction <- "Down"

# Combine up and down data
top_15_combined <- rbind(top_15_up, top_15_down)

# Dotplot
ggplot(top_15_combined, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  facet_wrap(~Direction, scales = "free_y") +
  scale_color_gradient(low = "#CC0033", high = "#0072B2") +
  theme_minimal() +
  labs(title = "GO Enrichment Analysis")
#save plot
ggsave("Combined_GO_MF_dotplot.png", width = 15, height = 9)

dotplot(top_15_combined, x="NES",decreasing = TRUE, showCategory=25, split=".sign", font.size = 8, label_format = 75) +
  scale_x_continuous(limit = c(-3, 3), breaks=seq(-3,3,1)) + facet_grid(.~.sign)+
  theme( 
    legend.text=element_text(size=9),
    legend.position = "right",
    legend.box = "vertical",
    text = element_text(size=12))

# Barplot
ggplot(top_15_combined, aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Direction, scales = "free_y") +
  scale_fill_gradient(low = "#CC0033", high = "#0072B2") +
  theme_minimal() +
  labs(title = "GO Enrichment Analysis")
# Save plot
ggsave("Combined_GO_barplot.png", width = 20, height = 12)



##test
# Assuming 'filtered_up_ordered' is your data frame and it has the columns 'Description', 'GeneRatio', 'ONTOLOGY'

# Filter top 15 items for each ONTOLOGY category
top_items_per_ontology <- filtered_up_ordered %>%
  group_by(ONTOLOGY) %>%
  top_n(10, GeneRatio) %>%
  ungroup()
# Define custom colors for each GO category


custom_colors <- c("BP" = "#F3C566",  # Blue for Biological Process
                   "CC" = "#025959",  # Green for Cellular Component
                   "MF" = "#F06060")  # Orange-Red for Molecular Function

# Generate the bar plot with the filtered top items
top_GO <- ggplot(top_items_per_ontology, aes(x = reorder(Description, GeneRatio), y = GeneRatio, fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +  scale_x_discrete(labels = function(x) str_wrap(x, width = 60)) + # Wrap labels at 30 characters
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  labs(title = "GO Enrichment Analysis")+  
  coord_flip()  # Flip the axes if needed

top_GO
ggsave("top_go.png", width = 8, height = 12)
# If 'Direction' is not a column in your data frame, remove the facet_wrap line

####

ggplot(top_items_per_ontology, aes(x = reorder(Description, -GeneRatio), y = GeneRatio, fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = custom_colors) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 30)) + # Wrap labels
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Adjust text angle and alignment
  labs(title = "GO Enrichment Analysis", x = NULL, y = "Gene Ratio") +
  coord_flip()  # Flip the axes
####
### A .#Visualization using ClusterProfiler embeded visualizationn method

up_GO <- dotplot(enrich.go.up.CC, showCategory = 15, color = "p.adjust", font.size = 6, label_format = 10)
up_GO
ggsave("up_GO_plot.png", up_GO, width = 8, height = 6, units = "in", dpi = 300)

bar_plot.up <- barplot(enrich.go.up.CC, showCategory = 50, , color = "p.adjust", font.size = 6, label_format = 10)
bar_plot.up

down_GO <- dotplot(enrich.go.down.CC, showCategory = 15, color = "p.adjust", font.size = 12, label_format = 10)
down_GO
bar_plot.down <- barplot(enrich.go.down.CC, showCategory = 50, , color = "p.adjust", font.size = 8, width = 0.8, label_format = 10) 
bar_plot.down

#Write results
ggsave("up_GO_plot.png", width = 20, height = 12)
####THE END #####################



#### Function to convert fractions to decimals 
fraction_to_decimal <- function(fraction_string) {
  parts <- unlist(strsplit(fraction_string, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
}

#apply fraction to decimal function to whatever df you want to visualize
result.enrich.go.up$GeneRatio <- sapply(result.enrich.go.up$GeneRatio, fraction_to_decimal)
result.enrich.go.down$GeneRatio <- sapply(result.enrich.go.down$GeneRatio, fraction_to_decimal)

# Convert the GO enrichment result into a data frame for easy customization
up_go_df <- as.data.frame(result.enrich.go.up)
go_top30 <- up_go_df %>%
  group_by(ONTOLOGY) %>%
  slice_max(order_by = -pvalue, n = 15) %>%
  ungroup()

go_top30 <- up_go_df %>%
  group_by(ONTOLOGY) %>%
  slice_max(order_by = Count, n = 15) %>%
  ungroup()

down_go_df <- as.data.frame(result.enrich.go.down)
go_top30 <- down_go_df %>%
  group_by(ONTOLOGY) %>%
  slice_max(order_by = Count, n = 15) %>%
  ungroup()

library(stringr)
library(ggplot2)
# Create the plot
########
p <- ggplot(go_top30, aes(x = GeneRatio, y = reorder(Description_wrapped, -p.adjust), size = Count, color = -p.adjust)) +
  geom_point(alpha = 0.8) +  # Set the transparency of the points
  scale_size(range = c(1, 5)) +  # Adjust bubble size
  scale_color_gradient(low = "blue", high = "red") +  # Set the color gradient
  facet_wrap(~ONTOLOGY, scales = "free_x", nrow = 1) +  # Facet by BP, CC, MF
  labs(x = "Gene Ratio", y = NULL, color = "-log10(p.adjust)", size = "Count") +  # Set axis labels
  scale_x_continuous(breaks = seq(0.02, 0.1, by = 0.5)) +  # Increase space between x-axis numbers
  theme_minimal() +  # Use minimal theme
  theme(
    panel.grid.major = element_line(color = "gray85"),  # Light grid lines
    panel.border = element_rect(color = "black", fill = NA),  # Black border
    strip.text = element_text(size = 12, face = "bold"),  # Bold titles for BP, CC, MF
    axis.text.y = element_text(size = 10),  # Adjust text size for GO terms
    plot.title = element_text(size = 14, hjust = 0.5),  # Centered title
    legend.position = "right"  # Place legend on the right
  ) + 
  ggtitle("A")  # Add "A" as a title in the top left corner
########
# Install ggtext if it's not already installed
if (!require(ggtext)) {
  install.packages("ggtext")
}
library(ggtext)
install.packages("stringr")
library(stringr)
# Wrap the GO term descriptions to a specific width (e.g., 30 characters)
go_top30$Description_wrapped <- str_wrap(go_top30$Description, width = 20)

# Create the plot
p <- ggplot(go_top30, aes(x = GeneRatio, y = reorder(Description_wrapped, -1*p.adjust), size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +  # Set the transparency of the points
  scale_size(range = c(1, 5)) +  # Adjust bubble size
  scale_color_gradient(low = "blue", high = "red") +  # Set the color gradient
  facet_wrap(~ONTOLOGY, scales = "free_x", nrow = 1) +  # Facet by BP, CC, MF
  labs(x = "Gene Ratio", y = NULL, color = "p.adjust", size = "Count") +  # Set axis labels
  theme_minimal() +  # Use minimal theme
  theme(
    panel.spacing = unit(2, "lines"),  # Increase the space between facets (adjust the number as needed)
    panel.grid.major = element_line(color = "gray85"),  # Light grid lines
    panel.border = element_rect(color = "black", fill = NA),  # Black border
    strip.text = element_text(size = 12, face = "bold"),  # Bold titles for BP, CC, MF
    axis.text.y = ggtext::element_markdown(size = 12),  # Adjust text size for GO terms, markdown for HTML rendering
    plot.title = element_text(size = 14, hjust = 0.5),  # Centered title
    legend.position = "right"  # Place legend on the right
  ) +  
  scale_y_discrete(labels = function(labels) {
    # Wrap the labels to limit the line length to a certain width (e.g., 40 characters)
    wrapped_labels <- str_wrap(labels, width = 40)
    
    # Extract qvalues for each label
    qvalue <- go_top30$qvalue[match(labels, go_top30$Description_wrapped)]
    
    # Color labels based on qvalue
    ifelse(qvalue < 0.1,
           paste0("<span style='color:#000000;'>", labels, "</span>"),  # Black for qvalue < 0.1
           paste0("<span style='color:#999999;'>", labels, "</span>"))# Grey otherwise
  }) +
  
  ggtitle("GO for UPregulated genes")  # Add "A" as a title in the top left corner

#Render the plot with special handling for HTML formatting
p + theme(axis.text.y = element_markdown(size = 10))  # Use element_markdown to apply the red color
print(p)

# Save the plot with more reasonable dimensions
ggsave("GO_DOWN_DEGS.png", height = 10, width = 12, dpi = 300)  # Save with 300 dpi for high quality


######TREEMAP######
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GOplot")

library(treemapify)    # For treemap visualization

go_updf <- as.data.frame(result.enrich.go.up)
go_up_top50 <- go_updf %>%
  group_by(ONTOLOGY) %>%
  slice_max(order_by = Count, n = 15) %>%
  ungroup()

go_downdf <- as.data.frame(result.enrich.go.down)
go_down_top50 <- go_downdf %>%
  group_by(ONTOLOGY) %>%
  slice_max(order_by = Count, n = 15) %>%
  ungroup()

# Combine the top 50 results for both up and downregulated DEGs
go_combined <- rbind(
  cbind(go_up_top50, regulation = "Upregulated"),
  cbind(go_down_top50, regulation = "Downregulated")
)



# Plotting the treemap
ggplot(go_combined, aes(area = Count, fill = -log10(p.adjust), label = Description, subgroup = regulation)) +
  geom_treemap() +
  geom_treemap_text(colour = "white", place = "centre", grow = TRUE, reflow = TRUE) +
  facet_wrap(~regulation, nrow = 2) +  # Separate upregulated and downregulated terms
  scale_fill_continuous(low =  "#084594", high = "#B30000") +  # Darker purple to dark red gradient
  # Adjust color scale
  theme_minimal() +
  labs(title = "Top 50 Biological Process GO Terms Enriched in DEGs",
       fill = "-log10(p.adjust)")


# Create the heatmap with a color gradient similar to your example
ggplot(go_combined, aes(area = Count, fill = p.adjust, label = Description, subgroup = regulation)) +
  geom_treemap() +
  geom_treemap_text(colour = "white", place = "centre", grow = TRUE, reflow = TRUE) +
  facet_wrap(~regulation, nrow = 2, scales = "free") +  # Facet to create the two sections
  scale_fill_gradient2(low = "#FFD700", mid = "#6A3D9A", high = "#1F78B4", midpoint = 0) +  # Purple to deep blue gradient
  theme_minimal() +
  labs(title = "Top 50 Biological Process GO Terms Enriched in DEGs",
       fill = "p.adjust") +
  theme(strip.text = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

####KEGG - ORA ######################
setwd ('/Users/macbook/Desktop/yomna/SoCoBio_DTP/Yihua_Wang/thesis data/yomna_RNA_Seq/Dec23/SUBSETS/WT-VS-KO/ENRICHMENT/ora/KEGG')

#####1st method - using enricher function - General KEGG enrichment using enricher##
kegg=read.gmt('/Users/macbook/Desktop/yomna/SoCoBio_DTP/Yihua_Wang/thesis data/yomna_RNA_Seq/DEGs/yomna-degs/A-DEG-results/KO_CON_VS_WT_CON/enrichement - wt-vs-ko/ClusterProfiler - KEGG/c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt')
orgdb <- org.Hs.eg.db
resKEGG=enricher(gene = down.degs.list1, TERM2GENE=kegg)
resKEGG_down=resKEGG@result

resKEGG_up=enricher(gene = up.degs.list1, TERM2GENE=kegg)
resKEGG_up=resKEGG_up@result

write.csv(resKEGG_down,"ericher.KEGG.downreg.csv")
write.csv(resKEGG_up,"ericher.KEGG.upreg.csv")

#######2nd method ## Usinng enrichKEGG which is more tailored for KEGG analysis######
orgdb <- org.Hs.eg.db
search_kegg_organism('hsa', by='kegg_code')
library(AnnotationDbi)  # For mapIds function

#mapping the gene symbols in gene list to their corresponding Entrez Gene IDs using the org.Hs.eg.db annotation package.
DEG.entrez_id_up = mapIds(x = org.Hs.eg.db,
                       keys =  up.degs.list1,
                       keytype = "SYMBOL",
                       column = "ENTREZID")

DEG.entrez_id_down = mapIds(x = org.Hs.eg.db,
                       keys =  down.degs.list1,
                       keytype = "SYMBOL",
                       column = "ENTREZID") 
# KEGG enrichment analysis
enrich.kegg.up <- enrichKEGG(gene = DEG.entrez_id_up,
                               organism = 'hsa', # for human, use 'hsa'
                               pAdjustMethod = 'BH', # method for adjusting p-values
                               qvalueCutoff = 0.05, # threshold for q-value
                               pvalueCutoff = 0.05) # threshold for p-value


enrich.kegg.down <- enrichKEGG(gene = DEG.entrez_id_down,
                               organism = 'hsa', # for human, use 'hsa'
                               pAdjustMethod = 'BH', # method for adjusting p-values
                               qvalueCutoff = 0.05, # threshold for q-value
                               pvalueCutoff = 0.05) # threshold for p-value

result.enrich.kegg.up <- enrich.kegg.up@result
result.enrich.kegg.down <- enrich.kegg.down@result

# Function to convert Entrez IDs to gene symbols in the geneID column
convert_entrez_to_symbols <- function(entrez_ids) {
  # Split the geneID column (Entrez IDs separated by "/")
  entrez_list <- strsplit(entrez_ids, "/")[[1]]
  
  # Map Entrez IDs to gene symbols
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = entrez_list,
                         keytype = "ENTREZID",
                         column = "SYMBOL",
                         multiVals = "first")
  
  # Combine the gene symbols into a single string, separated by "/"
  return(paste(gene_symbols, collapse = "/"))
}

# Apply the function to the geneID column in the KEGG enrichment results for upregulated genes
result.enrich.kegg.up$geneID_symbols <- sapply(result.enrich.kegg.up$geneID, convert_entrez_to_symbols)

# Apply the function to the geneID column in the KEGG enrichment results for downregulated genes
result.enrich.kegg.down$geneID_symbols <- sapply(result.enrich.kegg.down$geneID, convert_entrez_to_symbols)

# Function to convert fractions to decimals
fraction_to_decimal <- function(fraction_string) {
  parts <- unlist(strsplit(fraction_string, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
}

###apply fraction to decimal function to whatever df you want to visualize
result.enrich.kegg.up$GeneRatio <- sapply(result.enrich.kegg.up$GeneRatio, fraction_to_decimal)
result.enrich.kegg.down$GeneRatio <- sapply(result.enrich.kegg.down$GeneRatio, fraction_to_decimal)

write.csv(result.enrich.kegg.up,"enrich.kegg.up.csv")
write.csv(result.enrich.kegg.down,"enrich.kegg.down.csv")

KEGG_UP_df <- as.data.frame(result.enrich.kegg.up)
KEGG_down_df <- as.data.frame(result.enrich.kegg.down)

go_top30 <- KEGG_UP_df %>%
  slice_max(order_by = GeneRatio, n = 20) %>%
  ungroup()
kegg_top30_down <- KEGG_down_df %>%
  slice_max(order_by = GeneRatio, n = 20) %>%
  ungroup()

library(stringr)
# Install ggtext if it's not already installed
if (!require(ggtext)) {
  install.packages("ggtext")
}
library(ggtext)

# Wrap the GO term descriptions to a specific width (e.g., 30 characters)
go_top30$Description_wrapped <- str_wrap(go_top30$Description, width = 60)
kegg_top20_down$Description_wrapped <- str_wrap(kegg_top20_down$Description, width = 60)

# Create the plot
p <- ggplot(go_top30, aes(x = GeneRatio, y = reorder(Description_wrapped, -log10(p.adjust)), size = Count, color = -log10(p.adjust))) +
  geom_point(alpha = 0.8) +  # Set the transparency of the points
  scale_size(range = c(1, 5)) +  # Adjust bubble size
  scale_color_gradient(low = "blue", high = "red") +  # Set the color gradient
  #facet_wrap(~ONTOLOGY, scales = "free_x", nrow = 1) +  # Facet by BP, CC, MF
  labs(x = "Gene Ratio", y = NULL, color = "-log10(p.adjust)", size = "Count") +  # Set axis labels
  theme_minimal() +  # Use minimal theme
  theme(
    panel.grid.major = element_line(color = "gray85"),  # Light grid lines
    panel.border = element_rect(color = "black", fill = NA),  # Black border
    strip.text = element_text(size = 12, face = "bold"),  # Bold titles for BP, CC, MF
    axis.text.y = element_text(size = 10),  # Adjust text size for GO terms
    plot.title = element_text(size = 14, hjust = 0.5),  # Centered title
    legend.position = "right"  # Place legend on the right
  ) + 
  scale_y_discrete(labels = function(labels) {
    # Extract qvalues from the kegg_top20_down data
    qvalue <- go_top30$qvalue[match(labels, go_top30$Description_wrapped)]
    
    # Color labels based on qvalue
    ifelse(qvalue < 0.1,
           paste0("<span style='color:#000000;'>", labels, "</span>"),  # Black for qvalue < 0.1
           paste0("<span style='color:#999999;'>", labels, "</span>"))  # Grey otherwise
  }) +
  ggtitle("KEGG Plot for UPregulated Genes")


 # scale_y_discrete(labels = function(labels) {
    # Check and highlight the "cellular senescence" term in red
    #ifelse(labels %in% c("Cellular senescence", "HIF-1 signaling pathway", "Glycolysis / Gluconeogenesis", "Calcium signaling pathway"),
  #         paste0("<span style='color:red;'>", labels, "</span>"), 
   #        labels)
 # })+
  ggtitle("A")  # Add "A" as a title in the top left corner

# Render the plot with special handling for HTML formatting
p + theme(axis.text.y = element_markdown(size = 10))  # Use element_markdown to apply the red color
#print(p)
ggsave("KEGG_UP.png",height = 10, width = 12, dpi = 300 )



#you could filter your resultt based on p.adj
filtered_up <- result.enrich.kegg.up[result.enrich.kegg.up$p.adjust < 0.05, ]
filtered_down <- result.enrich.kegg.down[result.enrich.kegg.down$p.adjust < 0.05, ]
#you can reoder your filtered df based on Count
filtered_up_ordered <- filtered_up [order(filtered_up$Count, decreasing = T), ]
filtered_down_ordered <- filtered_down[order(filtered_down$Count, decreasing = T), ] 


####KEGG visualization #######################
library(ggplot2) 
################### You can combine up and down pathways in one object and visualize them
top_15_up <- head(result.enrich.kegg.up, 15)
top_15_down <- head(result.enrich.kegg.down, 15)

# Add Direction column
top_15_up$Direction <- "Up"
top_15_down$Direction <- "Down"

# Combine up and down data
top_15_combined <- rbind(top_15_up, top_15_down)

# Dotplot
kegg <- ggplot(top_15_combined, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  facet_wrap(~Direction, scales = "free_y") +
  scale_color_gradient(low = "#CC0033", high = "#0072B2") +
  theme_minimal() +
  labs(title = "KEGG Enrichment Analysis")

kegg
#save plot
ggsave("Combined_KEGG_dotplot.png", width = 15, height = 9)

# Barplot
kegg_bar <- ggplot(top_15_combined, aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Direction, scales = "free_y") +
  scale_fill_gradient(low = "#CC0033", high = "#0072B2") +
  theme_minimal() +
  labs(title = "KEGG Enrichment Analysis")

kegg_bar
# Save plot
ggsave("Combined_KEGG_barplot.png", width = 20, height = 12)


##Plot your data
#dotplot

up_GO_dotplot <- ggplot(top_15_up, aes(x = GeneRatio, y = reorder(Description, Count), color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "#CC0033", high = "#0072B2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(margin = margin(10, 10, 0, 0), size = 8),  # Adjust the right margin for x-axis labels
    axis.text.y = element_text(margin = margin(10, 10, 0, 0), size = 10))   # Adjust the right margin for y-axis labels

up_GO_dotplot
ggsave("KEGG_up_dotplot.png")
library(ggplot2)
library(stringr) # for str_wrap()

# Your original plot code here
# ...

# Add this to wrap labels at 30 characters
up_KEGG_dotplot <- up_GO_dotplot +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 40))

# Print the plot
up_KEGG_dotplot
#barplot
bar_plot_up_go <- ggplot(top_15_up, aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#CC0033", high = "#0072B2") +  # Customize colors
  labs(x = "Count", y = "Description", fill = "p.adjust") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(margin = margin(10, 10, 0, 0), size = 8),  # Adjust the right margin for x-axis labels
    axis.text.y = element_text(margin = margin(10, 10, 0, 0), size = 10))   # Adjust the right margin for y-axis labels

bar_plot_up_go
ggsave("KEGG_up_barplot.png")

####THE END #####################
####REACTOME - ORA  ###################

BiocManager::install("ReactomePA")
library(ReactomePA)

setwd('/Users/macbook/Desktop/yomna/SoCoBio_DTP/Yihua_Wang/thesis data/yomna_RNA_Seq/Dec23/SUBSETS/WT-VS-KO/ENRICHMENT/ora/REACTOM')

orgdb <- org.Hs.eg.db

#mapping the gene symbols in gene list to their corresponding Entrez Gene IDs using the org.Hs.eg.db annotation package.
DEG.entrez_id_up = mapIds(x = org.Hs.eg.db,
                          keys =  up.degs.list1,
                          keytype = "SYMBOL",
                          column = "ENTREZID")

DEG.entrez_id_down = mapIds(x = org.Hs.eg.db,
                            keys =  down.degs.list1,
                            keytype = "SYMBOL",
                            column = "ENTREZID") 
# Perform enrichment analysis with REACTOME pathways
result_reactome_up <- enrichPathway(DEG.entrez_id_up, organism = "human",
                                 pvalueCutoff = 0.05, pAdjustMethod = "BH",
                                 qvalueCutoff = 0.25, minGSSize = 10, maxGSSize = 500, readable = FALSE)

result_reactome_down <- enrichPathway(DEG.entrez_id_down, organism = "human",
                                 pvalueCutoff = 0.05, pAdjustMethod = "BH",
                                 qvalueCutoff = 0.25, minGSSize = 10, maxGSSize = 500, readable = FALSE)

result_reactome_up=result_reactome_up@result
result_reactome_down=result_reactome_down@result

write.csv(result_reactome_up,"result_reactome_up.csv")
write.csv(result_reactome_down,"result_reactome_down.csv")

# Function to convert fractions to decimals
fraction_to_decimal <- function(fraction_string) {
  parts <- unlist(strsplit(fraction_string, "/"))
  as.numeric(parts[1]) / as.numeric(parts[2])
}

###apply fraction to decimal function to whatever df you want to visualize
result_reactome_up$GeneRatio <- sapply(result_reactome_up$GeneRatio, fraction_to_decimal)
result_reactome_down$GeneRatio <- sapply(result_reactome_down$GeneRatio, fraction_to_decimal)

# Function to convert Entrez IDs to gene symbols in the geneID column
convert_entrez_to_symbols <- function(entrez_ids) {
  # Split the geneID column (Entrez IDs separated by "/")
  entrez_list <- strsplit(entrez_ids, "/")[[1]]
  
  # Map Entrez IDs to gene symbols
  gene_symbols <- mapIds(org.Hs.eg.db,
                         keys = entrez_list,
                         keytype = "ENTREZID",
                         column = "SYMBOL",
                         multiVals = "first")
  
  # Combine the gene symbols into a single string, separated by "/"
  return(paste(gene_symbols, collapse = "/"))
}

# Apply the function to the geneID column in the KEGG enrichment results for upregulated genes
result_reactome_up$geneID_symbols <- sapply(result_reactome_up$geneID, convert_entrez_to_symbols)

# Apply the function to the geneID column in the KEGG enrichment results for downregulated genes
result_reactome_down$geneID_symbols <- sapply(result_reactome_down$geneID, convert_entrez_to_symbols)


result_reactome_up_df <- as.data.frame(result_reactome_up)
result_reactome_down_df <- as.data.frame(result_reactome_down)

library(dplyr)

top20_up <- result_reactome_up_df %>%
  slice_max(order_by = Count, n = 20) %>%
  ungroup()
top20_down <- result_reactome_down_df %>%
  slice_max(order_by = Count, n = 20) %>%
  ungroup()
# Function to round p.adjust values and format scientific notation
format_p_adjust <- function(p.adjust) {
  formatted_value <- formatC(p.adjust, format = "e", digits = 2)
  return(formatted_value)
}


# Applying the function to a dataframe column
top20_up$p.adjust <- sapply(top20_up$p.adjust, format_p_adjust)
top20_down$p.adjust <- sapply(top20_down$p.adjust, format_p_adjust)

# Alternatively, with dplyr's mutate
top20_up <- top20_up %>%
  mutate(p.adjust = sapply(p.adjust, format_p_adjust))


library(stringr)
# Install ggtext if it's not already installed
if (!require(ggtext)) {
  install.packages("ggtext")
}
library(ggtext)

# Wrap the GO term descriptions to a specific width (e.g., 30 characters)
top20_up$Description_wrapped <- str_wrap(top20_up$Description, width = 20)
top20_down$Description_wrapped <- str_wrap(top20_down$Description, width = 20)

# Create the plot
p <- ggplot(top20_down, aes(x = GeneRatio, y = reorder(Description_wrapped, p.adjust), size = Count, color = as.numeric(p.adjust))) +
  geom_point(alpha = 0.8) +  # Set the transparency of the points
  scale_size(range = c(1, 5)) +  # Adjust bubble size
  scale_color_gradient(low = "blue", high = "red") +  # Set the color gradient
  #facet_wrap(~ONTOLOGY, scales = "free_x", nrow = 1) +  # Facet by BP, CC, MF
  labs(x = "Gene Ratio", y = NULL, color = "p.adjust", size = "Count") +  # Set axis labels
  theme_minimal() +  # Use minimal theme
  theme(
    panel.grid.major = element_line(color = "gray85"),  # Light grid lines
    panel.border = element_rect(color = "black", fill = NA),  # Black border
    strip.text = element_text(size = 12, face = "bold"),  # Bold titles for BP, CC, MF
    axis.text.y = element_text(size = 10),  # Adjust text size for GO terms
    plot.title = element_text(size = 14, hjust = 0.5),  # Centered title
    legend.position = "right", # Place legend on the right
    plot.margin = unit(c(1, 1, 2, 1), "lines")  # Increase bottom margin for longer labels
  ) +  
  scale_y_discrete(labels = function(labels) {
    # Extract qvalues for each label
    qvalue <- top20_down$qvalue[match(labels, top20_down$Description_wrapped)]
    
    # Color labels based on qvalue
    ifelse(qvalue < 0.1,
           paste0("<span style='color:#000000;'>", labels, "</span>"),  # Black for qvalue < 0.1
           paste0("<span style='color:#999999;'>", labels, "</span>"))  # Grey otherwise
  }) +
  ggtitle("REACTOME plot for Downregulated genes")  # Add "A" as a title in the top left corner

# Render the plot with special handling for HTML formatting
p + theme(axis.text.y = element_markdown(size = 10))  # Use element_markdown to apply the red color

#print(p)
ggsave("REACTOME_DOWN.png", height = 10, width = 12, dpi = 300)

#you could filter your result based on p.adj
filtered_up <- result_reactome_up[result_reactome_up$p.adjust < 0.05, ]
filtered_down <- result_reactome_down[result_reactome_down$p.adjust < 0.05, ]
#you can reoder your filtered df based on Count
filtered_up_ordered <- filtered_up [order(filtered_up$Count, decreasing = T), ]
filtered_down_ordered <- filtered_down[order(filtered_down$Count, decreasing = T), ] 

#####REACTOME visualization #########
#Visualization Using ggPlot2 R package###

top_15_up <- head(result_reactome_up, 15)
top_15_down <- head(result_reactome_down, 15)

# Add Direction column
top_15_up$Direction <- "Up"
top_15_down$Direction <- "Down"

# Combine up and down data
top_15_combined <- rbind(top_15_up, top_15_down)

# Dotplot
combined_reactome <- ggplot(top_15_combined, aes(x = GeneRatio, y = reorder(Description, GeneRatio), color = p.adjust, size = Count)) +
  geom_point() +
  facet_wrap(~Direction, scales = "free_y") +
  scale_color_gradient(low = "#CC0033", high = "#0072B2") +
  theme_minimal() +
  labs(title = "Reactome Enrichment Analysis")
combined_reactome

#save plot
ggsave("Combined_reactome_dotplot.png", width = 20, height = 9)

# Barplot
combined_bar <- ggplot(top_15_combined, aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Direction, scales = "free_y") +
  scale_fill_gradient(low = "#CC0033", high = "#0072B2") +
  theme_minimal() +
  labs(title = "Reactome Enrichment barplot Analysis")
combined_bar

# Save plot
ggsave("Combined_reactome_barplot.png", width = 20, height = 9)

#separate plots

# Dotplot
up_GO_dotplot <- ggplot(top_15_down, aes(x = GeneRatio, y = reorder(Description, Count), color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "#CC0033", high = "#0072B2") +
  theme_minimal() +
  theme()   # Adjust the right margin for y-axis labels

up_GO_dotplot
ggsave("down_reatome_dotplot.png", width = 10, height = 6)


up_reactome_dotplot <- ggplot(top_15_up, aes(x = GeneRatio, y = reorder(Description, Count), color = p.adjust, size = Count)) +
  geom_point() +
  scale_color_gradient(low = "#CC0033", high = "#0072B2") +
  theme_minimal() +
  theme()  

up_reactome_dotplot
ggsave("up_reactome_dotplot.png", width = 10, height = 6)

##test #to wrap long labels
library(ggplot2)
library(stringr) # for str_wrap()

# Your original plot code here
# ...

# Add this to wrap labels at 30 characters
up_reactome_dotplot <- up_reactome_dotplot +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 50))

# Print the plot
up_reactome_dotplot
##
#barplot
down_reactome_barplot <- ggplot(top_15_down, aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#CC0033", high = "#0072B2") +  # Customize colors
  labs(x = "Count", y = "Description", fill = "p.adjust") +
  theme_minimal() +
  theme()   # Adjust the right margin for y-axis labels

down_reactome_barplot
ggsave("down_reactome_barplot.png", width = 10, height = 6)


up_reactome_bar <- ggplot(top_15_up, aes(x = Count, y = reorder(Description, -p.adjust), fill = p.adjust)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Direction, scales = "free_y") +
  scale_fill_gradient(low = "#CC0033", high = "#0072B2") +
  theme_minimal() +
  theme()

up_reactome_bar
# Save plot
ggsave("up_reactome_barplot.png", width = 10, height = 6)


####THE END #####################