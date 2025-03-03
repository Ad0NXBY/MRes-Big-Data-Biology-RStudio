#DESEQ2======================================================================

#Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = Counts_gene[,3:11],
                              colData = data.frame(condition),
                              design = ~condition)

#DESeq2 analysis
dds <- DESeq(dds)
                  
#Transform the data for PCA
vsd <- vst(dds, blind = FALSE) #variance stabilizing transformation

#SCREE & PCA PLOT ====================================================================
##Scree plot --------------------------------------------------------------------------------
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
#DESeq2 results for comparisons
#plenti vs KO22
results_plenti_vs_KO22 <- results(dds, contrast = c("condition", "KO22", "plenti"))
summary(results_plenti_vs_KO22) #summary of results
#plenti vs KO23
results_plenti_vs_KO23 <- results(dds, contrast = c("condition", "KO23", "plenti"))
summary(results_plenti_vs_KO23) #summary of results
