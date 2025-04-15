# Install the VennDiagram package if not installed
if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram")
}

# Load the package
library(VennDiagram)

DEG_KO22 <- results_plenti_vs_KO22$GeneSymbol
DEG_KO23 <- results_plenti_vs_KO23$GeneSymbol

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list(KO22_vs_plenti = DEG_KO22, KO23_vs_plenti = DEG_KO23),
  category.names = c("KO22 vs Plenti", "KO23 vs Plenti"),
  filename = NULL,  # Set to NULL to return the plot as an object
  output = TRUE,
  fill = c("blue", "green"),
  alpha = 0.5,
  cex = 1.5,
  fontface = "bold",
  cat.cex = 1.2,
  cat.fontface = "bold",
  cat.col = c("blue", "green")
)

# Display the Venn diagram
grid.draw(venn.plot)
