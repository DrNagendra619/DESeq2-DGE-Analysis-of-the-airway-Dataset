### -------------------------------------------
### DESeq2 DGE Analysis (airway/Dexamethasone)
### -------------------------------------------

## -------------------------------------------
## Step 0: Install and Load Packages
## -------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(
  c("DESeq2", "pheatmap", "EnhancedVolcano", "airway", "apeglm"),
  update = FALSE, ask = FALSE
)

# Load libraries
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(airway)
library(apeglm)
library(ggplot2)

## -------------------------------------------
## Step 1: CONFIGURATION - Set Output Path
## -------------------------------------------
# Define the default folder where all results will be saved.
# We use forward slashes "/" for R paths, even on Windows.
output_path <- "D:/DOWNLOADS"

# Create the directory if it doesn't already exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
  print(paste("Created output directory:", output_path))
} else {
  print(paste("Output directory already exists:", output_path))
}


## -------------------------------------------
## Step 2: Load and Prepare Data
## -------------------------------------------
# Load the example dataset
data("airway")
countData <- assay(airway)
colData   <- as.data.frame(colData(airway))

# Create the DESeqDataSet object from the matrix
# The design formula ~dex tells DESeq2 to compare samples based on the "dex" column
dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ dex
)

# Pre-filter genes with very low counts (optional but recommended)
dds <- dds[rowSums(counts(dds)) > 10, ]

## -------------------------------------------
## Step 3: Run DESeq2 Analysis
## -------------------------------------------
# Run the main DESeq2 function
dds <- DESeq(dds)

# Check the names of the comparisons (coefficients)
print(resultsNames(dds))

# Automatically select the coefficient for the "dex" comparison
coef_name <- resultsNames(dds)[grep("dex", resultsNames(dds))]
cat("Using coefficient:", coef_name, "\n")

# Apply Log-Fold Change (LFC) shrinkage for better visualization and ranking
# 'apeglm' is a recommended method for shrinking LFC estimates
resLFC <- lfcShrink(dds, coef = coef_name, type = "apeglm")

## -------------------------------------------
## Step 4: Show Top Results
## -------------------------------------------
# Display the top genes, ordered by adjusted p-value (padj)
print("--- Top 10 Most Significant Genes ---")
head(resLFC[order(resLFC$padj), ], 10)

## -------------------------------------------
## Step 5: Volcano Plot (and Save)
## -------------------------------------------
# Create a volcano plot to visualize p-value vs. log2FoldChange
volcano_plot <- EnhancedVolcano(
  resLFC,
  lab = rownames(resLFC),
  x = "log2FoldChange",
  y = "pvalue",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Volcano Plot: Dexamethasone vs Untreated"
)

# Save the plot
ggsave(
  filename = file.path(output_path, "Volcano_Plot_Airway.png"),
  plot = volcano_plot,
  width = 10,
  height = 8
)

## -------------------------------------------
## Step 6: Heatmap of Top 20 DE Genes (and Save)
## -------------------------------------------
# Transform the count data for visualization (e.g., rlog or VST)
rld <- rlog(dds)

# Get the names of the top 20 most significant genes
top20 <- rownames(resLFC[order(resLFC$padj), ])[1:20]

# Define the file path for the heatmap
heatmap_file <- file.path(output_path, "Heatmap_Top20_Airway.png")

# Create a heatmap of the top 20 genes and save it to the file
pheatmap(
  assay(rld)[top20, ],
  scale = "row", # Scale expression per-gene (row)
  main = "Top 20 Differentially Expressed Genes",
  filename = heatmap_file,
  width = 8,
  height = 10
)

## -------------------------------------------
## Step 7: PCA Plot (and Save)
## -------------------------------------------
# Generate a Principal Component Analysis (PCA) plot
# This shows sample-to-sample distances
pca_plot <- plotPCA(rld, intgroup = "dex") +
  theme_bw() + # Add a cleaner theme
  ggtitle("PCA Plot: Airway Dataset")

# Save the plot
ggsave(
  filename = file.path(output_path, "PCA_Plot_Airway.png"),
  plot = pca_plot,
  width = 8,
  height = 6
)

## -------------------------------------------
## Step 8: Save Results
## -------------------------------------------
# Define the full file path for the results CSV
results_csv_file <- file.path(output_path, "DESeq2_results_airway.csv")

# Save the full results table to a CSV file
write.csv(as.data.frame(resLFC), results_csv_file)

message("---")
message(paste("✅ Done! All results and plots saved to:", output_path))
message("---")
