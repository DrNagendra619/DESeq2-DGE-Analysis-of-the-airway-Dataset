# DESeq2-DGE-Analysis-of-the-airway-Dataset
DESeq2 DGE Analysis of the 'airway' Dataset
# 🔬 DESeq2 DGE Pipeline: Airway Cells (Dexamethasone Treatment)

This R script automates a comprehensive and standard bioinformatics workflow for **Differential Gene Expression (DGE)** analysis using the built-in **`airway`** RNA-Seq dataset. The analysis compares gene expression in human airway smooth muscle cells treated with **Dexamethasone** against **Untreated** controls.

The pipeline utilizes the robust **`DESeq2`** package, implements **`apeglm`** for accurate Log-Fold Change (LFC) estimation, and generates a suite of essential Quality Control (QC) and results visualizations.

## 🚀 Key Features

* **Standard DGE Workflow:** Demonstrates the core, best-practice steps for analyzing RNA-Seq count data using `DESeq2`.
* **Built-in Data Use:** Utilizes the readily available `airway` Bioconductor dataset for immediate, runnable execution.
* **LFC Shrinkage:** Applies **`apeglm`** to generate stabilized Log-Fold Change estimates, improving the reliability of gene ranking and visualization.
* **QC Transformation:** Uses **Regularized Log (rlog)** transformation for accurate distance and clustering plots.
* **Integrated Visualization:** Generates three essential plots: **Volcano Plot**, **PCA Plot**, and **Heatmap** of the top DE genes.
* **Reproducible Output:** Saves the full, shrunken DGE results to a CSV file.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose | |
| :--- | :--- | :--- | :--- |
| **Dataset** | `airway` (Built-in) | RNA-Seq data on Dexamethasone treatment in human airway cells. |
| **DGE Tool** | `DESeq2` | Statistical method optimized for count data. |
| **LFC Estimation** | `apeglm` | Provides robust and biologically meaningful LFCs. |
| **Comparison** | Dexamethasone vs. Untreated | Identifies the transcriptional response to the corticosteroid. |
| **QC Transformation** | **Regularized Log (rlog)** | Prepares data for global QC plots (PCA, Heatmaps). |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

The script automatically installs and loads the necessary Bioconductor and CRAN packages:
* `DESeq2`
* `pheatmap`
* `EnhancedVolcano`
* `airway`
* `apeglm`
* `ggplot2`

### ⚙️ Execution

1.  **Download** the `DESeq2 DGE Analysis of the 'airway' Dataset.R` file.
2.  **Optional:** The output path is set to `D:/DOWNLOADS` (Step 1). You can change this path if needed.
3.  **Execute** the script in your R environment:
    ```R
    source("DESeq2 DGE Analysis of the 'airway' Dataset.R")
    ```

---

## 📁 Output Files (3 Plots + 1 CSV)

All output files are saved to the specified `output_path` (default: `D:/DOWNLOADS`).

### Statistical Results

| Filename | Type | Description |
| :--- | :--- | :--- |
| `DESeq2_results_airway.csv` | CSV | Full table of DGE results, including **shrunken log2FoldChange** (LFC) and **adjusted p-value** (padj). |

### Visualization and QC Plots

| Filename | Analysis Stage | Description |
| :--- | :--- | :--- |
| `Volcano_Plot_Airway.png` | Results | **Volcano Plot** showing shrunken $\log_2 \text{Fold Change}$ vs. $P_{\text{value}}$, highlighting DEGs. |
| `PCA_Plot_Airway.png` | QC / Results | **Principal Component Analysis (PCA)** plot demonstrating strong sample clustering by Dexamethasone treatment. |
| `Heatmap_Top20_Airway.png` | Results | **Heatmap** visualizing the expression profiles of the **Top 20 Most Significant DE Genes** across all samples. |
