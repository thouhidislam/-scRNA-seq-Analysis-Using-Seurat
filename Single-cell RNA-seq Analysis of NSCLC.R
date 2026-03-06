############################################################
# Project: Single-cell RNA-seq Analysis of NSCLC
# Author: Thouhid Islam
# Purpose: Process and analyze NSCLC 10X Genomics data
# Tools: Seurat, tidyverse, hdf5r
# Note: Fully reproducible workflow
############################################################

############################
# 1. Setup Environment
############################

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install necessary packages only if not installed
# install.packages("Seurat")
# install.packages("tidyverse")
# install.packages("hdf5r")

library(Seurat)
library(tidyverse)
library(hdf5r)

############################
# 2. Load Raw Single-Cell Data
############################

# Define file path for the 10X Genomics H5 matrix
data_file <- "C:/......./20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5"


# Load the dataset into R
raw_counts <- Read10X_h5(filename = data_file)

# Extract only the gene expression counts
gene_counts <- raw_counts$`Gene Expression`

############################
# 3. Initialize Seurat Object
############################

# Create a Seurat object for downstream analysis
seurat_obj <- CreateSeuratObject(
  counts = gene_counts,
  project = "NSCLC",
  min.cells = 3,       # Keep genes expressed in at least 3 cells
  min.features = 200   # Keep cells with at least 200 detected genes
)

############################
# 4. Quality Control Metrics
############################

# Calculate mitochondrial gene percentage to check cell quality
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)

# Visualize distributions of key metrics
VlnPlot(
  seurat_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3
)

# Explore relationship between total counts and detected features
FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
) + geom_smooth(method = "lm")

# Filter out low-quality cells based on gene count and mitochondrial content
seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 2500 &
    percent.mt < 5
)

############################
# 5. Normalize the Data
############################

# Apply log-normalization to make gene expression comparable across cells
seurat_obj <- NormalizeData(seurat_obj)

############################
# 6. Identify Highly Variable Genes
############################

# Select 2000 variable genes for dimensionality reduction
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

# Inspect top 10 variable genes
top10_genes <- head(VariableFeatures(seurat_obj), 10)

# Visualize variable genes with labels
plot_var <- VariableFeaturePlot(seurat_obj)
LabelPoints(plot = plot_var, points = top10_genes, repel = TRUE)

############################
# 7. Scale Data for PCA
############################

# Scale the entire dataset before running PCA
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)

############################
# 8. Principal Component Analysis
############################

# Run PCA to reduce dimensionality
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))

# Print top genes contributing to the first 5 PCs
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize PCA heatmap for first dimension
DimHeatmap(seurat_obj, dims = 1, cells = 500, balanced = TRUE)

# Determine optimal number of PCs
ElbowPlot(seurat_obj)

############################
# 9. Clustering Cells
############################

# Build a K-nearest neighbor graph
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)

# Cluster cells at multiple resolutions for flexibility
seurat_obj <- FindClusters(seurat_obj, resolution = c(0.1, 0.3, 0.5, 0.7, 1))

# Use resolution 0.5 as main clustering result
Idents(seurat_obj) <- "RNA_snn_res.0.5"

############################
# 10. Dimensional Reduction and Visualization
############################

# Compute UMAP embedding for visualization
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:15
)

# Visualize clusters on the UMAP projection
DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "RNA_snn_res.0.5",
  label = TRUE
)

############################
# 11. Save Processed Seurat Object
############################

saveRDS(seurat_obj, file = "C:/Users/HP/Desktop/start_seurat/seurat_nsclc_final.rds")

############################
# 12. Record Session Info
############################

# Save package versions and R session info for reproducibility
sessionInfo()
