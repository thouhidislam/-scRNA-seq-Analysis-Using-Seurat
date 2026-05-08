# =========================================================
# Reading Single-Cell RNA-seq Matrices in Multiple Formats
# and Converting Them into Seurat Objects
# =========================================================

# Load libraries
library(Seurat)
library(SeuratDisk)

# ---------------------------------------------------------
# 1. Read .RDS format
# ---------------------------------------------------------

# Load Seurat object stored as .rds
rds_obj <- readRDS('ependymal_cells.rds')

# ---------------------------------------------------------
# 2. Read 10X CellRanger .HDF5 (.h5) format
# ---------------------------------------------------------

# Import 10X Genomics HDF5 matrix
hdf5_obj <- Read10X_h5(
  filename = "20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
  use.names = TRUE,
  unique.features = TRUE
)

# Convert to Seurat object
seurat_hdf5 <- CreateSeuratObject(counts = hdf5_obj)

# ---------------------------------------------------------
# 3. Read .mtx matrix format
# ---------------------------------------------------------

# Import matrix, feature, and barcode files
mtx_obj <- ReadMtx(
  mtx = "raw_feature_bc_matrix/matrix.mtx.gz",
  features = "raw_feature_bc_matrix/features.tsv.gz",
  cells = "raw_feature_bc_matrix/barcodes.tsv.gz"
)

# Convert to Seurat object
seurat_mtx <- CreateSeuratObject(counts = mtx_obj)

# ---------------------------------------------------------
# 4. Read .loom format
# ---------------------------------------------------------

# Connect to loom file
loom_obj <- Connect(
  filename = "adult-hem-organs-10X-bone-marrow.loom",
  mode = 'r'
)

# Convert to Seurat object
seurat_loom <- as.Seurat(loom_obj)

# ---------------------------------------------------------
# 5. Read .h5ad (AnnData) format
# ---------------------------------------------------------

# Step 1: Convert .h5ad to .h5Seurat
Convert(
  "adata_SS2_for_download.h5ad",
  dest = "h5seurat",
  overwrite = TRUE
)

# Step 2: Load .h5Seurat file
seurat_anndata <- LoadH5Seurat(
  "adata_SS2_for_download.h5seurat"
)

# =========================================================
# Courtesy
# =========================================================
# This script was prepared for educational purposes while
# learning single-cell RNA-seq data handling in R.
#
# Special courtesy and learning credit to:
# https://www.youtube.com/watch?v=3xcTpqQzUwQ&list=PLJefJsd1yfhagnkss5B1YCsHaH0GWQfFT&index=2
# =========================================================