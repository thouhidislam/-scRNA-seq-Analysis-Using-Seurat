
# Single-Cell RNA-seq Analysis of NSCLC Using Seurat

## Project Overview
This project demonstrates a complete workflow for analyzing single-cell RNA sequencing (scRNA-seq) data from a Non-Small Cell Lung Cancer (NSCLC) dataset.  
The workflow uses the Seurat R package to perform quality control, normalization, dimensionality reduction, clustering, and visualization of single-cell data.

## Key Steps
1. **Data Loading** – Import 10X Genomics H5 matrix into R.  
2. **Quality Control** – Filter low-quality cells based on gene count and mitochondrial content.  
3. **Normalization** – Log-normalize counts for comparability across cells.  
4. **Variable Feature Selection** – Identify genes with high variability for downstream analysis.  
5. **Scaling & PCA** – Scale the data and perform Principal Component Analysis.  
6. **Clustering** – Build a KNN graph and cluster cells at multiple resolutions.  
7. **Visualization** – Generate UMAP plots with cluster annotations.  
8. **Save Results** – Save the processed Seurat object for future analysis.

## Data
- Input: `20k_NSCLC_DTC_3p_nextgem_donor_1_count_sample_feature_bc_matrix.h5`  
- Source: 10X Genomics single-cell dataset (NSCLC).  

> **Note:** Raw data files are not included due to size constraints. Please place the H5 file in the `data/` directory before running the script.

## Tools & Packages
- **R** (v4.x recommended)  
- **Seurat** – For single-cell analysis  
- **tidyverse** – Data manipulation and visualization  
- **hdf5r** – Reading 10X H5 format files  

## Outputs
- Processed Seurat object: `results/seurat_nsclc_final.rds`  
- Cluster plots and QC plots: save in `figures/` folder  

## Author
Your Name – Bioinformatics | Single-Cell RNA Sequencing  
Email: your.email@example.com  

## License
This project is licensed under the MIT License. See `LICENSE` file for details
