#!/usr/bin/env Rscript

# Setup ========================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
  library(ggplot2)
})

# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e9)

# Parse command line arguments
args <- commandArgs(trailingOnly=TRUE)

bin_size <- args[1]
n_sketch_cells <- args[2]
cluster_res <- args[3]
cluster_npcs <- args[4]
cluster_ndims <- 1:cluster_npcs
sample_name <- args[5]

# check for valid bin size
if (!any(bin_size == c(2,4,8,16))) {
  stop("Invalid bin size provided. Should be one of 2, 4, 8, or 16")
}

# get assay name for given bin size
if (bin_size == 16) {
  assay_name <- "Spatial.016um"
} else {
  assay_name <- paste0("Spatial.00",bin_size,"um")
}

# import data and convert to BPcells format on-disk ============================
file_path <- "./outs/"
sdata <- Load10X_Spatial(file_path, bin.size = bin_size)

count_mat_int <- convert_matrix_type(sdata[[assay_name]]$counts, type = "double")

mat_dir <- paste0("BPcells_", assay_name)


write_matrix_dir(mat = count_mat_int, dir = mat_dir, overwrite = TRUE)
counts.mat <- open_matrix_dir(dir = mat_dir)

sdata[[assay_name]]$counts <- counts.mat

rm(counts.mat, count_mat_int)

# Create sketch assay with subset of cells =====================================
sdata <- NormalizeData(sdata)
sdata <- FindVariableFeatures(sdata)
sdata <- SketchData(
  object = sdata,
  ncells = n_sketch_cells,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# Perform clustering on sketch subset ==========================================
DefaultAssay(sdata) <- "sketch"
sdata <- FindVariableFeatures(sdata)
sdata <- ScaleData(sdata)
sdata <- RunPCA(sdata, npcs = cluster_npcs)
sdata <- FindNeighbors(sdata, dims = cluster_ndims)
sdata <- FindClusters(sdata, resolution = cluster_res)

# compute UMAP embedding =======================================================
sdata <- RunUMAP(sdata, dims = cluster_ndims, return.model = T)
sketch_clusters_plot <- DimPlot(sdata, label = T, label.size = 3, reduction = "umap") + NoLegend()
sketch_spatial_plot <- SpatialDimPlot(sdata, group.by = "seurat_clusters") + NoLegend()
# project clustering onto full dataset =========================================
sdata <- ProjectData(
  object = sdata,
  assay = assay_name,
  full.reduction = "pca.full",
  sketched.assay = "sketch",
  sketched.reduction = "pca",
  umap.model = "umap",
  dims = cluster_ndims,
  refdata = list(cluster_full = "seurat_clusters")
)
# now that we have projected the full dataset, switch back to analyzing all cells
DefaultAssay(sdata) <- assay_name
full_clusters_plot <- DimPlot(sdata, label = T, label.size = 3, reduction = "full.umap", group.by = "cluster_full", alpha = 0.1) + NoLegend()
full_spatial_plot <- SpatialDimPlot(sdata, group.by = "cluster_full") + NoLegend()

# generate plots with umap embedding ===========================================
plot_filename <- paste0(sample_name,"_",bin_size, "um_clusters_UMAP.pdf")
pdf(plot_filename, useDingbats = FALSE, width = 11)
(sketch_clusters_plot | full_clusters_plot) /
  (sketch_spatial_plot | full_spatial_plot)
dev.off()

# write clusters to file =======================================================
library(tidyverse)
clusters <- sdata@meta.data |> 
  as.data.frame() |>
  rownames_to_column("spot_id") |> 
  as_tibble() |> 
  mutate(spot_id = str_replace(spot_id, "s_", "")) |> 
  select(spot_id, cluster_full)

out_fn <-  paste0(sample_name,"_",bin_size,"um_clusters.csv.gz")
write_csv(clusters, out_fn)

