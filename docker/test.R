# test package library and install
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE))
  BiocManager::install("EnsDb.Hsapiens.v86")
if (!requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE))
  BiocManager::install("EnsDb.Mmusculus.v79")
if (!requireNamespace("scater", quietly = TRUE))
  BiocManager::install("scater")
if (!requireNamespace("bluster", quietly = TRUE))
  BiocManager::install("bluster")
if (!requireNamespace("GenomeInfoDb", quietly = TRUE))
  BiocManager::install("GenomeInfoDb")
if (!requireNamespace("GenomeInfoDb", quietly = TRUE))
  BiocManager::install("GenomeInfoDb")
if (!requireNamespace("IRanges", quietly = TRUE))
  BiocManager::install("IRanges")
if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")
if (!requireNamespace("rtracklayer", quietly = TRUE))
  BiocManager::install ("rtracklayer")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
if (!requireNamespace("RColorBrewer", quietly = TRUE))
  install.packages('RColorBrewer')
if (!requireNamespace("reticulate", quietly = TRUE))
  install.packages("reticulate")
if (!requireNamespace("plyr", quietly = TRUE))
  install.packages('plyr')
if (!requireNamespace("dsb", quietly = TRUE))
  install.packages('dsb')
if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages('Seurat')
if (!requireNamespace("Signac", quietly = TRUE))
  install.packages("Signac")
if (!requireNamespace("cluster", quietly = TRUE))
  install.packages("cluster")
if (!requireNamespace("igraph", quietly = TRUE))
  install.packages ("igraph")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages ("ggplot2")
if (!requireNamespace("Matrix", quietly = TRUE))
  install.packages ("Matrix") 
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!requireNamespace("tinytex", quietly = TRUE))
  install.packages("tinytex") 
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages ( "tidyverse" ) 
if (!requireNamespace("devtools", quietly = TRUE))
  library(devtools)
if (!requireNamespace("MAESTRO", quietly = TRUE))
  install_github("liulab-dfci/MAESTRO")
if (!requireNamespace("data.table", quietly = TRUE))
  install.packages("data.table")
if (!requireNamespace("igraph", quietly = TRUE))
  install.packages ("igraph")
if (!requireNamespace("parallel", quietly = TRUE))
  install.packages("parallel")
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr") 
if (!requireNamespace("Hmisc", quietly = TRUE))
  install.packages("Hmisc")
if (!requireNamespace("CellChat", quietly = TRUE))
  devtools::install_github("sqjin/CellChat")
if (!requireNamespace("patchwork", quietly = TRUE))
  devtools::install_github("thomasp85/patchwork")
# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
if (!requireNamespace("Seurat", quietly = TRUE))
  install.packages('Seurat')
if (!requireNamespace("plyr", quietly = TRUE))
  install.packages('plyr')
if (!requireNamespace("dsb", quietly = TRUE))
  install.packages('dsb')
if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")
if (!requireNamespace("RColorBrewer", quietly = TRUE))
  install.packages('RColorBrewer')
if (!requireNamespace("reticulate", quietly = TRUE))
  install.packages("reticulate")
if (!requireNamespace("CellChat", quietly = TRUE))
  devtools::install_github("sqjin/CellChat")
if (!requireNamespace("patchwork", quietly = TRUE))
  devtools::install_github("thomasp85/patchwork")
if (!requireNamespace("scater", quietly = TRUE))
  BiocManager::install("scater")
  
library(MAESTRO)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Mmusculus.v79)
library(scater)
library(Seurat)
library(Signac)
library(cluster)
library(bluster)
library(GenomeInfoDb)
library(igraph)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(Matrix)
library(dplyr)
library(tinytex)
library(tidyverse)
library(rtracklayer)
library(reticulate)
library(dplyr)
library(parallel)
library(igraph)
library(data.table)
library(Hmisc)
library(dplyr)
library(Seurat)
library(patchwork)
library(cluster)

## Setup path and environment
source("/deepmaps/scRNA_scATAC1.r")
# set python environment
Sys.setenv(RETICULATE_PYTHON = "/home/user/miniconda/bin/python")
use_python("/home/user/miniconda/bin/python")
py_config()
lisa_path <- "/deepmaps/lisa_output/"
jaspar_path <- "/home/user/miniconda/lib/python3.8/site-packages/lisa/data/"


# Read matched scRNA and scATAC data and quality control
lymph_obj <- ReadData(h5Path = "/deepmaps/docker/data/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5", data_type = "scRNA_scATAC", min_cell = 0.001, dataFormat = "h5")
ATAC_gene_peak <- CalGenePeakScore(peak_count_matrix = lymph_obj@assays$ATAC@counts,organism = "GRCh38")

# Calculate gene active score (integeation)
velo <-"/deepmaps/docker/data/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv"
GAS_obj <- calculate_GAS_v1(ATAC_gene_peak = ATAC_gene_peak, obj = lymph_obj, method = "velo", veloPath = velo)
GAS <- GAS_obj[[1]]
lymph_obj <- GAS_obj[[2]]

# Perform heterogeneous graph transformer (HGT) model on GAS matrix
HGT_result <- run_HGT(GAS = as.matrix(GAS),result_dir='/deepmaps', data_type='scRNA_scATAC', envPath=NULL, lr=0.2, epoch=30, n_hid=128, n_heads=16)

# Cluster the cells
cell_hgt_matrix <- HGT_result[['cell_hgt_matrix']]
rownames(cell_hgt_matrix) <- colnames(GAS)

lymph_obj <- lymph_obj[, colnames(GAS)]
cell_hgt_matrix <- cell_hgt_matrix[colnames(GAS),]

HGT_embedding <-
  CreateDimReducObject(embeddings = cell_hgt_matrix,
                       key = "HGT_",
                       assay = "RNA")
lymph_obj@reductions[['HGT']] <- HGT_embedding
lymph_obj <-
  FindVariableFeatures(lymph_obj, selection.method = "vst", nfeatures = 2000)
lymph_obj <- ScaleData(lymph_obj, features = VariableFeatures(lymph_obj))
lymph_obj <-
  RunUMAP(
    lymph_obj,
    reduction = 'HGT',
    dims = 1:ncol(cell_hgt_matrix),
    reduction.name = "umap.rna",
    reduction.key = "rnaUMAP_"
  )
lymph_obj <-
  FindNeighbors(lymph_obj,
                reduction = "HGT",
                dims = 1:ncol(cell_hgt_matrix))
lymph_obj <- FindClusters(lymph_obj, resolution = 1)
graph.out <- as.factor(lymph_obj$seurat_clusters)

DefaultAssay(lymph_obj) <- "RNA"
png("plot.png")
DimPlot(lymph_obj, reduction = 'umap.rna')
dev.off()

# Calculate cell cluster active gene modules and run LISA for TF infer.
co <- get_gene_module(obj = lymph_obj, GAS = GAS, att = HGT_result[['attention']],method = 'SFP' )

dir.create(lisa_path, showWarnings = F)
write_GM(co = co, lisa_path = lisa_path)

# Run LISA on generated gene modules
system(
  paste0(
    "/home/user/miniconda/bin/python /deepmaps/run_lisa.py --path ",
    lisa_path,
    " --species ",
    "hg38"
  )
)

# Filter gene with no accessible peak in promoter
gene_peak_pro <- AccPromoter(obj = lymph_obj, gene_peak = ATAC_gene_peak, GAS = GAS, species = 'hg38')

pre_regulon_res <- Calregulon(GAS = GAS, co = co,gene_peak_pro = gene_peak_pro, species = "hg38", jaspar_path = jaspar_path, lisa_path = lisa_path)
BA_score <- pre_regulon_res[[1]]
ct_regulon_v1 <- pre_regulon_res[[2]]
TFinGAS<- pre_regulon_res[[3]]

# Combine same TFs
peak_TF <- uni(gene_peak_pro = gene_peak_pro, BA_score = BA_score)
head(peak_TF[1:3,1:10])

# Calculate regulatory Intensive (RI) score in cell level and infer cell type active regulon
RI_C <- RI_cell(obj = lymph_obj, ct_regulon = ct_regulon_v1, GAS = GAS, gene_peak_pro = gene_peak_pro, peak_TF = peak_TF, graph.out = graph.out)
head(RI_C[1:3,1:10])

# Calculate regulon active score1 (RAS1)
regulon_res <- calRAS(RI_C = RI_C, ct_regulon = ct_regulon_v1, graph.out = graph.out, TFinGAS = TFinGAS)
RAS_CT <- regulon_res[[1]]
RI_CT <- regulon_res[[2]]
ct_regulon <- regulon_res[[3]]
RAS_C1 <- regulon_res[[4]]
head(ct_regulon[1:3])

# Calculate RAS2 same topology in a row
RAS_C2 <- CalRAS2(ct_regulon, graph.out)
head(RAS_C2[1:3,1:10])

# Calculate master TF
masterTF <- masterFac(ct_regulon = ct_regulon, RI_CT = RI_CT)
TF_cen <- masterTF[[1]]
gene_cen <- masterTF[[2]]
network <- masterTF[[3]]

# Calculate cel type specific regulons
DR <-
  calDR(
    ct_regulon = RAS_C2,
    graph.out = graph.out,
    only.pos = T,
    lfcThres = 0.25,
    pvalThres = 0.05
  )
DR[1:2]