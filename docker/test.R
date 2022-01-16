## Load function
source("/deepmaps/scRNA_scATAC1.r")
# set python environment
Sys.setenv(RETICULATE_PYTHON = "/home/user/miniconda/bin/python")
use_python("/home/user/miniconda/bin/python")
py_config()

lymph_obj <- ReadData(h5Path = "/deepmaps/docker/data/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5", data_type = "scRNA_scATAC", min_cell = 0.001, dataFormat = "h5")
ATAC_gene_peak <- CalGenePeakScore(peak_count_matrix = lymph_obj@assays$ATAC@counts,organism = "GRCh38")

velo <-"/deepmaps/docker/data/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv"
GAS_obj <- calculate_GAS_v1(ATAC_gene_peak = ATAC_gene_peak, obj = lymph_obj, method = "velo", veloPath = velo)
GAS <- GAS_obj[[1]]
lymph_obj <- GAS_obj[[2]]

HGT_result <- run_HGT(GAS = as.matrix(GAS),result_dir='/deepmaps', data_type='scRNA_scATAC', envPath=NULL, lr=0.2, epoch=30, n_hid=128, n_heads=16)

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

DefaultAssay(lymph_obj) <- "RNA"
png("plot.png")
DimPlot(lymph_obj, reduction = 'umap.rna')
dev.off()