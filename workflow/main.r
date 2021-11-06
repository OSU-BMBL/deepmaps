# Setup
source("/scratch/deepmaps/code/deepmaps.R")
Sys.setenv(RETICULATE_PYTHON = "/home/wan268/.conda/envs/hgt1/bin/python")
use_python("/home/wan268/.conda/envs/hgt1/bin/python")
py_config()

jaspar_path <- "/scratch/deepmaps/jaspar/"

################## Params
base_dir <- "/home/wan268/hgt/RNA_ATAC/"
res <- 1
species <- "hg38"
obj_basename <- "lymph_14k"
anno <- qs::qread(paste0(base_dir, "hg38_annotations.qsave"))
data_path <-
  paste0(base_dir,
         "lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5")

frag_path <-
  paste0(base_dir, "lymph_node_lymphoma_14k_atac_fragments.tsv.gz")

velo_path <-
  paste0(base_dir,
         "velo/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.csv")

lisa_path <- paste0(base_dir, obj_basename, "/")


################## Step 1: read file -> Generate Seurat obj & GAS
obj <-
  Read10Xdata(h5Path = data_path,
              fragmentsPath = frag_path,
              annoObj = anno)

obj <- filterCell(obj)

GAS_path <-
  paste0(base_dir,
         "lym.txt")

GAS <- read.table(GAS_path, sep = " ")

#calulate peak_gene relation
ATAC_gene_peak <-
  CalGenePeakScore(peak_count_matrix = obj@assays$ATAC@counts,
                   organism = "GRCh38")

#filter gene if no accessible peak in the promoter
gene_peak_pro <-
  AccPromoter(obj, ATAC_gene_peak, GAS, species = species)

#GAS <-
#  calculate_GAS(ATAC_gene_peak, obj, method = "velo", veloPath = velo_path)
write.table(GAS, "GAS_test_output.txt")

################## Step 2: Run HGT in bash -> att & hgt_matrix
# After HGT is done:

cell_hgt_matrixPath <-
  paste0(
    base_dir,
    "7/cell/n_batch50_batch_size_110sample_depth_4_nheads_16_nlayers_2_sample_width_8_lr_0.2_n_hid_128_epoch_30"
  )
attPath <-
  paste0(
    base_dir,
    "7/att/n_batch50_batch_size_110sample_depth_4_nheads_16_nlayers_2_sample_width_8_lr_0.2_n_hid_128_epoch_30"
  )

# Double check if files exist
file.exists(cell_hgt_matrixPath)
file.exists(attPath)

################## Step 3: clustering & generate gene modules
#qs::qsave(obj, paste0(base_dir, obj_basename, ".qsave"))

colnames(GAS) <- str_replace_all(colnames(GAS), "\\.", "-")

cell_hgt_matrix <- read.table(cell_hgt_matrixPath)
cell_hgt_matrix <- as.matrix(cell_hgt_matrix)
rownames(cell_hgt_matrix) <- colnames(GAS)
att <- read.csv(attPath)

obj <- obj[, colnames(GAS)]
GAS <- GAS[, colnames(obj)]
cell_hgt_matrix <- cell_hgt_matrix[colnames(GAS),]

HGT_embedding <-
  CreateDimReducObject(embeddings = cell_hgt_matrix,
                       key = "HGT_",
                       assay = "RNA")
obj@reductions[['HGT']] <- HGT_embedding
obj <-
  FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <-
  RunUMAP(
    obj,
    reduction = 'HGT',
    dims = 1:ncol(cell_hgt_matrix),
    reduction.name = "umap.rna",
    reduction.key = "rnaUMAP_"
  )
obj <-
  FindNeighbors(obj,
                reduction = "HGT",
                dims = 1:ncol(cell_hgt_matrix))
obj <- FindClusters(obj, resolution = res)

DefaultAssay(obj) <- "ATAC"

#DimPlot(obj, reduction="umap.rna", group.by = "seurat_clusters")
colnames(obj@assays$ATAC@data) <- str_replace(colnames(obj@assays$ATAC@data),"-","\\.")

#infer gene modules

m <- get_gene_module(obj, cell_hgt_matrix, att, GAS)
co <- m[[1]]
graph.out <- m[[2]]

#write gene modules
dir.create(lisa_path, showWarnings = F)
write_GM(co, lisa_path)

################## Step 4: Run LISA

system(
  paste0(
    "/home/wan268/.conda/envs/lisa/bin/python /scratch/deepmaps/code/run_lisa.py --path ",
    lisa_path,
    " --species ",
    species
  )
)

################## Step 5: regulon

m <-
  Calregulon(
    GAS,
    co,
    gene_peak_pro,
    species = species,
    jaspar_path = jaspar_path,
    lisa_path = lisa_path
  )
BA_score <- m[[1]]
ct_regulon <- m[[2]]
TFinGAS<-m[[3]]
##calculate RI
m1<-uni(gene_peak_pro,BA_score)
peak_TF<-m1[[2]]
RI_C <-
  RI_cell(obj, ct_regulon, GAS, gene_peak_pro, peak_TF, graph.out)
#RI_C <- TG_cell
##infer ct_Regulon
m <- calRAS(RI_C, ct_regulon, graph.out, TFinGAS)
RAS <- m[[1]]
RI_CT <- m[[2]]
ct_regulon <- m[[3]]
RAS_C <- m[[4]]

##calculate RAS(2) same topology in a row
RAS_C1 <- CalRAS2(ct_regulon, graph.out)

##calculate master TF
masterTF <- masterFac(ct_regulon, RI_CT)
TF_cen <- masterTF[[1]]
gene_cen <- masterTF[[2]]
network <- masterTF[[3]]

##calculate DR
DR <- calDR(RAS_C1, graph.out, FindALLMarkers = T)


################## Step 6: save
save.image(paste0(base_dir, obj_basename, "_1102.rdata"))
qs::qsave(obj, "/scratch/deepmaps/data/lymph_obj_1102.qsave")

#save.image(paste0("/scratch/deepmaps/data/lymph_obj_1102.rdata"))
#load(paste0(base_dir, obj_basename, "_1102.rdata"))
#obj <- qs::qread(paste0("/scratch/deepmaps/data/lymph_obj_1102.qsave"))

case_result <- list(
  RAS=RAS,
  RAS_C=RAS_C,
  GAS=GAS,
  RI_CT=RI_CT,
  RI_C=RI_C,
  DR=DR,
  ct_regulon=ct_regulon,
  masterTF=masterTF,
  graph.out = graph.out
)
qs::qsave(case_result, "/scratch/deepmaps/data/lymph_case_result_1102.qsave")
