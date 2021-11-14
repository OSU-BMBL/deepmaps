# Setup
source("/scratch/deepmaps/code/deepmaps.R")
Sys.setenv(RETICULATE_PYTHON = "/home/wan268/.conda/envs/hgt1/bin/python")
use_python("/home/wan268/.conda/envs/hgt1/bin/python")
py_config()

jaspar_path <- "/scratch/deepmaps/jaspar/"

################## Params
jobid <- "lymph"
result_dir <- "/scratch/deepmaps/data/7/"
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

GAS_path <- ""

################## Step 1: read file -> Generate Seurat obj & GAS
obj <-
  Read10Xdata(h5Path = data_path,
              fragmentsPath = frag_path,
              annoObj = anno, 
              min_cell = 0.001)

obj <- filterCell(obj)

#calulate peak_gene relation
ATAC_gene_peak <-
  CalGenePeakScore(peak_count_matrix = obj@assays$ATAC@counts,
                   organism = "GRCh38")

if (file.exists(GAS_path)) {
  GAS <- read.table(GAS_path, sep = " ")
} else {
  GAS <-
    calculate_GAS_v1(ATAC_gene_peak, obj, method = "velo", veloPath = velo_path)
}

#write.table(GAS, paste0(result_dir, "GAS.txt"))

#filter gene if no accessible peak in the promoter
gene_peak_pro <-
  AccPromoter(obj, ATAC_gene_peak, GAS, species = species)

################## Step 2: Run HGT

system(
  paste0(
    "/scratch/deepmaps/code/hgt.sh"
  )
)


################## Step 2: Run HGT in bash -> att & hgt_matrix
# After HGT is done:

cell_hgt_matrixPath <-
  paste0(
    "/scratch/deepmaps/data/7/cell/n_batch50_batch_size_110sample_depth_4_nheads_16_nlayers_2_sample_width_8_lr_0.2_n_hid_128_epoch_30"
  )
attPath <-
  paste0(
    "/scratch/deepmaps/data/7/att/n_batch50_batch_size_110sample_depth_4_nheads_16_nlayers_2_sample_width_8_lr_0.2_n_hid_128_epoch_30"
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

DefaultAssay(obj) <- "RNA"

#DimPlot(obj, reduction="umap.rna", group.by = "seurat_clusters")
#colnames(obj@assays$ATAC@data) <- str_replace(colnames(obj@assays$ATAC@data),"-","\\.")

#infer gene modules

co <- get_gene_module(obj, cell_hgt_matrix, att, GAS)
graph.out <- Idents(obj)
#graph.out <- droplevels(graph.out)
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
ct_regulon_v1 <- m[[2]]
TFinGAS<-m[[3]]
##calculate RI
m1<-uni(gene_peak_pro,BA_score)
peak_TF<-m1[[2]]
RI_C <-
  RI_cell(obj, ct_regulon_v1, GAS, gene_peak_pro, peak_TF, graph.out)
#RI_C <- TG_cell
##infer ct_Regulon
m <- calRAS(RI_C, ct_regulon_v1, graph.out, TFinGAS)
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
DR <-
  calDR(
    RAS_C1,
    graph.out,
    only.pos = T,
    lfcThres = 0.25,
    pvalThres = 0.05
  )

DR_all <-
  calDR(
    RAS_C1,
    graph.out,
    only.pos = F,
    lfcThres = 0,
    pvalThres = 1
  )

################## Step 6: save

case_result <- list(
  RAS=RAS,
  RAS_C=RAS_C1,
  GAS=GAS,
  RI_CT=RI_CT,
  RI_C=RI_C,
  DR=DR,
  DR_all = DR_all,
  ct_regulon=ct_regulon,
  masterTF=masterTF,
  graph.out = graph.out,
  TF_cen = TF_cen,
  gene_cen = gene_cen,
  network = network
)

#obj <- qs::qread(paste0("/scratch/deepmaps/data/lymph_obj_1102.qsave"))
qs::qsave(case_result, "/scratch/deepmaps/data/lymphoma_14k_case_result_1110.qsave")
save.image(paste0("/scratch/deepmaps/data/lymph_case_result_1110.rdata"))
qs::qsave(obj, "/scratch/deepmaps/data/lymphoma_14k_obj.qsave")
DimPlot(obj)
#load(paste0("/scratch/deepmaps/data/lymph_case_result_1102.rdata"))
