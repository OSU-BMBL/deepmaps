# Setup

Sys.setenv(RETICULATE_PYTHON = "/home/wan268/.conda/envs/hgt1/bin/python")
use_python("/home/wan268/.conda/envs/hgt1/bin/python")
py_config()
source("/scratch/deepmaps/code/deepmaps.R")
jaspar_path <- "/scratch/deepmaps/jaspar/"

################## Params
base_dir <- "/home/wan268/hgt/RNA_ATAC/"
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

#calulate peak_gene relation
ATAC_gene_peak <-
  CalGenePeakScore(peak_count_matrix = obj@assays$ATAC@counts,
                   organism = "GRCh38")

#filter gene if no accessible peak in the promoter
gene_peak_pro <-
  AccPromoter(obj, ATAC_gene_peak, GAS, species = "human")

GAS <-
  calculate_GAS(ATAC_gene_peak, obj, method = "velo", veloPath = velo_path)
write.table(GAS, "GAS_test_output.txt")

################## Step 2: Run HGT in bash -> att & hgt_matrix
# After HGT is done:
GAS_path <-
  paste0(base_dir,
         "lymph_node_lymphoma_14k_filtered_feature_bc_matrix.txt")

cell_hgt_matrixPath <-
  paste0(base_dir,
         "7_7/cell/930_n_hid_52_nheads_13_nlayers_2_lr_0.3")
attPath <-
  paste0(base_dir,
         "7_7/att/gas_n_hid_52_nheads_13_nlayers_2_lr_0.3")

# Double check if files exist
file.exists(cell_hgt_matrixPath)
file.exists(attPath)

################## Step 3: clustering & generate gene modules

#obj <- qs::qread(paste0(base_dir, obj_basename, ".qsave"))
GAS <- read.table(GAS_path, sep = " ")
colnames(GAS) <- str_replace_all(colnames(GAS), "\\.", "-")

cell_hgt_matrix <- read.table(cell_hgt_matrixPath)
cell_hgt_matrix <- as.matrix(cell_hgt_matrix)
att <- read.csv(attPath)
rownames(cell_hgt_matrix) <- colnames(GAS)

obj <- obj[, colnames(GAS)]
GAS <- GAS[, colnames(obj)]
cell_hgt_matrix <- cell_hgt_matrix[colnames(GAS), ]

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
obj <- FindClusters(obj, resolution = 0.5)

#infer gene modules
m <- get_gene_module(obj, cell_hgt_matrix, att, GAS)
graph.out <- m[[1]]
co <- m[[2]]

#write gene modules
dir.create(lisa_path, showWarnings = F)
write_GM(co, lisa_path)

################## Step 4: Run LISA

system(
  paste0(
    "/home/wan268/.conda/envs/lisa/bin/python /scratch/deepmaps/code/run_lisa.py --path ",
    lisa_path
  )
)

################## Step 5: regulon

m <-
  Calregulon(
    GAS,
    co,
    gene_peak_pro,
    speices = "hg38",
    jaspar_path = jaspar_path,
    lisa_path = lisa_path
  )
BA_score <- m[[1]]
ct_regulon <- m[[2]]

##calculate RI
RI_C <-
  RI_cell(obj, ct_regulon, GAS, gene_peak_pro, BA_score, graph.out)

##infer ct_Regulon
m <- calRAS(RI_C, ct_regulon, graph.out)
RAS <- m[[1]]
RI_CT <- m[[2]]
ct_re <- m[[3]]
RAS_C <- m[[4]]
ct_regulon <- ct_re

##calculate RAS(2) same topology in a row
RAS_C1 <- CalRAS2(ct_regulon, graph.out)

##calculate master TF
m <- masterFac(ct_regulon, RI_CT)
TF_cen <- m[[1]]
gene_cen <- m[[2]]
network <- m[[3]]

##calculate DR
DR <- calDR(RAS_C1, graph.out, FindALLMarkers = T)

##calculate VR

library(iterators)
library(foreach)
library(doParallel)

cl <- makeCluster(12)
registerDoParallel(cl)

graph.out0 <-
  (graph.out)[which(graph.out %in% as.integer(substring(colnames(RI_CT), 3, nchar(colnames(
    RI_CT
  )))))]
#names(graph.out0)<-names(graph.out)
tfsmatrix <- cal_tfmatrixs(graph.out0, graph.out, RI_C, ct_regulon)
print(str(tfsmatrix))
#cl<-makeCluster(detectCores()-2)
#registerDoParallel(cl)
clusterExport(cl, "tfsmatrix")
clusterExport(cl, "graph.out0")
clusterExport(cl, c("cal_sp", "Score_matrix", "cal_meanmatrix", "cal_clust"))

VR <-
  foreach(i = names(tfsmatrix), .combine = 'rbind') %dopar% cal_sp(i)
rownames(VR) <- names(tfsmatrix)

################## Step 6: save
save.image(paste0(base_dir, obj_basename, ".rdata"))
#load(paste0(base_dir, obj_basename, ".rdata"))
qs::qsave(obj, "/scratch/deepmaps/data/lymph_obj.qsave")
