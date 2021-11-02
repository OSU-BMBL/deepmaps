# Setup
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/wan268/.conda/envs/hgt1/bin/python")
use_python("/home/wan268/.conda/envs/hgt1/bin/python")
py_config()
source("/scratch/deepmaps/code/deepmaps.R")
jaspar_path <- "/scratch/deepmaps/jaspar/"

################## Params 1
base_dir <- "/home/wan268/hgt/RNA_ATAC/"

obj_basename <- "lymph_14k"
lisa_path <- paste0(base_dir, obj_basename, "/")
GAS_path <-
  paste0(base_dir,
         "lymph_node_lymphoma_14k_filtered_feature_bc_matrix.txt")
cell_hgt_matrixPath <-
  paste0(base_dir,
         "7_7/cell/930_n_hid_52_nheads_13_nlayers_2_lr_0.3")

attPath <-
  paste0(base_dir,
         "7_7/att/gas_n_hid_52_nheads_13_nlayers_2_lr_0.3")

file.exists(cell_hgt_matrixPath)
file.exists(attPath)
##################

################## Params 2
base_dir <- "/home/wan268/hgt/RNA_ATAC/"

obj_basename <- "pbmc_10k_unsort"
lisa_path <- paste0(base_dir, obj_basename, "/")
GAS_path <-
  paste0(base_dir,
         "pbmc_unsorted_10k_filtered_feature_bc_matrix.txt")
cell_hgt_matrixPath <-
  paste0(base_dir,
         "3_3/cell/930_n_hid_96_nheads_16_nlayers_2_lr_0.4")

attPath <-
  paste0(base_dir,
         "3_3/att/gas_n_hid_96_nheads_16_nlayers_2_lr_0.4")

##################

obj <- qs::qread(paste0(base_dir, obj_basename, ".qsave"))
GAS <- read.table(GAS_path, sep = " ")
colnames(GAS) <- str_replace_all(colnames(GAS), "\\.", "-")

cell_hgt_matrix <- read.table(cell_hgt_matrixPath)
cell_hgt_matrix <- as.matrix(cell_hgt_matrix)
att <- read.csv(attPath)
rownames(cell_hgt_matrix) <- colnames(GAS)

obj <- obj[, colnames(GAS)]
GAS <- GAS[, colnames(obj)]
cell_hgt_matrix <- cell_hgt_matrix[colnames(GAS),]
#att <-

#
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
m <- gene_module(obj, cell_hgt_matrix, att, GAS)
graph.out <- m[[1]]
co <- m[[2]]

#write gene modules
dir.create(lisa_path, showWarnings = F)
write_GM(co, lisa_path)

###
#
###

system(
  paste0(
    "/home/wan268/.conda/envs/lisa/bin/python /scratch/deepmaps/code/run_lisa.py --path ",
    lisa_path
  )
)

ATAC_gene_peak <-
  CalGenePeakScore(peak_count_matrix = obj@assays$ATAC@counts)

#filter gene if no accessible peak in the promoter
gene_peak_pro <-
  AccPromoter(obj, ATAC_gene_peak, GAS, species = "human")

#infer candidate regulon
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
#m <- uni(gene_peak_pro, BA_score)
#mat <- m[[1]]
#mat1 <- m[[2]]
#peak_TF <- mat1
#gene_peak <- gene_peak_pro

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

#regulon information
#rm(peak)
save.image(paste0(base_dir, obj_basename, ".rdata"))
#load(paste0(base_dir, obj_basename, ".rdata"))
qs::qsave(obj, "/scratch/deepmaps/data/lymph_obj.qsave")



library(RColorBrewer)
display.brewer.all(type = "div")
barplot(lengths(co))
barplot(table(table(unlist(co))))

a <-
  unlist(strsplit(names(ct_regulon), "_"))[seq(2, length(unlist(strsplit(
    names(ct_regulon), "_"
  ))), 2)]
#the times of gene number in each ct
barplot(table(a), col = brewer.pal(11, "BrBG")[7])
#the number of gene in each regulon
barplot(lengths(ct_regulon), col = brewer.pal(11, "BrBG")[7])
b <-
  unlist(strsplit(names(ct_regulon), "_"))[seq(1, length(unlist(strsplit(
    names(ct_regulon), "_"
  ))), 2)]
#the time of TF appear in different CT
barplot(table(b), col = brewer.pal(11, "BrBG")[7])
barplot(table(table(b)), col = brewer.pal(11, "BrBG")[7])
#enrichment analysis
library(enrichR)
dbs <- c("KEGG_2019_Human")
res <- list()
for (i in 1:length(ct_regulon)) {
  genes <- unlist(ct_regulon[[i]])
  enrichr_res <- enrichr(genes, dbs)
  res[[names(ct_regulon[i])]] <-
    enrichr_res[[1]][enrichr_res$KEGG_2019_Human$Adjusted.P.value < 0.05, ]
  print(enrichr_res[[1]][enrichr_res$KEGG_2019_Human$Adjusted.P.value <
                           0.05, ])
  #this_pval <- enrichr_res$KEGG_2019_Human$Adjusted.P.value[1]
  #old_pvalue <- append(old_pvalue, this_pval)
}
zz <- list()
for (i in (1:length(res))) {
  zz <- append(zz, nrow(res[[i]]))
  
}
#the number of enriched regulons
length(unlist(zz)[unlist(zz) > 0])
#the number of enriched pathway
sum(unlist(zz)[unlist(zz) > 0])
#the name of enriched  regulon with CT
enr <- names(res[unlist(zz) > 0])
#the name of enriched regulon without CT
enr2 <-
  unlist(strsplit(enr, "_"))[seq(1, length(unlist(strsplit(enr, "_"))), 2)]
ct <-
  unlist(strsplit(enr, "_"))[seq(2, length(unlist(strsplit(enr, "_"))), 2)]
#unique enriched pathway
term <- list()
for (i in (1:length(res[enr]))) {
  #print(res[enr][[i]]$Term)
  term <- c(term, c(res[enr][[i]]$Term))
}
Uterm <- unique(unlist(term))
DR <- DR[DR$p_val_adj < 0.05, ]
for (i in unique(ct)) {
  e1 <- enr[i == ct]
  
}