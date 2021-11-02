# Setup
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "/home/wan268/.conda/envs/hgt1/bin/python")
use_python("/home/wan268/.conda/envs/hgt1/bin/python")
py_config()
source("/scratch/deepmaps/code/deepmaps.R")



################## Params 1
base_dir <- "/home/wan268/hgt/RNA_ATAC/"
data_path <-
  paste0(base_dir,
         "lymph_node_lymphoma_14k_filtered_feature_bc_matrix.h5")
frag_path <-
  paste0(base_dir, "lymph_node_lymphoma_14k_atac_fragments.tsv.gz")
anno <- qs::qread(paste0(base_dir, "hg38_annotations.qsave"))

output_basename <- "lymph_14k"
velo <- ""
##################

################## Params 2
base_dir <- "/home/wan268/hgt/RNA_ATAC/"
data_path <-
  paste0(base_dir, "pbmc_unsorted_10k_filtered_feature_bc_matrix.h5")
frag_path <-
  paste0(base_dir, "pbmc_unsorted_10k_atac_fragments.tsv.gz")
anno <- qs::qread(paste0(base_dir, "hg38_annotations.qsave"))

output_basename <- "pbmc_10k_unsort"
velo <- ""

GAS_path <-
  paste0(base_dir,
         "pbmc_unsorted_10k_filtered_feature_bc_matrix.txt")
##################

raw_obj <-
  Read10Xdata(h5Path = data_path,
              fragmentsPath = frag_path,
              annoObj = anno)

obj <-
  readmatrix(raw_obj@assays$RNA@counts,
             raw_obj@assays$ATAC@counts,
             min_cell = 0.1)

#data preprocess
obj <- filterCell(obj)

#calulate peak_gene relation
ATAC_gene_peak <-
  CalGenePeakScore(peak_count_matrix = obj@assays$ATAC@counts,
                   organism = "GRCh38")

qs::qsave(obj, paste0(base_dir, output_basename, ".qsave"))



#integration
#GAS <- calculate_GAS(ATAC_gene_peak, obj, method = "velo", veloPath = velo)
#write.table(GAS, GAS_path)


#STEP 2

GAS <- read.table(GAS_path, sep = " ")

#infer gene modules
cell_hgt_matrixPath <-
  paste0(
    base_dir,
    "3/cell/n_batch50_batch_size_110sample_depth_5_nheads_13_nlayers_2_sample_width_11_lr_0.3_n_hid_78_redution_AE_rf_0.0_factor_0.5_pacience_5_layertype_hgt_loss_kl_optimizer_adamw_dropout_0.0"
  )

attPath <-
  paste0(
    base_dir,
    "3/att/n_batch50_batch_size_110sample_depth_5_nheads_13_nlayers_2_sample_width_11_lr_0.3_n_hid_78_redution_AE_rf_0.0_factor_0.5_pacience_5_layertype_hgt_loss_kl_optimizer_adamw_dropout_0"
  )

m <- gene_module(obj, cell_hgt_matrixPath, attPath, GAS)
graph.out <- m[[1]]
co <- m[[2]]

#write gene modules
write_GM(co)

###
#python run lisa.py
###

#filter gene if no accessible peak in the promoter
gene_peak_pro <-
  AccPromoter(obj, ATAC_gene_peak, GAS, species = "human")

#infer candidate regulon
m <- Calregulon(GAS, co, gene_peak_pro, speices = "hg38")
BA_score <- m[[1]]
ct_regulon <- m[[2]]
m <- uni(gene_peak_pro, BA_score)
mat <- m[[1]]
mat1 <- m[[2]]
peak_TF <- mat1
gene_peak <- gene_peak_pro

##calculate RI
RI_C <-
  RI_cell(obj, ct_regulon, GAS, gene_peak_pro, peak_TF, graph.out)

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

cl <- makeCluster(19)
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
    enrichr_res[[1]][enrichr_res$KEGG_2019_Human$Adjusted.P.value < 0.05,]
  print(enrichr_res[[1]][enrichr_res$KEGG_2019_Human$Adjusted.P.value <
                           0.05,])
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
DR <- DR[DR$p_val_adj < 0.05,]
for (i in unique(ct)) {
  e1 <- enr[i == ct]
  
}