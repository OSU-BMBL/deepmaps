library(dplyr)
library(Seurat)
library(patchwork)
library(cluster)
myArgs <- commandArgs(trailingOnly = TRUE)
GAS <- read.table(myArgs[1])
l <- read.csv(myArgs[2])
#GAS<-read.table("/fs/ess/PCON0022/Xiaoying/RNA/Klein/process.txt")
#l<-read.csv('/fs/ess/PCON0022/wxy/RNA_test_data/Raw_data/12.Klein/Klein_cell_lable.csv')
label <- l$Cluster
names(label) <- l$Cell
#label<-label[colnames(GAS)]
kidney0 <- CreateSeuratObject(GAS)
cell_dir <- myArgs[3]
#print(dim(kidney0))
#print(length(label))
#cell_dir<-'/fs/ess/PCON0022/wxy/raw/cell'
#gene_dir<-'/fs/ess/PCON0022/wxy/raw/gene'
#result_dir<-'/fs/ess/PCON0022/wxy/raw/print1.txt'
silandplot <- function(cell_hgt_matrix) {
  rownames(cell_hgt_matrix) <- colnames(GAS)
  label <- label[colnames(GAS)]
  #rownames(gene_hgt_matrix) <- rownames(GAS)
  HGT_embedding <-
    CreateDimReducObject(embeddings = cell_hgt_matrix,
                         #loadings = gene_hgt_matrix ,
                         key = "HGT_",
                         assay = "RNA")
  kidney <- kidney[, colnames(GAS)]
  kidney@reductions[['HGT']] <- HGT_embedding
  kidney <-
    FindVariableFeatures(kidney, selection.method = "vst", nfeatures = 2000)
  kidney <- ScaleData(kidney, features = VariableFeatures(kidney))
  #kidney <- RunPCA(kidney)
  kidney <-
    RunUMAP(
      kidney,
      reduction = 'HGT',
      dims = 1:nhid,
      reduction.name = "umap.rna",
      reduction.key = "rnaUMAP_"
    )
  kidney <- FindNeighbors(kidney, reduction = "HGT", dims = 1:nhid)
  #rownames(kidney@graphs$graph.name)
  kidney <- FindClusters(kidney, resolution = 0.5)
  #DimPlot(kidney, reduction = 'umap.rna')
  #Idents(kidney)
  
  ARI <-
    igraph::compare(Idents(kidney), as.factor(label), method = "adjusted.rand")
  sil1 <-
    silhouette(as.numeric(Idents(kidney)), dist(cell_hgt_matrix))
  sil2 <-
    silhouette(as.numeric(Idents(kidney)),
               dist(kidney@reductions$umap.rna@cell.embeddings))
  graph.out <- Idents(kidney)
  
  #ARI<-igraph::compare(as.numeric(graph.out),as.factor(label),method="adjusted.rand")
  sil <- list()
  sil[['sil1']] <- sil1
  sil[['sil2']] <- sil2
  sil[['ctn']] <- length(unique(Idents(kidney)))
  sil[['ARI']] <- ARI
  return(sil)
}
file = myArgs[4]
#file = "n_batch64_batch_size_64sample_depth_4_nheads_8_nlayers_4_sample_width_8_lr_0.5_n_hid_64_redution_AE_layertype_hgt_loss_kl"
nhid = unlist(strsplit(file, "n_hid_"))[2]
nhid = unlist(strsplit(nhid, "_redution_"))[1]
kidney <- kidney0
cell_hgt_matrix <- read.table(file.path(cell_dir, file))
cell_hgt_matrix <- as.matrix(cell_hgt_matrix)
sil <- silandplot(cell_hgt_matrix)
sil1 <- summary(sil[['sil1']])$avg.width
sil2 <- summary(sil[['sil2']])$avg.width
ct_number <- sil[['ctn']]
ARI <- sil[['ARI']]
cat(ARI, sil1)
