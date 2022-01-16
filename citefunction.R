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


# Creat a seuart object
# input:
# h5Path: address of  h5 file to read. (When dataFormat is 'h5', it couldn't be NULL.)
# rna_matrix: a RNA matrix where the rows correspond to genes and the columns correspond to cells. (When dataFormat is 'matrixs', it couldn't be NULL.)
# atac_matrix: a matrix where the rows correspond to peaks and the columns correspond to cells. (When data_type is 'scRNA_scATAC'and dataFormat is 'matrixs', it couldn't be NULL.)
# adt_matrix: a matrix where the rows correspond to proteins and the columns correspond to cells. (When data_type is 'CITE' and dataFormat is 'matrixs', it couldn't be NULL.)
# data_type: the type of your input data('CITE' or 'scRNA_scATAC')
# dataFormat: the format of your input data ('matrixs' or 'h5')
# min_cell: the peak / gene will be removed if the value in the gene / peak with more than min_cell cell is equal to zero
# nmad: a numeric scalar, specifying the minimum number of MADs away from median required for a value to be called an outlier
# gene.filter: if do gene filtering
# gene.filter: if do cell filtering
# output:
# - obj: a seuart object for cite-seq data after normalizing and scaling (When data_type is 'CITE'); a seuart object for scRNA-seq and scATAC-seq (When data_type is 'scRNA_scATAC')

ReadData <- function(h5Path = NULL, rna_matrix = NULL, atac_matrix = NULL, adt_matrix = NULL, data_type = NULL, dataFormat = NULL, min_cell=0.1, nmad=3, gene.filter=TRUE, cell.filter=TRUE){
  if (data_type=='CITE'){
    # obtain RNA and ADT matrixs
    if (dataFormat=='matrixs'){
      rna <- rna_matrix
      adt <- adt_matrix
    }else if(dataFormat=='h5'){
      h5 <- Read10X_h5(h5Path)
      rna <- h5$`Gene Expression`
      adt <- h5$`Antibody Capture`

    }
    # gene filtering
    if (gene.filter==TRUE){
      binaryrna <- rna
      binaryadt <- adt
      binaryrna[binaryrna>0] <-1
      binaryadt[binaryadt>0] <-1
      rna <- rna[which(rowSums(binaryrna) > ncol(binaryrna)*min_cell),]
      adt <- adt[which(rowSums(binaryadt) > ncol(binaryadt)*min_cell),]
    }

    #Setup a Seurat object, add the RNA and protein data
    obj <- CreateSeuratObject(counts = rna)
    obj [["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
    adt_assay <- CreateAssayObject(counts = adt)
    obj[["ADT"]] <- adt_assay
    # cell filtering
    if (cell.filter==TRUE){
      adtn<-isOutlier(
        obj$nCount_ADT,
        nmads = nmad,
        log = F,
        type = "both"
      )
      rnan<-isOutlier(
        obj$nCount_RNA,
        nmads = nmad,
        log = F,
        type = "both"
      )
      mito<-isOutlier(
        obj$percent.mt,
        nmads = nmad,
        log = F,
        type = "both"
      )
      obj <-
        AddMetaData(obj, adtn, col.name = "adtn")
      obj <-
        AddMetaData(obj, rnan, col.name = "rnan")
      obj <-
        AddMetaData(obj, mito, col.name = "mito")
      obj<-subset(
        x = obj,
        subset = adtn == F &
          rnan == F &
          mito == F
      )
    }
    # normalization and scaling
    DefaultAssay(obj) <- 'RNA'
    obj <- NormalizeData(obj) %>% FindVariableFeatures()
    all.genes <- rownames(obj)
    obj <- ScaleData(obj, features = all.genes)
    DefaultAssay(obj) <- 'ADT'
    VariableFeatures(obj) <- rownames(obj[["ADT"]])
    obj <- NormalizeData(obj, normalization.method = 'CLR', margin = 2) %>%
      ScaleData()
    DefaultAssay(obj) <- 'RNA'
  }else if (data_type=='scRNA_scATAC'){
    if (gene.filter==FALSE){
      min_cell = 0
    }
    if (dataFormat == "h5") {
      obj <- Read10Xdata(h5Path = h5Path, min_cell = min_cell)
    }else{
      obj <- readmatrix(rna_matrix = rna_matrix, atac_matrix = atac_matrix, min_cell = min_cell)
    }
    if (cell.filter==TRUE){
      obj <- filterCell(obj, nmad= nmad)
    }
  }
  return(obj)
}


# splice togethor the RNA matrix and ADT matrix and do normalization
# input:
# obj: a seuart object obtained from creatobject fuction
# output:
# GAS: the spliced and normalized matrix as the input of HGT model
CLR <- function(obj){
  m1 <- obj@assays$RNA@counts[obj@assays$RNA@var.features,]
  m2 <- obj@assays$ADT@counts
  m3 <- m2
  # Add 'p_' to the name of the protein
  rownames(m3) <- paste0('p_',rownames(m2))
  m <- NormalizeData(rbind( m1, m3), normalization.method = 'CLR', margin = 2)
  return(m)
}


# Cluster cells based on HGT embedding
# input:
# cell_hgt_matrixPath: the path of cell_embedding matrix
# nhid: hyper-parameter of HGT model
# resolution: resolution of cell clustering
# obj: a seurat object obtained from creatobject fuction
# output:
# obj: a seurat object containing cell clustering result
HGT_cluster <- function(obj, cell_hgt_matrix, nhid, resolution){
  cell_hgt_matrix<-as.matrix(cell_hgt_matrix)
  rownames(cell_hgt_matrix) <-colnames(rna)
  HGT_embedding <-
    CreateDimReducObject(
      embeddings = cell_hgt_matrix,
      key = "HGT_",
      assay = "RNA"
    )
  obj<-obj[,colnames(rna)]
  obj@reductions[['HGT']] <- HGT_embedding
  for(j in 1:ncol(obj@reductions$HGT@cell.embeddings)){
    obj[[paste0('HGT_',j)]] <- obj@reductions$HGT@cell.embeddings[,paste0('HGT_',j)]
  }
  obj <- FindNeighbors(obj, reduction = "HGT",dims=1:nhid)
  obj <- FindClusters(obj, resolution = resolution)

  return(obj)
}


# Draw a heatmap of marker proteins or genes
# required packages: dsb, ComplexHeatmap
# input:
# obj: a seurat object obtained from HGT_cluster fuction
# marker: character of markers(the order of markers corresponding to ctorder)
# ctorder: the order of celltypes to display
# assays: RNA or ADT
# output:
# a heatmap of marker proteins or genes

MarkerHeatmap <- function(obj,marker,ctorder,assays){
  sortmarker <- 1:length(marker)
  names(sortmarker) <- marker
  if (assays == 'RNA'){
    rnamarker <- intersect(rownames(obj@assays$RNA@counts),marker)
    rnamarker <- rnamarker[order(sortmarker[rnamarker])]
    rnam <- AverageExpression(obj,assays = 'RNA', slot='counts',group.by = 'cell_type',features = rnamarker)
    rnam <- rnam[[1]][unique(rnamarker),ctorder]
    t <- DSBNormalizeProtein(rnam,rnam,use.isotype.control =FALSE,denoise.counts =FALSE)
  }else if(assays == 'ADT'){
    adtmarker <- intersect(rownames(obj@assays$ADT@counts),marker)
    adtmarker <- adtmarker[order(sortmarker[adtmarker])]
    adtm <- AverageExpression(obj,assays = 'ADT', slot='counts',group.by = 'cell_type',features = adtmarker)
    adtm <- adtm[[1]][unique(adtmarker),ctorder]
    t <- DSBNormalizeProtein(adtm,adtm,use.isotype.control =FALSE,denoise.counts =FALSE)
  }
  Heatmap(t(t),c("white","blue"),cluster_rows = F,cluster_columns = F,
          heatmap_legend_param = list( title = "normalized.expression",title_position = "leftcenter-rot" ))
}


# take a subset of obj and preprocess data again
# input:
# obj: a seurat object obtained from HGT_cluster fuction
# I: cell names subset
# output:
# obj: a subset of obj after preprocessing

subobject <- function(obj,I){
  obj <- obj[,I]
  DefaultAssay(obj) <- 'RNA'
  obj <- NormalizeData(obj) %>% FindVariableFeatures()
  all.genes <- rownames(obj)
  obj <- ScaleData(obj, features = all.genes)
  DefaultAssay(obj) <- 'ADT'
  VariableFeatures(obj) <- rownames(obj[["ADT"]])
  obj <- NormalizeData(obj, normalization.method = 'CLR', margin = 2) %>%
    ScaleData()
  return(obj)
}


# find features mostly associated with a given feature based on the expression of a type of features and two types of cells.
# input:
# feature: a given feature (gene or protein)
# ident: an expression matrix of a type of features and two types of cells.
# output:
# co: features mostly associated with the given feature
Findcofeatures <- function(feature,ident){
  y <- ident[feature,]
  co <- c()
  for (f in rownames(ident)){
    co <- rbind(co, c(f,cor(ident[f,],y,method = "spearman")))
    if (is.na(cor(ident[f,],y,method = "spearman"))){
      print(f)
    }
  }
  co <- co[order(co[,2],decreasing =TRUE),]
  print(which(co[,1]==feature))
  co <- co[which(co[,1]!=feature),]
  co <- co[which(co[,2]>0),]
  return(co)
}


# find all features mostly associated with top DE features based on the expression of two types of cells.
# input:
# rna.markers: the result of DE genes analysis between two cell types from FindMarkers function
# rna_ident: an RNA expression matrix of two types of cells after discarding zero expression genes.
# adt.markers: the result of DE proteins analysis between two cell types from FindMarkers function
# adt_ident: an ADT abundance matrix of two types of cells after discarding zero abundance proteins.
# output:
# allcofeatures: features mostly associated with the given features
Findallcofeatures <- function(rna.markers,rna_ident,adt.markers,adt_ident){
  allcofeatures <- list()
  options(scipen = 100)
  corna1 <- list()
  corna2 <- list()
  coadt1 <- list()
  coadt2 <- list()
  for (i in 1:10){
    feature <- rownames(rna.markers)[i]
    print(feature)
    co <- Findcofeatures(feature,rna_ident)
    co <- co[which(co[,1] %in% rownames(rna.markers)[1:10]==FALSE),]
    if (rna.markers[feature,"avg_log2FC"]>0){
      corna1[[feature]] <- co[1:5,1]
    }else{
      corna2[[feature]] <- co[1:5,1]
    }
  }
  for (i in 1:3){
    feature <- rownames(adt.markers)[i]
    print(feature)
    co <- Findcofeatures(feature,adt_ident)
    co <- co[which(co[,1] %in% rownames(adt.markers)[1:3]==FALSE),]
    if (adt.markers[feature,"avg_log2FC"]>0){
      coadt1[[feature]] <- co[1:5,1]
    }else{
      coadt2[[feature]] <- co[1:5,1]
    }
  }
  allcofeatures[['corna1']] <- corna1
  allcofeatures[['corna2']] <- corna2
  allcofeatures[['coadt1']] <- coadt1
  allcofeatures[['coadt2']] <- coadt2

  return(allcofeatures)
}


# calculate subset of obj under two cell types, the DE features and top associated features
# input:
# obj: a seurat object obtained from HGT_cluster fuction
# ident.1, ident.2: two cell types
# output:
# obj_co: a list of the subset of obj and selected features
subobjco <- function(obj, ident.1, ident.2){
  obj_co <- list()
  # DEG analysis between ident.1 and ident.2
  DefaultAssay(obj) <- 'RNA'
  rna.markers <-  FindMarkers(obj, ident.1 = ident.1 ,ident.2 = ident.2)
  DefaultAssay(obj) <- 'ADT'
  adt.markers <-  FindMarkers(obj, ident.1 = ident.1,ident.2 = ident.2)

  I0 <- colnames(obj)[which(obj@meta.data$cell_type %in% c(ident.1,ident.2))]
  obj0 <- subobject(obj,I0)

  # expression data for rnas and adts of two cell types
  rna_ident  <- obj0@assays$RNA@data
  adt_ident  <- obj0@assays$ADT@data
  rna_ident <- rna_ident[which(rowSums(rna_ident)!=0),]
  adt_ident <- adt_ident[which(rowSums(adt_ident)!=0),]

  cofeatures <- Findallcofeatures(rna.markers = rna.markers, rna_ident = rna_ident, adt.markers = adt.markers, adt_ident = adt_ident)
  corna1 <- cofeatures[['corna1']]
  corna2 <- cofeatures[['corna2']]
  coadt1 <- cofeatures[['coadt1']]
  coadt2 <- cofeatures[['coadt2']]
  obj_co[['cofeatures']] <- cofeatures
  obj_co[['obj']] <- obj0
  obj_co[['I']] <- I0
  obj_co[['rna.markers']] <- rna.markers
  obj_co[['adt.markers']] <- adt.markers
  return(obj_co)
}


# draw heatmap on correlation matrix between selected features
# required packages: RColorBrewer
# input:
# obj_co: a list obtained from subbmco function
# output:
# heatmap: heatmap on correlation matrix between selected features
#library(RColorBrewer)
DEfeaturesheatmap <- function(obj_co){
  obj0 <- obj_co$obj
  corna1 <- obj_co$cofeatures$corna1
  corna2 <- obj_co$cofeatures$corna2
  coadt1 <- obj_co$cofeatures$coadt1
  coadt2 <- obj_co$cofeatures$coadt2
  corna <- c(corna1,corna2)
  coadt <- c(coadt1,coadt2)
  # mm - correlations matrix as input of Heatmap
  m <-rbind(obj0@assays$RNA@scale.data[unique(c(unlist(corna), names(corna))),], obj0@assays$ADT@scale.data[unique(c(unlist(coadt), names(coadt))),])
  rownames(m) <- c(unique(c(unlist(corna), names(corna))),unique(c(unlist(coadt), names(coadt))))
  mm <- cor(t(m))
  # lab - the coressponding cel types of features
  lab <- list()
  for (f in rownames(m)){
    if (f %in% unique(c(unlist(corna1), names(corna1),unlist(coadt1), names(coadt1)))){
      lab[[f]] <- color1
    }else{
      lab[[f]] <- color2
    }
  }
  p <- heatmap(mm,RowSideColors=unlist(lab),col=colorRampPalette(c(rep("royalblue",1),'white','#F6DBDB',rep("firebrick3",2)))(56))
  return(p)
}


# draw a line representing a given feature
# input：
# f: a given feature
# dff: a dataframe containing expression of selected features and a dimension values of HGT embedding.
# color : color of a line
# output:
# draw a line representing a given feature
line <- function(f,dff,color){
  df <- data.frame(cbind(dff[,f],dff$index))
  model1=loess(X1 ~ X2,data=df,span=1 )
  model <- stats::predict(model1)
  lines(model, x=df$X2, col=color)
}


# draw the lines by a loess smoothing function based on the corresponding embedding and scaled gene expressions in cells
# input:
# obj: a seurat object obtained from HGT_cluster fuction
# obj_co: a list obtained from subbmco function
# HGT.dim: a dimension of HGT embedding
# color1, color2: two colors coressponding to two cell types
# output:
# a plot show the relationship between HGT embedding and feature expression
DEfeaturesloess <- function(obj, obj_co, HGT.dim = n, color1, color2){
  obj0 <- obj_co$obj
  I0 <- obj_co$I
  corna1 <- obj_co$cofeatures$corna1
  corna2 <- obj_co$cofeatures$corna2
  coadt1 <- obj_co$cofeatures$coadt1
  coadt2 <- obj_co$cofeatures$coadt2
  corna <- c(corna1,corna2)
  coadt <- c(coadt1,coadt2)

  m <-obj0@assays$RNA@scale.data[unique(c(unlist(corna), names(corna))),]
  df0 <- data.frame(t(m))
  colnames(df0) <- rownames(m)
  df0$index <-unlist(obj[[paste0('HGT_', HGT.dim)]][I0,])
  df0  <- df0[order(df0$index),]
  m <-obj0@assays$ADT@scale.data[unique(c(unlist(coadt), names(coadt))),]
  df1 <- data.frame(t(m))
  colnames(df1) <- rownames(m)
  df1$index <-unlist(obj[[paste0('HGT_', HGT.dim)]][I0,])
  df1  <- df1[order(df1$index),]

  rna.markers <- obj_co$rna.markers
  adt.markers <- obj_co$adt.markers
  par(pin = c(3,3))
  y=seq(-2,2,0.1)
  x=seq(min(df0$index),max(df0$index),length.out=length(y))
  #x=seq(-0.1,0.05,length.out=length(y))
  p <- plot(x,y,col="white",xlab = paste0('HGT_', HGT.dim), ylab="scaled experssion",type="l",main="Loess Smoothing and Prediction")
  for (j in 1:10){
    feature <- rownames(rna.markers)[j]
    if (rna.markers[feature,'avg_log2FC']>0){
      color <- color1
    }else{
      color <- color2
    }
    line (feature,df0,color)
    for (cof in unlist(corna[[feature]])[1:4]){
      line (cof,df0,color)
    }
  }
  for (j in 1:3){
    feature <- rownames(adt.markers)[j]
    if (adt.markers[feature,'avg_log2FC']>0){
      color <- color1
    }else{
      color <- color2
    }
    line (feature,df1,color)
    for (cof in unlist(coadt[[feature]])[1:4]){
      line (cof,df1,color)
    }
  }
  return(p)
}


# required package -- reticulate
# input:
# GAS: the spliced and normalized matrix obtained from CLR function
# result_dir: The address for storing the models and optimization results(Type:str)
# epoch:(Type:int)
# lr: learning rate(Type:float)
# n_hid: Number of hidden dimension(Type:int)
# n_heads: Number of attention head(Type:int)
# cuda: 0 use GPU0 else cpu(Type:int)
# data_type: 'CITE', 'scRNA_scATAC', or 'multipleRNA'
# envPath: The address for environment to use if use.env is TRUE(Type:str)
# output:
# HGT_result: a list containing requried results of HGT model as follows:
# parameters: given parameters from user --epoch, lr, n_hid, n_heads, cuda
# cell_hgt_matrix: cell embedding matrix
# feature_hgt_matrix : gene embedding matrix and protein embedding matrix when data_type is 'CITE';
# attention: attention meassage for features and cells
# data_type: 'CITE', 'scRNA_scATAC', or 'multipleRNA'
# result_dir: The address for storing the models and optimization results
# GAS: the spliced and normalized matrix obtained from CLR function

run_HGT <- function(GAS,result_dir,data_type,envPath=NULL,lr=NULL, epoch=NULL, n_hid=NULL, n_heads=NULL,cuda=0){
  if (data_type == 'CITE') {
    if (is.null(lr)){lr = 0.1}
    if (is.null(epoch)){epoch = 50}
    if (is.null(n_hid)){n_hid = 104}
    if (is.null(n_heads)){n_heads = 13}
  }
  if (data_type == 'scRNA_scATAC') {
    if (is.null(lr)){lr = 0.1}
    if (is.null(epoch)){epoch = 100}
    if (is.null(n_hid)){n_hid = 128}
    if (is.null(n_heads)){n_heads = 16}
  }
  if (data_type == 'multipleRNA') {
    if (is.null(lr)){lr = 0.1}
    if (is.null(epoch)){epoch = 100}
    if (is.null(n_hid)){n_hid = 104}
    if (is.null(n_heads)){n_heads = 13}
  }
  print(epoch)
  if (!is.null(envPath)){use_condaenv(envPath)}
  list_in <- assign("list_in", list(lr=lr, epoch=epoch, n_hid=n_hid, n_heads=n_heads, result_dir=result_dir, cuda=cuda, data_type=data_type, cell_gene=GAS, gene_name=rownames(GAS), cell_name=colnames(GAS)), envir = .GlobalEnv)
  source_python('./arg.py')
  cell_hgt_matrix <- py$cell_matrix
  gene_hgt_matrix <- py$gene_matrix
  attention <- py$df2
  rownames(cell_hgt_matrix) <- list_in$cell_name
  rownames(gene_hgt_matrix) <- list_in$gene_name
  HGT_result <- list()
  HGT_result[['parameters']] <- data.frame(lr,epoch,n_hid,n_heads,cuda)
  HGT_result[['GAS']] <- GAS
  HGT_result[['cell_hgt_matrix']] <- cell_hgt_matrix
  HGT_result[['feature_hgt_matrix']] <- gene_hgt_matrix
  HGT_result[['attention']] <- attention
  HGT_result[['result_dir']] <- result_dir
  HGT_result[['data_type']] <- data_type
  return(HGT_result)
}
