## ----------------------------------------------------------------------------------------------------------------
library(bluster)
library(GenomeInfoDb)
library(igraph)
library(cluster)
library(GenomicRanges)
library(IRanges)
library(ggplot2)
library(Matrix)
library(dplyr)
library(data.table)
library(BioQC)
library(tinytex)
library(tidyverse)
library(JASPAR2020)
library(rtracklayer)
library(motifmatchr)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(EnsDb.Mmusculus.v79)
library(Seurat)
library(scater)
library(Signac)
library(MAESTRO)
library(reticulate)
## ----------------------------------------------------------------------------------------------------------------

#knitr::purl('deepmaps.rmd')



## ----------------------------------------------------------------------------------------------------------------
# Read 10X data
##' Read matched scRNA + scATAC data from H5 file
#input:
# 1 - h5Path: the path of h5 file
# 2 - min_cell: the peak / gene will be removed if the value in the gene / peak with more than min_cell cell is equal to zero
#output:
#a seurat object
Read10Xdata <-
  function(h5Path,
           annoObj = NULL,
           fragmentsPath = NULL,
           hintPath = NULL,
           min_cell = 0.01) {
    inputdata.10x <- Read10X_h5(h5Path)
    rna_counts <- inputdata.10x$`Gene Expression`
    atac_counts <- inputdata.10x$Peaks
    grange.counts <-
      StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
    grange.use <-
      seqnames(grange.counts) %in% standardChromosomes(grange.counts)
    atac_counts <- atac_counts[as.vector(grange.use),]
    chrom_assay <- CreateChromatinAssay(
      counts = atac_counts,
      sep = c(":", "-"),
      min.cells = ncol(atac_counts) * min_cell,
      fragments = fragmentsPath,
      annotation = anno
      #min.feature = 300,
    )
    tmp_obj <- CreateSeuratObject(counts = chrom_assay,
                                  assay = "ATAC",)
    exp_assay <-
      CreateAssayObject(counts = rna_counts,
                        min.cells = ncol(rna_counts) * min_cell)
    tmp_obj[["RNA"]] <- exp_assay
    DefaultAssay(tmp_obj) <- "RNA"
    tmp_obj[["percent.mt"]] <-
      PercentageFeatureSet(tmp_obj, pattern = "^MT-")
    return (tmp_obj)
  }


# Read data with matrix format
##' Read matched scRNA+scATAC data from matrix format
#input:
# 1 - rna_matrix: an expression matrix with gene * cell
# 2 - atac_matrix: an accessibility matrix with peak * cell
# 3 - min_cell: the peak/gene will be removed if the value in the gene/peak with more than min_cell cell is equal to zero
#output:
#a seurat object
readmatrix <- function(rna_matrix, atac_matrix, min_cell = 0.1) {
  rna_matrix <-
    rna_matrix[, intersect(colnames(atac_matrix), colnames(rna_matrix))]
  atac_matrix <-
    atac_matrix[, intersect(colnames(atac_matrix), colnames(rna_matrix))]
  grange.counts <-
    StringToGRanges(rownames(atac_matrix), sep = c(":", "-"))
  grange.use <-
    seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_matrix <- atac_matrix[as.vector(grange.use),]
  atac_matrix <-
    atac_matrix[lengths(strsplit(gsub(":", "-", rownames(atac_matrix)) , split = "-")) ==
                  3, ]
  rna_matrix <- rna_matrix[unique(rownames(rna_matrix)), ]
  #cell_type<-rna_cell$V7[-1][grepl("*RNA*",rna_gene$V2[-1])]
  #min_cell=0.01
  chrom_assay <- CreateChromatinAssay(
    counts = atac_matrix,
    sep = c(":", "-"),
    #genome = annota,
    #fragments = fragments,
    min.cells = ncol(atac_matrix) * min_cell,
    #min.feature = 300,
    #annotation = annotations
  )
  obj <- CreateSeuratObject(counts = chrom_assay,
                            assay = "ATAC")
  exp_assay <-
    CreateAssayObject(counts = rna_matrix,
                      min.cells = ncol(rna_matrix) * min_cell)
  obj[["RNA"]] <- exp_assay
  DefaultAssay(obj) <- "RNA"
  obj[["percent.mt"]] <-
    PercentageFeatureSet(obj, pattern = "^mt-")
  return(obj)
}


## ----------------------------------------------------------------------------------------------------------------
# Filter abnormal cells
##' Filter abnormal cells
#input:
# 1 - obj: a seurat object
# 2 - nmad: a numeric scalar, specifying the minimum number of MADs away from median required for a value to be called an outlier
# output:
# a seurat object

filterCell <- function(obj, nmad = 3) {
  atac <- isOutlier(obj$nCount_ATAC,
                    nmads = nmad,
                    log = F,
                    type = "both")
  rna <- isOutlier(obj$nCount_RNA,
                   nmads = nmad,
                   log = F,
                   type = "both")
  
  mito <- isOutlier(obj$percent.mt,
                    nmads = nmad,
                    log = F,
                    type = "both")
  obj <-
    AddMetaData(obj, atac, col.name = "atac")
  obj <-
    AddMetaData(obj, rna, col.name = "rna")
  obj <-
    AddMetaData(obj, mito, col.name = "mito")
  obj <- subset(x = obj,
                subset = atac == F &
                  rna == F &
                  mito == F)
  return(obj)
  
}


## ----------------------------------------------------------------------------------------------------------------
##' Calculate gene active score matrix
# input:
# 1 - peak_count_matrix: a peak_count matrix from scATAC-seq with peak * cell which return from filterCell function
# 2 - organism: species type GRCh38 / GRCm38
# output:
# a gene * peak matrix, the elements represent the regulatory potential for peak to gene

CalGenePeakScore <-
  function(peak_count_matrix, organism = "GRCh38") {
    pbmc_peak <- peak_count_matrix
    n <- nrow(pbmc_peak)
    dia <- diag(n)
    rownames(dia) <- rownames(pbmc_peak)
    colnames(dia) <- 1:ncol(dia)
    gene_peak <-
      ATACCalculateGenescore(dia,
                             organism = organism,
                             decaydistance = 10000,
                             model = "Enhanced")
    colnames(gene_peak) <- rownames(peak_count_matrix)
    return (gene_peak)
  }


## ----------------------------------------------------------------------------------------------------------------
##' Calculate gene active score matrix
#input:
# 1 - ATAC_gene_peak: a matrix with gene * peak which return from CalGenePeakScore fucntion
# 2 - obj: a seurat object after data preprocessing which return from filterCell function
# 3 - method: the method to integrate scRNA-seq and scATAC-seq velo (velocity) / WNN (weighted nearest neighbor)
# 4 - veloPath: if use velocity method, the veloPath should be provided
#output:
# GAS matrix with gene * peak, the elements represent the gene activity score in each cell
# a gene * peak matrix, the elements represent the regulatory potential for peak to gene

calculate_GAS <-
  function(ATAC_gene_peak,
           obj,
           method = "velo",
           veloPath = NULL) {
    peak_count <- obj@assays$ATAC@counts
    gene_count <- obj@assays$RNA@counts
    peak_count[peak_count > 0] = 1
    WA <- ATAC_gene_peak %*% peak_count
    WA <- WA[which(rowSums(as.matrix(WA)) > 0), ]
    gene_count <-
      gene_count[which(rowSums(as.matrix(gene_count)) > 0), ]
    commongene <-
      intersect(x = rownames(WA), y = rownames(gene_count))
    WA <- as.matrix(WA)
    WA <- WA[commongene, ]
    gene_count <- gene_count[commongene, ]
    gene_rowsum <- rowSums(gene_count)
    peak_rowsum <- rowSums(WA)
    norm_gene_count <- gene_count / rowSums(gene_count)
    norm_WBinary <- WA / rowSums(WA)
    #norm_gene_count<-NormalizeData(CreateSeuratObject(counts = gene_count))$ RNA@data
    gene_count <- norm_gene_count
    #norm_WBinary<-NormalizeData(CreateSeuratObject(counts = WA))$RNA@data
    peak_count <- norm_WBinary
    print(str(peak_count))
    if (method == "velo") {
      velo <- read.csv(veloPath, header = TRUE)
      ##remove duplicated rows with same gene
      if (length(which(duplicated.default(velo[, 1]))) > 0) {
        velo <- velo[-which(duplicated.default(velo[, 1]) == T),]
      }
      
      ##Filter rows if gene name equals to NA
      if (length(which(is.na(velo[, 1]))) > 0) {
        velo <- velo[-which(is.na(velo[, 1])),]
      }
      rownames(velo) <- velo[, 1]
      velo <- velo[,-1]
      velo <- as.matrix(velo)
      #colnames(velo) <- gsub("-", ".", colnames(velo))
      colnames(velo) <- gsub("\\.", "-", colnames(velo))
      rna <- gene_count
      atac <- peak_count
      velo <-
        velo[intersect(rownames(rna), rownames(velo)), intersect(colnames(rna), colnames(velo))]
      rna <-
        rna[intersect(rownames(rna), rownames(velo)), intersect(colnames(rna), colnames(velo))]
      atac <-
        atac[intersect(rownames(rna), rownames(velo)), intersect(colnames(rna), colnames(velo))]
      gene_rowsum <- gene_rowsum[rownames(rna)]
      peak_rowsum <- peak_rowsum[rownames(rna)]
      
      genes <- dim(velo)[1]
      cells <- dim(velo)[2]
      # rank matrix
      rank_cell <- velo
      rank_gene <- velo
      rank_cell <- apply(velo, 2, rank)
      rank_gene <- t(apply(velo, 1, rank))
      
      rank_cell[velo > 0] = genes - rank_cell[velo > 0]
      rank_cell[velo < 0] = rank_cell[velo < 0] - 1
      rank_cell[velo == 0] = 0
      rank_gene[velo > 0] = cells - rank_gene[velo > 0]
      rank_gene[velo < 0] = rank_gene[velo < 0] - 1
      rank_gene[velo == 0] = 0
      
      # number of positive/negative for each gene/cell
      cell_posi_num <- colSums(velo > 0)
      cell_nega_num <- colSums(velo < 0)
      gene_posi_num <- rowSums(velo > 0)
      gene_nega_num <- rowSums(velo < 0)
      
      # weights
      weights <-
        ((rank_cell ^ 2 + rank_gene ^ 2) / ((t((t(velo > 0)) * cell_posi_num + (t(velo <
                                                                                    0)) * cell_nega_num
        )) ^ 2 + ((velo > 0) * gene_posi_num + (velo < 0) * gene_nega_num
        ) ^ 2 + (velo == 0))) ^ 0.5
      weights[velo < 0] = weights[velo < 0] * (-1)
      
      # GAS
      GAS <-
        rna * gene_rowsum + ((1 + weights) * atac) * ((1 + weights) *
                                                        peak_rowsum)
    }
    if (method == "wnn") {
      obj <- obj[, colnames(gene_count)]
      DefaultAssay(obj) <- "RNA"
      obj <-
        FindVariableFeatures(obj,
                             selection.method = "vst",
                             nfeatures = 2000)
      obj <-
        ScaleData(obj, features = VariableFeatures(obj))
      obj <- RunPCA(obj)
      # ATAC analysis
      # We exclude the first dimension as this is typically correlated with sequencing depth
      DefaultAssay(obj) <- "ATAC"
      obj <- RunTFIDF(obj)
      obj <- FindTopFeatures(obj, min.cutoff = 'q0')
      obj <- RunSVD(obj)
      obj <-
        RunUMAP(
          obj,
          reduction = 'lsi',
          dims = 2:50,
          reduction.name = "umap.atac",
          reduction.key = "atacUMAP_"
        )
      obj <-
        FindMultiModalNeighbors(obj,
                                reduction.list = list("pca", "lsi"),
                                dims.list = list(1:50, 2:50))
      GAS <-
        gene_count * gene_rowsum * obj$RNA.weight + peak_count * obj$ATAC.weight *
        peak_rowsum
    }
    
    return(GAS)
  }


## ----------------------------------------------------------------------------------------------------------------
# CT active gene modules calculation
# input:
# 1 - GAS: a gene active matrix with gene * cell which return from calculate_GAS function
# 2 - cell_hgt_matrixPath: the path of cell-embedding matrix which otains from HGT function
# 3 - attPath: the path of attention matrix with gene-cell * head which obtain from HGT function
#output:
# 1 - co (variable 1): a biological gene module. a list with name CT-i and active gene list in CT-i
# 2 - graph.out (variable 2): the cluster result with a factor format

get_gene_module <-
  function(obj, cell_hgt_matrix, att, GAS, cutoff = 1.6) {
    graph.out <- Idents(obj)
    nhead <- ncol(att)
    gene_name <- rownames(GAS)[att$gene + 1]
    cell_name <- colnames(GAS)[att$cell + 1]
    
    att$ct <- graph.out[cell_name]
    att$gene_name <- gene_name
    att$cell_name <- cell_name
    mod <- function(x) {
      return(sqrt(sum(c(x ^ 2))))
    }
    nor <- function(x) {
      return((x - min(x)) / (max(x) - min(x)))
    }
    
    att[, 4:nhead] <- nor(att[, 4:nhead])
    attention <-
      aggregate(x = as.list(att[, 4:nhead]),
                by = list(att$ct, att$gene_name),
                mean)
    #att[,4:nhead]<-1-att[,4:nhead]
    weight <- apply(att[, 4:nhead], 1, mod)
    df <-
      data.frame(
        'node1' = att$gene_name,
        'node2' = att$cell_name,
        'weight' = weight,
        'ct' = att$ct
      )
    attention <-
      aggregate(x = df$weight, by = list(df$ct, df$node1), mean)
    co <- list()
    for (i in (0:(length(unique(att$ct)) - 1))) {
      t <-
        mean(attention[attention$Group.1 == i,]$x) + 1.6 * sd(attention[attention$Group.1 ==
                                                                          i,]$x)
      co[[paste('ct', i, sep = "_")]] <-
        attention[attention$Group.1 == i,]$Group.2[attention[attention$Group.1 ==
                                                               i,]$x > t]
    }
    m <- list()
    m[[1]] <- co
    m[[2]] <- graph.out
    return (m)
    
  }


## ----------------------------------------------------------------------------------------------------------------
# gene module save
#input:
# 1 - co: the active gene module from get_gene_module function
# 2 - lisa_path: the path of active gene module to save
#result
# write gene module to the lisa_path

write_GM <- function(co, lisa_path) {
  if (length(dir(path = lisa_path, pattern = ".csv")) >
      0) {
    system(paste0("rm ", lisa_path, "*.csv "))
    
  }
  if (length(dir(path = lisa_path, pattern = ".txt")) >
      0) {
    system(paste0("rm ", lisa_path, "*.txt "))
  }
  
  for (j in (1:length(co))) {
    if (length(unique(co[[j]])) < 20 |
        length(unique(co[[j]])) > 20000) {
      next
    } else{
      ct <- unlist(strsplit(names(co[j]), split = "_"))[1]
      
      write.table(
        co[[j]],
        paste0(lisa_path, names(co[j]), ".txt"),
        quote = F,
        sep = "\t",
        row.names = F,
        col.names = F
      )
    }
  }
}


## ----------------------------------------------------------------------------------------------------------------
# Filter gene with no accessible peak in promoter
# input:
# 1 - obj: a seurat object which return from filterCell function
# 2 - gene_peak: a matrix with gene * peak from scATAC-seq which return from filterCell function
# 3 - GAS: the GAS matrix with gene * cell which return calculate_GAS function
# 4 - species: human / mouse
#output:
# a matrix with gene * peak. The gene with no accessible peak will be removed

AccPromoter <- function(obj, gene_peak, GAS, species = "hg38") {
  peak_cell <- obj@assays$ATAC@counts
  if (species == "hg38") {
    gene.ranges <- genes(EnsDb.Hsapiens.v86)
  } else{
    gene.ranges <- genes(EnsDb.Mmusculus.v79)
  }
  
  gene.use <-
    seqnames(gene.ranges) %in% standardChromosomes(gene.ranges)[standardChromosomes(gene.ranges) !=
                                                                  "MT"]
  gene.ranges <- gene.ranges[as.vector(gene.use)]
  gene.ranges <-
    gene.ranges[gene.ranges$gene_name %in% rownames(GAS)]
  genebodyandpromoter.coords <-
    Extend(x = gene.ranges,
           upstream = 2000,
           downstream = 0)
  #str(genebodyandpromoter.coords)
  x <- as.data.frame(genebodyandpromoter.coords@ranges)
  peaks <-
    GRanges(
      seqnames = paste("chr", genebodyandpromoter.coords@seqnames, sep = ""),
      ranges = IRanges(start = , x$start,
                       width = x$width)
    )
  
  peak_name <-
    colnames(gene_peak)[lengths(strsplit(gsub(":", "-", colnames(gene_peak)) , split = "-")) ==
                          3]
  peak_name <-
    do.call(what = rbind, strsplit(gsub(":", "-", peak_name) , split = "-"))
  peak_name <- as.data.frame(peak_name)
  names(peak_name) <- c("chromosome", 'start', 'end')
  peak_name <- GenomicRanges::makeGRangesFromDataFrame(peak_name)
  #str(peaks)
  over <- findOverlaps(peak_name, peaks)
  str(over)
  promoter_gene <-
    genebodyandpromoter.coords$gene_name[unique(over@to)]
  str(promoter_gene)
  gene_peak <- gene_peak[promoter_gene, ]
  
  return(gene_peak)
}


## ----------------------------------------------------------------------------------------------------------------
##' infer ct active regulons
# input:
# 1 - GAS: the GAS matrix with gene * cell which return calculate_GAS function
# 2 - co: a list of bio network which reture from gene_ function
# 3 - gene_peak_pro: the matrix with gene * peak which return AccPromoter function
# 4 - species: human / mouse (human = "hg38", mouse = "mm10" )
# 5 - humanPath: if species == human, the TF binding RData absolute path of hg38 should be provided
# 6 - mousePath: if species == mouse, the TF binding RData absolute path of mm10 should be provided
#output:
# 1 - BA_score a TF binding affinity matrix with TF * peak, the elements in the matrix is the binding power of TF to peak
# 2 - ct_regulon: candidate cell type active regulon
# 3 - TFinGAS: TRUE / FALSE, if number of the intersection of candidate TF from LISA and gene in GAS > 50 TFinGAS will be true, else it will be false

Calregulon <-
  function(GAS,
           co,
           gene_peak_pro,
           species = "hg38",
           jaspar_path = "/scratch/deepmaps/jaspar",
           lisa_path = "/home/wan268/hgt/RNA_ATAC/lymph_14k/") {
    if (species == "hg38") {
      tfbs_df <- qs::qread(paste0(jaspar_path, "hg38_lisa_500.qsave"))
      tfbs_df <- tfbs_df[1:(nrow(tfbs_df) - 1),]
    }
    else {
      tfbs_df <- readRDS("/fs/ess/PCON0022/wxy/mm10.rds")
      tfbs_df[tfbs_df$V6 > 500, ]
    }
    
    BA_score <-
      matrix(0, ncol(gene_peak_pro), length(unique(tfbs_df$V4)))
    colnames(BA_score) <- unique(tfbs_df$V4)
    rownames(BA_score) <- colnames(gene_peak_pro)
    gene_TF <-
      matrix(0, nrow(gene_peak_pro), length(unique(tfbs_df$V4)))
    colnames(gene_TF) <- unique(tfbs_df$V4)
    rownames(gene_TF) <- rownames(gene_peak_pro)
    
    peak <- tfbs_df[, 1:3]
    colnames(peak) <- c("chromosome", 'start', 'end')
    peak <- GenomicRanges::makeGRangesFromDataFrame(peak)
    
    ct_subregulon <- list()
    ct_regulon <- list()
    coexp_tf <- list()
    No <- 1
    for (i in (1:length(co))) {
      if (length(co[[i]]) > 0) {
        co[[i]] <- intersect(co[[i]], rownames(gene_peak_pro))
        a <- which(gene_peak_pro[co[[i]], ] > 0, arr.ind = T)
        op <- colnames(gene_peak_pro)[unname(a[, 'col'])]
        peak_name <-
          op[lengths(strsplit(gsub(":", "-", op) , split = "-")) == 3]
        peak_name <-
          do.call(what = rbind, strsplit(gsub(":", "-", peak_name) , split = "-"))
        peak_name <- as.data.frame(peak_name)
        names(peak_name) <- c("chromosome", 'start', 'end')
        peak_name <-
          GenomicRanges::makeGRangesFromDataFrame(peak_name)
        over <- findOverlaps(peak_name, peak)
        #print(i)
        p <- op[over@from]
        pp <- tfbs_df$V5[over@to] / 100
        df <- data.frame(p, pp, tfbs_df$V4[over@to])
        hh <- df[!duplicated(df[, -2]), ]
        for (k1 in (1:nrow(hh))) {
          BA_score[hh[k1, ]$p, hh[k1, ]$tfbs_df.V4.over.to.] <- hh[k1, ]$pp
        }
        
        gene_TF <- gene_peak_pro %*% BA_score
        TF <- unique(tfbs_df$V4[over@to])
        
        if (length(co[[i]]) < 20000 & length(co[[i]]) > 20) {
          tf <-
            read.csv(paste(lisa_path,
                           names(co[i]),
                           ".txt.csv",
                           sep = ""))
          tf_pval_0.05 <- unique(tf[, 3][tf$summary_p_value < 0.05])
          TF <- intersect(unique(TF), tf_pval_0.05)
        }
        #print(length(TF))
        
        #gene_TF[co[[i]],pp]<-gene_peak_pro[co[[i]],p] %*% BA_score[p,TF]
        coexp_tf[[names(co[i])]] <- TF
        #print(length(TF))
        if (length(TF) > 0) {
          for (k in 1:length(TF)) {
            if (TF[k] %in% rownames(GAS)) {
              a <- unlist(strsplit(names(co[i]), "_"))
              a <- paste0(a[1], a[2])
              h <- paste(TF[k], a, sep = "_")
              ct_subregulon[[h]] <-
                co[[i]][gene_TF[co[[i]], TF[k]] > 0]
            }
          }
        }
      }
    }
    
    m <- list()
    m[[1]] <- BA_score
    ct_subregulon <- ct_subregulon[lengths(ct_subregulon) > 10]
    m[[2]] <- ct_subregulon
    
    return(m)
  }


## ----------------------------------------------------------------------------------------------------------------
# combine same TF
#input:
# 1 - gene_peak_pro: a matrix with gene * peak. The gene with no accessible peak will be removed which return from AccPromoter function
# 2 - BA_score: a TF binding affinity matrix with TF * peak, the elements in the matrix is the binding power of TF to peak which returen from Calregulon function
#output:
# 1 - peak_TF: a matrix with peak * TF without repeat TF

uni <- function(gene_peak_pro, BA_score) {
  gene_TF <- gene_peak_pro %*% BA_score
  rownames(gene_TF) <- rownames(gene_peak_pro)
  mat <-
    matrix(0, nrow = length(rownames(gene_TF)), ncol = length(unique(colnames(gene_TF))))#peak_TF
  mat1 <-
    matrix(0, nrow = length(rownames(BA_score)), ncol = length(unique(colnames(BA_score))))#gene_TF
  rownames(mat1) <- rownames(BA_score)
  colnames(mat1) <- unique(colnames(BA_score))
  #mat1: peak_TF score
  for (x in unique(colnames(BA_score))) {
    if (is.null(nrow(gene_TF[, colnames(BA_score) == x]))) {
      #mat<-rbind(mat, unlist(gene_peak_matrix[rownames(gene_peak_matrix)==x,]))
      mat1[, x] <-
        unname(unlist(BA_score[, colnames(BA_score) == x]))
    } else{
      #mat<-rbind(mat,unlist(colSums(gene_peak_matrix[rownames(gene_peak_matrix)==x,])))
      mat1[, x] <-
        unname(unlist(rowSums(BA_score[, colnames(BA_score) == x])))
    }
  }
  
  
  
  rownames(mat) <- rownames(gene_TF)
  colnames(mat) <- unique(colnames(gene_TF))
  
  for (x in unique(colnames(gene_TF))) {
    if (is.null(nrow(gene_TF[, colnames(gene_TF) == x]))) {
      print("111")
      #mat<-rbind(mat, unlist(gene_peak_matrix[rownames(gene_peak_matrix)==x,]))
      mat[, x] <- unname(unlist(gene_TF[, colnames(gene_TF) == x]))
    } else{
      #mat<-rbind(mat,unlist(colSums(gene_peak_matrix[rownames(gene_peak_matrix)==x,])))
      mat[, x] <-
        unname(unlist(rowSums(gene_TF[, colnames(gene_TF) == x])))
    }
  }
  m = list()
  m[[1]] <- mat
  m[[2]] <- mat1
  return (m)
}



## ----------------------------------------------------------------------------------------------------------------
# Regulatory Intensive (RI) score in cell level
# input:
# 1 - obj: a seurat object which return from filterCell function
# 2 - ct_regulon: cell type active regulon
# 3 - GAS: the GAS matrix with gene * cell which return calculate_GAS function
# 4 - gene_peak_pro: the matrix with gene * peak which return from AccPromoter function
# 5 - peak_TF: the matrix with peak * TF which return from uni function
# 6 - graph.out: a factor variable. The predict cell cluster which return from get_gene_module function
#output:
# 1 - RI_C: a regulatory intensive matrix with TF-gene pair * cell, the element means the intensity of TF to gene in each cell

RI_cell <-
  function(obj,
           ct_regulon,
           GAS,
           gene_peak_pro,
           peak_TF,
           graph.out) {
    v <- vector()
    for (i in (1:length(ct_regulon))) {
      v <-
        append(v, paste(unlist(strsplit(
          names(ct_regulon[i]), "_"
        ))[1], ct_regulon[[i]], sep = "_"))
    }
    peak_cell <- obj@assays$ATAC@counts
    peak_cell[peak_cell > 0] = 1
    peak_cell <- (peak_cell[, colnames(GAS)])
    #graph.out<-Idents(obj)
    #TG_cell <- foreach(i=1:length(ct_regulon), .packages='Matrix',.combine='c') %dopar%{
    TG_cell <- matrix(0, length(unique(v)), length(graph.out))
    rownames(TG_cell) <- unique(v)
    t <-
      unlist(strsplit(unique(v), "_"))[seq(1, length(unlist(strsplit(unique(
        v
      ), "_"))), 2)]
    g <-
      unlist(strsplit(unique(v), "_"))[seq(2, length(unlist(strsplit(unique(
        v
      ), "_"))), 2)]
    t <- unique(t)
    g <- unique(g)
    t1 <-
      unlist(strsplit(unique(v), "_"))[seq(1, length(unlist(strsplit(unique(
        v
      ), "_"))), 2)]
    g1 <-
      unlist(strsplit(unique(v), "_"))[seq(2, length(unlist(strsplit(unique(
        v
      ), "_"))), 2)]
    #TfGene_cell <- foreach(i=1:length(graph.out), .packages='Matrix',.combine='c') %dopar%
    #  {
    gene_peak <- gene_peak_pro
    for (j in (1:length(graph.out))) {
      hhh <-
        (peak_cell[, j] * gene_peak[g, ]) %*% (peak_cell[, j] * peak_TF[, t])
      #hhh<-gene_peak[g,] %*%peak_TF[,t]
      bb = data.frame('gene' = g1, 'tf' = t1)
      bb$tf = as.factor(bb$tf)
      bb$gene = as.factor(bb$gene)
      levels(bb$gene) = 1:length(levels(bb$gene))
      levels(bb$tf) = 1:length(levels(bb$tf))
      bb <- as.matrix(bb)
      bb <- apply(bb, 2, as.numeric)
      TG_cell[, j] <- hhh[bb]
      #print(j)
    }
    return(TG_cell)
  }

# calculate regulon active score in cell level / cell type level
# input:
# 1 - RI_C: a regulatory intensive matrix with TF-gene pair * cell, the element means the intensity of TF to gene in each cell
# 2 - ct_regulon: cell type active regulon
# 3 - graph.out: a factor variable. The cell cluster which return from get_gene_module function
# 4 - TFinGAS: TRUE / FALSE, if number of the intersection of candidate TF from LISA and gene in GAS > 50 TFinGAS will be true, else it will be false which return from Calregulon function
#output:
# 1 - RAS: regulon active score in cell type level
# 2 - RI_CT: regulatory intensive score in cell type level
# 3 - ct-regulon: cell type active regulon
# 4 - RAS_C: a matrix regulon-CT * cell, regulatory active score in cell level

calRAS <- function(RI_C, ct_regulon, graph.out) {
  a <- unlist(strsplit(names(ct_regulon), "_"))
  CT <- unique(a[seq(2, length(a), 2)])
  TF <- unique(a[seq(1, length(a), 2)])
  RI_CT <- matrix(0, nrow(RI_C), length(CT))
  rownames(RI_CT) <- rownames(RI_C)
  colnames(RI_CT) <- CT
  for (i in CT) {
    a <- graph.out == as.numeric(substring(i, 3, nchar(i)))
    RI_CT[, i] <- rowSums((RI_C[, a])) / length(graph.out[a])
  }
  g <-
    unlist(strsplit(rownames(RI_C), "_"))[seq(2, length(unlist(strsplit(
      rownames(RI_C), "_"
    ))), 2)]
  RAS_E <- RI_C * as.matrix(GAS[g, ])
  colnames(RAS_E) <- colnames(RI_C)
  #RAS_C TF(regulon)*cell
  RAS_C <-
    as(matrix(0, nrow = length(TF), ncol = length(graph.out)), "sparseMatrix")
  
  rownames(RAS_C) <- TF
  colnames(RAS_C) <- names(graph.out)
  for (i in (1:length(ct_regulon))) {
    tf <- unlist(strsplit(names(ct_regulon[i]), "_"))[1]
    #ct<-unlist(strsplit(names(ct_regulon[i]),"_"))[2]
    RAS_C[tf, ] <-
      colMeans(RAS_E[paste(tf, ct_regulon[[i]], sep = "_"), ])
  }
  
  #RAS TF*CT
  RAS <-
    as(matrix(0, nrow = length(TF), ncol = length(CT)), "sparseMatrix")
  rownames(RAS) <- TF
  colnames(RAS) <- CT
  for (i in (1:length(ct_regulon))) {
    tf <- unlist(strsplit(names(ct_regulon[i]), "_"))[1]
    ct <- unlist(strsplit(names(ct_regulon[i]), "_"))[2]
    if (sum(GAS[tf, graph.out == substring(ct, 3, nchar(ct))]) == 0) {
      RAS[tf, ct] <- 0
    }
    else{
      RAS[tf, ct] <-
        mean(RAS_E[paste(tf, ct_regulon[[i]], sep = "_"), graph.out == substring(ct, 3, nchar(ct))])
      
    }
  }
  ct_re <- list()
  for (i in (1:length(ct_regulon))) {
    tf <- unlist(strsplit(names(ct_regulon[i]), "_"))[1]
    ct <- unlist(strsplit(names(ct_regulon[i]), "_"))[2]
    re <-
      names(RI_CT[paste(tf, ct_regulon[[i]], sep = "_"), ct][RI_CT[paste(tf, ct_regulon[[i]], sep =
                                                                           "_"), ct] > 0])
    if (length(re) > 0) {
      ct_re[[names(ct_regulon[i])]] <-
        unlist(strsplit(re, "_"))[seq(2, length(unlist(strsplit(re, "_"))), 2)]
    }
  }
  ct_re <- ct_re[lengths(ct_re) > 10]
  
  m <- list()
  m[[1]] <- RAS
  m[[2]] <- RI_CT[rowSums(RI_CT) > 0, ]
  m[[3]] <- ct_re
  m[[4]] <- RAS_C
  return (m)
}

# calculate regulon active score in cell level / cell type level
# input:
# 1 - ct_regulon: a matrix with gene * peak from scATAC-seq which return from get_gene_module function
# 2 - graph.out: a factor variable. The predict cell cluster which return from get_gene_module function
#output:
# 1 - RAS_C1: a matrix with regulon-CT * cell. Regulon active score in cell type level

CalRAS2 <- function(ct_regulon, graph.out) {
  RAS_C1 <-
    as(matrix(0, nrow = length(ct_regulon), ncol = length(graph.out)), "sparseMatrix")
  rownames(RAS_C1) <- names(ct_regulon)
  colnames(RAS_C1) <- names(graph.out)
  for (i in (1:length(ct_regulon))) {
    tf <- unlist(strsplit(names(ct_regulon[i]), "_"))[1]
    ct <- unlist(strsplit(names(ct_regulon[i]), "_"))[2]
    g <- ct_regulon[[i]]
    RAS_C1[names(ct_regulon[i]), ] <-
      colMeans(RI_C[paste(tf, g, sep = "_"), ])
  }
  return(RAS_C1)
}




## ----------------------------------------------------------------------------------------------------------------
#find master TF and master gene, build a cell type gene regulatory network
#input:
# 1 - ct_regulon: cell type active regulon
# 2 - RI_CT: regulatory intensive score in cell type level which return from calRAS funciton
#output:
# 1 - TF_cen: TF centrality score in each CT
# 2 - gene_cen: gene centrality score in each CT
# 3 - network: adjcent matrix in a CT

masterFac <- function(ct_regulon, RI_CT) {
  TF_CT <- unlist(strsplit(names(ct_regulon), split = "_"))
  TF <- unique(TF_CT[seq(1, length(TF_CT), 2)])
  CT <- unique(TF_CT[seq(2, length(TF_CT), 2)])
  TF_R <- TF_CT[seq(1, length(TF_CT), 2)]
  CR_R <- TF_CT[seq(2, length(TF_CT), 2)]
  TF_cen <- list()
  network <- list()
  gene_cen <- list()
  for (i in (CT)) {
    TF_CT <- unlist(strsplit(names(ct_regulon[CR_R == i]), split = "_"))
    TF_CT <- unique(TF_CT[seq(1, length(TF_CT), 2)])
    gene_CT <- unique((union(unlist(ct_regulon[CR_R == i]), TF_CT)))
    adj <- matrix(0, nrow = length(gene_CT), ncol = length(gene_CT))
    
    rownames(adj) <- gene_CT
    colnames(adj) <- gene_CT
    for (k in (1:length(ct_regulon[CR_R == i]))) {
      j <- ct_regulon[CR_R == i][k]
      adj1 <-
        matrix(0, nrow = length(gene_CT), ncol = length(gene_CT))
      rownames(adj1) <- gene_CT
      colnames(adj1) <- gene_CT
      RI_CT[paste(unlist(strsplit(names(j), split = "_"))[1], unlist(unname(j)), sep =
                    "_"), i]
      adj1[unlist(strsplit(names(j), split = "_"))[1], unlist(unname(j))] <-
        RI_CT[paste(unlist(strsplit(names(j), split = "_"))[1], unlist(unname(j)), sep =
                      "_"), i]
      adj1[unlist(unname(j)), unlist(strsplit(names(j), split = "_"))[1]] <-
        RI_CT[paste(unlist(strsplit(names(j), split = "_"))[1], unlist(unname(j)), sep =
                      "_"), i]
      adj = adj1 + adj
    }
    network[[i]] <- adj
    g <- graph.adjacency(adj, mode = "undirected", weighted = T)
    cen <- igraph::evcent(g)$vector
    TF_cen[[i]] <- cen[TF_CT]
    gene_cen[[i]] <- cen[setdiff(names(cen), TF_CT)]
    
  }
  m <- list()
  m[[1]] <- TF_cen
  m[[2]] <- gene_cen
  m[[3]] <- network
  return (m)
}



## ----------------------------------------------------------------------------------------------------------------
#find master TF and master gene, build a cell type gene regulatory network
#input:
# 1 - RAS_C1: a matrix with regulon-CT * cell. Regulon active score in cell type level
# 2 - graph.out: a factor variable. The predict cell cluster which return from get_gene_module function
#output:
# 1 - DR: differnent regulons.

calDR <-
  function(RAS_C1,
           graph.out,
           FindALLMarkers = T,
           ident.1 = 1,
           ident.2 = c(2, 3, 4)) {
    pbmc <- CreateSeuratObject(counts = RAS_C1)
    Idents(pbmc) <- graph.out
    if (FindALLMarkers == T) {
      #pbmc <- NormalizeData(pbmc)
      #pbmc <- ScaleData(pbmc, features = rownames(pbmc))
      #Idents(pbmc)<-graph.out
      DR <-
        FindAllMarkers(pbmc, only.pos = T, logfc.threshold = 0.25)
      DR <- DR[DR$p_val < 0.05, ]
      m <-
        unlist(strsplit(DR$gene, "-"))[seq(2, length(unlist(strsplit(DR$gene, "-"))), 2)]
      m <-
        unlist(strsplit(m, "ct"))[seq(2, length(unlist(strsplit(m, "ct"))), 2)]
      DR <- DR[as.numeric(m) == DR$cluster, ]
    } else{
      DR <-
        FindMarkers(pbmc,
                    ident.1 = ident.1,
                    ident.2 = ident.2,
                    min.pct = 0.25)
    }
    
    return (DR)
  }





## ----------------------------------------------------------------------------------------------------------------
# tf list (cell_gene matrix)
cal_tfmatrixs <- function(graph.out0,
                          graph.out,
                          TG_cell,
                          ct_regulon) {
  cells <- names(graph.out0)
  colnames(TG_cell) <- names(graph.out)
  TGs <- rownames(TG_cell)
  tfsmatrix <- list()
  TF_CT <- unlist(strsplit(names(ct_regulon), split = "_"))
  tflist <- unique(TF_CT[seq(1, length(TF_CT), 2)])
  for (tf in tflist) {
    genes <- unique(TGs[grep(tf, TGs)])
    x <- TG_cell[genes, cells]
    colnames(x) <- cells
    tfsmatrix[[tf]] <- t(x)
  }
  regulons <- c()
  for (tf in names(tfsmatrix)) {
    m <- tfsmatrix[[tf]]
    if (all(m == 0) == F) {
      regulons <- c(regulons, tf)
    }
  }
  tfsmatrix <- tfsmatrix[names(tfsmatrix) %in% regulons]
  return(tfsmatrix)
}


# input: regulon matrix (CT*gene)
# output: CT cluster
cal_clust <- function(m) {
  a <-
    data.frame(m)[which(apply(m, 1, function(row)
      all(row == 0)) == F), ]
  if (nrow(a) > 2) {
    b <- dist(a)
    c <- hclust(b)
    clusts <- c()
    for (i in 2:(nrow(a) - 1)) {
      d <- cutree(c, k = i , h = NULL)
      sil1 <- summary(silhouette(d, b))$avg.width
      names(sil1) <- i
      clusts <- c(clusts, sil1)
    }
    k <- names(which(clusts == max(clusts)))
    hc <- cutree(c, k , h = NULL)
    names(hc) <- rownames(a)
    if (length(which(apply(m, 1, function(row)
      all(row == 0)) == T)) > 0) {
      d <-
        c((max(hc) + 1):(max(hc) + length(which(
          apply(m, 1, function(row)
            all(row == 0)) == T
        ))))
      names(d) <-
        names(which(apply(m, 1, function(row)
          all(row == 0)) == T))#[-1]
      hc <- c(hc, d)
    }
    
  } else{
    hc <- c(1:nrow(m))
    names(hc) <- rownames(m)
  }
  
  return(hc)
}

# center
cal_meanmatrix <- function(m, mclust) {
  rownames(m) <- mclust
  c <- c()
  for (ct in unique(mclust)) {
    c <- rbind(c, apply(subset(m, rownames(m) == ct), 2, mean))
  }
  return(c)
}

Score_matrix <- function(meanm) {
  n <- length(which(apply(meanm, 1, function(row)
    all(row == 0)) == F))
  if (n == 1) {
    mScore <- sum(dist(meanm)) / (0.5 * (log(ncol(meanm))))
  } else{
    mScore <- sum(dist(meanm)) / (choose(n, 2) * (log(ncol(meanm))))
  }
  return(mScore)
  
}


cal_sp <- function(regulon) {
  library(cluster)
  m0 <- tfsmatrix[[regulon]]
  m <- cal_meanmatrix(m0, graph.out0)
  mclust <- cal_clust(m)
  mScore <- Score_matrix(cal_meanmatrix(m, mclust))
  sams <- c()
  for (j in 1:100) {
    m0clust <- sample(graph.out0, length(graph.out0))
    sams <-
      cbind(sams, Score_matrix(cal_meanmatrix(cal_meanmatrix(m0, m0clust), mclust)))
  }
  pvalue <- length(which(mScore < sams)) / length(sams)
  sp <- c(mScore, pvalue, mclust)
  names(sp) <- c('VS', 'pvalue', 'mclust')
  #rownames(sp)<-regulon
  print(pvalue)
  print(mScore)
  return(sp)
}
