inputdata.10x <- Read10X_h5("/home/wan268/hgt/newdata/GSM5265318_SJTALL005006_filtered_feature_bc_matrix.h5")
RNA <- inputdata.10x$`Gene Expression`
ATAC<- inputdata.10x$Peaks
pbmc <- CreateSeuratObject(counts = pbmc.rna@assays$RNA@counts)
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 1)
genelist<-c('CD3E','CD4','CD8A')
DotPlot(pbmc[pbmc@assays$RNA@var.features,], assay='RNA', features = genelist, dot.scale = 8)