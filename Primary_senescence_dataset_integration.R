# Script for integrating Primary senescence single-cell RNA-seq datasets
# Last updated: 2024-04-04
# Description: QU vs IR integration, batch correction via Harmony, clustering, and marker identification
# R version: 4.3.3
# Required packages: Seurat, harmony, ggplot2, SeuratObject

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
})

QU_path <- "../sc_cellranger/Output/QU/outs/filtered_feature_bc_matrix"
QU.data <- Read10X(data.dir = QU_path)
QU <- CreateSeuratObject(counts = QU.data, project = "QU", min.cells = 3, min.features = 200)
QU
QU <- RenameCells(QU, add.cell.id = "QU")

IR_path <- "../sc_cellranger/Output/IR/outs/filtered_feature_bc_matrix"
IR.data <- Read10X(data.dir = IR_path)
IR <- CreateSeuratObject(counts = IR.data, project = "IR", min.cells = 3, min.features = 200)
IR <- RenameCells(IR, add.cell.id = "IR")

Seurat <- merge(QU, y = IR)

Seurat[["percent.mt"]] <- PercentageFeatureSet(Seurat, pattern = "^MT-")
Seurat_V1 <- subset(x = Seurat,subset = (nFeature_RNA > 500)& (nCount_RNA >500) &(percent.mt <30))

counts <- Seurat_V1@assays$RNA@counts
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
Seurat_V2 <- CreateSeuratObject(filtered_counts, meta.data = Seurat_V1@meta.data)

s.genes <- Seurat::cc.genes$s.genes
g2m.genes <- Seurat::cc.genes$g2m.genes

# g2m, s gene expression-based scoring
Seurat_V2 <- CellCycleScoring(Seurat_V2,
                              g2m.features = g2m.genes,
                              s.features = s.genes,
                              set.ident = TRUE)

Seurat_V2 <- NormalizeData(Seurat_V2)
Seurat_V2 <- FindVariableFeatures(Seurat_V2)
Seurat_V2 <- ScaleData(Seurat_V2)
Seurat_V2 <- RunPCA(Seurat_V2,features = c(s.genes,g2m.genes))

Seurat_list <- SplitObject(Seurat_V2, split.by="orig.ident")
Seurat_list <- lapply(X = Seurat_list, 
                      FUN = SCTransform, 
                      method = "glmGamPoi", 
                      return.only.var.genes = FALSE,
                      vars.to.regress = c("nCount_RNA","percent.mt"))
head(Seurat_list[[1]]@meta.data)


var.features <- SelectIntegrationFeatures(object.list = Seurat_list, nfeatures = 3000)

Seurat_cominbed.1 <- merge(x = Seurat_list[[1]], y = Seurat_list[2:length(Seurat_list)], merge.data=TRUE)
VariableFeatures(Seurat_cominbed.1) <- var.features
Seurat_cominbed.1 <- ScaleData(Seurat_cominbed.1)
Seurat_cominbed.1 <- RunPCA(Seurat_cominbed.1, verbose = FALSE)
P37.SI.2<- Seurat_cominbed.1
pct <- P37.SI.2[["pca"]]@stdev/sum(P37.SI.2[["pca"]]@stdev)*100
cumu <- cumsum(pct)
co1 <- which(cumu >90 & pct <5)[1]
co2 <- sort(which((pct[1:length(pct)-1] -pct[2:length(pct)]) >0.1), decreasing = T)[1] +1
pcs <- min(co1,co2)
co1
co2
pcs
plot_df <- data.frame(pct = pct,
                      cumu = cumu,
                      rank = 1:length(pct))
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank >pcs)) +
  geom_text() +
  geom_vline(xintercept = 90, color = "gray") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "gray") +
  theme_bw()

PCA <- 30
P37.SI.2@reductions$pca.2 <- P37.SI.2@reductions$pca
P37.SI.2@reductions$pca.2@cell.embeddings <- P37.SI.2@reductions$pca@cell.embeddings[,1:PCA]
P37.SI.2@reductions$pca.2@feature.loadings <- P37.SI.2@reductions$pca@feature.loadings[,1:PCA]

P37.SI.2 <- RunHarmony(P37.SI.2, assay.use="SCT", group.by.vars = "orig.ident",reduction = "pca.2")
P37.SI.2 <- RunUMAP(P37.SI.2, reduction = "harmony", dims = 1:PCA)
P37.SI.2 <- FindNeighbors(P37.SI.2, reduction = "harmony", dims = 1:PCA) %>% FindClusters(resolution = c(0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8))

QI.combined_umap <- P37.SI.2
QI.combined_umap@meta.data$seurat_clusters <- QI.combined_umap@meta.data$SCT_snn_res.0.3
Idents(QI.combined_umap) <- 'seurat_clusters'

# Create a folder
FN <- "PCA30_mito30"
folder_path <- file.path("./Output/QU_vs_IR", FN)
dir.create(folder_path, showWarnings = T, recursive = TRUE)

saveRDS(QI.combined_umap, file = file.path(folder_path, paste0(FN, "_QI_harmony.RDS")))

QI_umap_CCS <- CellCycleScoring(QI.combined_umap, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

integ_all.markers <- FindAllMarkers(QI.combined_umap, min.pct = 0.25, logfc.threshold = 0, assay = "RNA")
write.table(integ_all.markers, file = file.path(folder_path, paste0(FN, "_harmony_QU_vs_IR_genelist.txt")), sep = "\t")