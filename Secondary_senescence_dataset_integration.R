# Script for integrating Secondary senescence single-cell RNA-seq datasets
# Last updated: 2025-04-04
# Description: QCMT vs SCMT integration, batch correction via Harmony, clustering, and marker identification
# R version: 4.3.3
# Required packages: Seurat, harmony, ggplot2, SeuratObject

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
})

QCMT_path <- "../sc_cellranger/Output/QCMT/outs/filtered_feature_bc_matrix"
QCMT.data <- Read10X(data.dir = QCMT_path)
QCMT <- CreateSeuratObject(counts = QCMT.data, project = "QCMT", min.cells = 3, min.features = 200)
QCMT
QCMT <- RenameCells(QCMT, add.cell.id = "QCMT")

SCMT_path <- "../sc_cellranger/Output/SCMT/SCMT_RE1_RE2_output/outs/filtered_feature_bc_matrix"
SCMT.data <- Read10X(data.dir = SCMT_path)
SCMT <- CreateSeuratObject(counts = SCMT.data, project = "SCMT", min.cells = 3, min.features = 200)
SCMT <- RenameCells(SCMT, add.cell.id = "SCMT")

Seurat <- merge(QCMT, y = SCMT)

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
P37.SI.2 <- FindNeighbors(P37.SI.2, reduction = "harmony", dims = 1:PCA) %>% FindClusters(resolution = c(0.2,0.25,0.26,0.27,0.28,0.29,0.3,0.35,0.4,0.5,0.6,0.7,0.8))

DimPlot(P37.SI.2, group.by = c("orig.ident"),pt.size = 1)

QCSC.combined_umap <- P37.SI.2

QCSC.combined_umap@meta.data$seurat_clusters <- QCSC.combined_umap@meta.data$SCT_snn_res.0.3
Idents(QCSC.combined_umap) <- 'seurat_clusters'

QCSC_umap_CCS <- CellCycleScoring(QCSC.combined_umap, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Create a folder 
FN <- "PCA30_mito30"
folder_path <- file.path("./Output/QCMT_vs_SCMT", FN)
dir.create(folder_path, showWarnings = T, recursive = TRUE)

saveRDS(QCSC.combined_umap, file = file.path(folder_path, paste0(FN, "_QCSC_harmony.RDS")))

QCSC_umap_CCS <- CellCycleScoring(QCSC.combined_umap, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

integ_all.markers <- FindAllMarkers(QCSC.combined_umap, min.pct = 0.25, logfc.threshold = 0, assay = "RNA")
write.table(integ_all.markers, file = file.path(folder_path, paste0(FN, "_harmony_QCMT_vs_SCMT_genelist.txt")), sep = "\t")