# GEO GSE292762 : Primary senescence and Secondary senescence scRNA-seq Analysis
🔗 [View dataset on GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292762)

R scripts enabling the main steps of the single-cell RNA-Seq analysis of primary senescence and secondary senescence data.


## QU vs IR: Primary Senescence scRNA-seq Integration
This repository contains Seurat R scripts for integrating and analyzing QU and IR samples from primary senescence single-cell RNA-seq experiments.

## QCMT vs SCMT: Secondary Senescence scRNA-seq Integration
This repository contains Seurat R scripts for integrating and analyzing QCMT and SCMT samples from secondary senescence single-cell RNA-seq experiments.

### Workflow Overview
- Load and QC filtering
- Gene filtering (expressed in ≥10 cells)
- Cell cycle scoring
- SCTransform + Harmony integration
- Dimensionality reduction (PCA, UMAP)
- Clustering at multiple resolutions
- Marker gene detection

### Requirements
- R (v4.3.3)
- Seurat (v4.4.0)
- harmony
- ggplot2

### 💾 Output Files
- Clustered Seurat object: `*_QI_harmony.RDS`, `*_QCSC_harmony.RDS`
- Marker gene list: `*_harmony_QU_vs_IR_genelist.txt`, `*_harmony_QCMT_vs_SCMT_genelist.txt`

### 📌 Notes
- Paths are relative and may need to be adjusted depending on your system.
- Please install packages before running scripts.
