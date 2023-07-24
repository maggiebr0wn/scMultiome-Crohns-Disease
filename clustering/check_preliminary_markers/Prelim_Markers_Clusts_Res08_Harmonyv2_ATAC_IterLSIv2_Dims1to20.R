#!/usr/bin/env Rscript

#library(ArchR)
library(ggplot2)
library(parallel)
library(scales)
library(dplyr)
library(Seurat)
library(stringr)


# 07-26-2022
# Load ArchR object, FindMarkers for top 8 cluster sets

#addArchRThreads(threads = 6)
#addArchRGenome("hg38")

### Load ArchR Object ###
#setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022")
#projIBD <- loadArchRProject("Filtered-Final-projIBD-ATACandGEX-07-2022")

### Load Seurat Object ###
seurat_IBD <- readRDS("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/Filtered-Final-projIBD-ATACandGEX-07-2022/IBD_seurat.RDS")

# cluster_list <- c("Clusts_Res1_ATAC_IterLSIv3_Dims1to20", 
#                   "Clusts_Res05_Combo_IterLSIv3_Dims1to20", 
#                   "Clusts_Res05_Harmonyv1_ATAC_IterLSIv3_Dims1to30", 
#                   "Clusts_Res08_Harmonyv2_ATAC_IterLSIv2_Dims1to20", 
#                   "Clusts_Res05_Harmonyv3_ATAC_IterLSIv2_Dims1to30", 
#                   "Clusts_Res05_Harmonyv2_Combo_IterLSIv3_Dims1to20", 
#                   "Clusts_Res05_Harmonyv2_Combo_IterLSIv3_Dims1to20", 
#                   "Clusts_Res05_Harmonyv3_Combo_IterLSIv3_Dims1to20")

cluster <- "Clusts_Res08_Harmonyv2_ATAC_IterLSIv2_Dims1to20"

gene_markers <- c(
  "MS4A1", "PAX5", "IRF8", "CD27", "CD24", "TBX21", "ZEB2", "ITGAX", "CXCR5", "TCL1A", # Naive/Mem/DN B cells
  "JCHAIN", "PRDM1", "IRF4", "XBP1", "CD19", "SDC1", "CD38", # ASC B cells
  "MYBL2", "MKI67", # proliferating cells
  "CD3D", # T cells
  "GZMH", "GZMB", # cytotoxic T cells
  "CD8A", # CD8+ T cells
  "IL7R", "CCR7", "S100A4", # Naive CD4+ T Cell, Memory CD4+
  "VCAN", # Monocytes
  "CD14", "LYZ", # CD14+ Mono
  "FCGR3A", "MS4A7", # FCGR3A+ Mono
  "GNLY", "NKG7", # NK
  "FCER1A", "CST3", #DC
  "PPBP" # Platelet
)

setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/Filtered-Final-projIBD-ATACandGEX-07-2022/CLUSTS_UMAPS/Prelim_Markers_WRS")

### GEX markers ###
DefaultAssay(seurat_IBD) <- "RNA"

seurat_IBD <- SetIdent(seurat_IBD, value = cluster)

markersGEX <- FindAllMarkers(
  seurat_IBD, assay = "RNA"
)

# keep output
filename <- paste(cluster, "_GEX_Markers.csv", sep = "")
write.csv(markersGEX, file = filename)

markersGEX %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

levels(seurat_IBD) <- eval(parse(text=(paste("str_sort(unique(seurat_IBD$", cluster, "), numeric = TRUE)", sep = ""))))

pdf(file = paste(cluster, "_GEX_Markers_dotplot", ".pdf", sep = ""), height = 6, width = 22)
DotPlot(seurat_IBD, features = rownames(seurat_IBD[rownames(seurat_IBD) %in% top5$gene,]), 
        assay = "RNA") + RotatedAxis() + scale_colour_gradient2(low = "grey46", mid = "dodgerblue3", high = "#FC4E07")
DotPlot(seurat_IBD, features = rownames(seurat_IBD[rownames(seurat_IBD) %in% top5$gene,]), 
        assay = "ATAC") + RotatedAxis() + scale_colour_gradient2(low = "grey46", mid = "darkmagenta", high = "gold")
dev.off()


### ATAC Markers ###
DefaultAssay(seurat_IBD) <- "ATAC"

seurat_IBD <- SetIdent(seurat_IBD, value = cluster)

markersATAC <- FindAllMarkers(
  seurat_IBD, assay = "ATAC"
)

# keep output
filename <- paste(cluster, "_ATAC_Markers.csv", sep = "")
write.csv(markersATAC, file = filename)

markersATAC %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC -> top5

levels(seurat_IBD) <- eval(parse(text=(paste("str_sort(unique(seurat_IBD$", cluster, "), numeric = TRUE)", sep = ""))))

pdf(file = paste(cluster, "_ATAC_Markers_dotplot", ".pdf", sep = ""), height = 6, width = 22)
DotPlot(seurat_IBD, features = rownames(seurat_IBD[rownames(seurat_IBD) %in% top5$gene,]), 
        assay = "RNA") + RotatedAxis() + scale_colour_gradient2(low = "grey46", mid = "dodgerblue3", high = "#FC4E07")
DotPlot(seurat_IBD, features = rownames(seurat_IBD[rownames(seurat_IBD) %in% top5$gene,]), 
        assay = "ATAC") + RotatedAxis() + scale_colour_gradient2(low = "grey46", mid = "darkmagenta", high = "gold")
dev.off()
