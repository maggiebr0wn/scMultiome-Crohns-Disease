# 07-06-2022

library(Seurat)
library(ArchR)
library(ggplot2)
library(stringr)
library(parallel)
library(textworks)

addArchRThreads(threads = 10)
addArchRGenome("hg38")

# This script adds paired GEX data
# Filters out low quality GEX cells

##############################
### 1.0 Load ArchR Project ###
##############################

setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/")

projIBD <- loadArchRProject("Filtered-projIBD-ATAC-07-2022")
# numberOfCells(1): 97593
# medianTSS(1): 17.184
# medianFrags(1): 6882

###############################
### 2.0 Load scRNA-seq Data ###
###############################


setwd("/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/")
seRNA <- import10xFeatureMatrix(
  input = c(
    "Sample_0/outs/filtered_feature_bc_matrix.h5",
    "Sample_1/outs/filtered_feature_bc_matrix.h5",
    "Sample_2/outs/filtered_feature_bc_matrix.h5",
    "Pool_2/outs/filtered_feature_bc_matrix.h5",
    "Pool_3/outs/filtered_feature_bc_matrix.h5",
    "Pool_4/outs/filtered_feature_bc_matrix.h5",
    "Pool_5/outs/filtered_feature_bc_matrix.h5",
    "Pool_6/outs/filtered_feature_bc_matrix.h5",
    "Pool_7/outs/filtered_feature_bc_matrix.h5",
    "Pool_8/outs/filtered_feature_bc_matrix.h5",
    "Pool_9/outs/filtered_feature_bc_matrix.h5",
    "Pool_10/outs/filtered_feature_bc_matrix.h5",
    "Pool_11/outs/filtered_feature_bc_matrix.h5"),
  names = c("Sample_0", "Sample_1", "Sample_2", "Pool_2", "Pool_3", "Pool_4", "Pool_5", "Pool_6",
            "Pool_7", "Pool_8", "Pool_9", "Pool_10", "Pool_11")
)

# > dim(seRNA)
# [1]  36577 106296

##############################
### 3.0 Add scRNA-seq Data ###
##############################

# Subset scRNA cells to only include ones in scATAC
seRNA_alpha <- seRNA[, colnames(seRNA) %in% projIBD$cellNames]
# > dim(seRNA_alpha)
# [1] 36577 93907

projIBD <- addGeneExpressionMatrix(input = projIBD, seRNA = seRNA_alpha, force = TRUE)

# Subset scATAC cells to only include ones in scRNA
setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/")

filtered_projIBD <- subsetArchRProject(
  ArchRProj = projIBD,
  cells = colnames(seRNA_alpha),
  outputDirectory = "Filtered-Partially-projIBD-ATACandGEX-07-2022",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)
# numberOfCells(1): 93907
# medianTSS(1): 17.317
# medianFrags(1): 7337

Multiome_RNA_QC <- read.delim("Multiome_RNA_QC_commonPools_n31.txt")

# Update with separately QC'd GEX cells

cellIDS <- Multiome_RNA_QC$Cells
cellIDS <- str_replace(cellIDS, "M1_", "Sample_0#")
cellIDS <- str_replace(cellIDS, "M2_", "Sample_1#")
cellIDS <- str_replace(cellIDS, "M3_", "Sample_2#")
cellIDS <- str_replace(cellIDS, "M4_", "Pool_2#")
cellIDS <- str_replace(cellIDS, "M5_", "Pool_3#")
cellIDS <- str_replace(cellIDS, "M6_", "Pool_4#")
cellIDS <- str_replace(cellIDS, "M7_", "Pool_5#")
cellIDS <- str_replace(cellIDS, "M8_", "Pool_6#")
cellIDS <- str_replace(cellIDS, "M9_", "Pool_7#")
cellIDS <- str_replace(cellIDS, "M10_", "Pool_8#")
cellIDS <- str_replace(cellIDS, "M11_", "Pool_9#")
cellIDS <- str_replace(cellIDS, "M12_", "Pool_10#")
cellIDS <- str_replace(cellIDS, "M13_", "Pool_11#")

idxFiltered_Pass <- which(filtered_projIBD$cellNames %in% cellIDS)
cellsPass <- filtered_projIBD$cellNames[idxFiltered_Pass]
length(cellsPass)

final_filtered_projIBD <- subsetArchRProject(
  ArchRProj = filtered_projIBD,
  cells = cellsPass,
  outputDirectory = "Filtered-Final-projIBD-ATACandGEX-07-2022",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

