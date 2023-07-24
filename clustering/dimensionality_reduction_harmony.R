# 07-07-2022

library(Seurat)
library(ArchR)
library(ggplot2)
library(parallel)
library(textworks)

# Try various methods/parameters of dimensionality reductions

addArchRThreads(threads = 10)
addArchRGenome("hg38")

setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022")
projIBD <- loadArchRProject("Filtered-Final-projIBD-ATACandGEX-07-2022")
# numberOfCells(1): 81280
# medianTSS(1): 17.689
# medianFrags(1): 8369.5

### LSI and Harmony ###
# ATAC Only
# ATAC  RNA
# LSIv1, 2, 3, Dims 1:15, 1:20, 1:30
# W/ and W/o Harmony: corr cutoff 0.25, 0.5 0.75

# Add Batch
cM <- confusionMatrix(projIBD$Sample, projIBD$Sample)
labelOld <- rownames(cM)

remapClust <- c("Sample_0" = "Batch1",
                "Sample_1" = "Batch2", "Sample_2" = "Batch2", "Pool_2" = "Batch2", "Pool_3" = "Batch2",
                "Pool_4" = "Batch3", "Pool_5" = "Batch3", "Pool_6" = "Batch3", "Pool_7" = "Batch3",
                "Pool_8" = "Batch3", "Pool_9" = "Batch4", "Pool_10" = "Batch4", "Pool_11" = "Batch4")


remapClust <- remapClust[names(remapClust) %in% labelOld]
labelNew <- mapLabels(labelOld, oldLabels = names(remapClust), newLabels = remapClust)
projIBD$Batch <- mapLabels(projIBD$Sample, newLabels = labelNew, oldLabels = labelOld)


##################
### 1.0 LSI v1 ###
##################

### ATAC 1:15 ### --> interactive :)

projIBD <- addIterativeLSI(
     ArchRProj = projIBD,
     force = TRUE,
     useMatrix = "TileMatrix",
     name = "ATAC_IterLSIv1_Dims1to15",
     iterations = 3,
     clusterParams = list(
         resolution = 0.2,
         sampleCells = 10000,
         n.start = 10
       ),
     depthCol = "nFrags",
     varFeatures = 25000,
     dimsToUse = 1:15
   )

### RNA 1:15 ### --> interactive :)

projIBD <- addIterativeLSI(
     ArchRProj = projIBD,
     force = TRUE,
     useMatrix = "GeneExpressionMatrix",
     name = "RNA_IterLSIv1_Dims1to15",
     iterations = 3,
     clusterParams = list(
         resolution = 0.2,
         sampleCells = 10000,
         n.start = 10
       ),
     depthCol = "Gex_nUMI",
     varFeatures = 2500,
     firstSelection = "variable",
     binarize = FALSE,
     saveIterations = FALSE,
     dimsToUse = 1:15
   )

### Combo 1:15 ### --> interactive :)

projIBD <- addCombinedDims(
     projIBD,
     reducedDims = c("ATAC_IterLSIv1_Dims1to15", "RNA_IterLSIv1_Dims1to15"),
     name =  "Combo_IterLSIv1_Dims1to15")

### ATAC 1:20 ### --> interactive :)

projIBD <- addIterativeLSI(
     ArchRProj = projIBD,
     force = TRUE,
     useMatrix = "TileMatrix",
     name = "ATAC_IterLSIv1_Dims1to20",
     iterations = 3,
     clusterParams = list(
         resolution = 0.2,
         sampleCells = 10000,
         n.start = 10
       ),
     depthCol = "nFrags",
     varFeatures = 25000,
     dimsToUse = 1:20
   )

### RNA 1:20 ### --> interactive :)

projIBD <- addIterativeLSI(
     ArchRProj = projIBD,
     force = TRUE,
     useMatrix = "GeneExpressionMatrix",
     name = "RNA_IterLSIv1_Dims1to20",
     iterations = 3,
     clusterParams = list(
         resolution = 0.2,
         sampleCells = 10000,
         n.start = 10
       ),
     depthCol = "Gex_nUMI",
     varFeatures = 2500,
     firstSelection = "variable",
     binarize = FALSE,
     saveIterations = FALSE,
     dimsToUse = 1:20
   )

### Combo 1:20 ### --> interactive :)

projIBD <- addCombinedDims(
     projIBD,
     reducedDims = c("ATAC_IterLSIv1_Dims1to20", "RNA_IterLSIv1_Dims1to20"),
     name =  "Combo_IterLSIv1_Dims1to20")

### ATAC 1:30 ### --> interactive :)

projIBD <- addIterativeLSI(
     ArchRProj = projIBD,
     force = TRUE,
     useMatrix = "TileMatrix",
     name = "ATAC_IterLSIv1_Dims1to30",
     iterations = 3,
     clusterParams = list(
         resolution = 0.2,
         sampleCells = 10000,
         n.start = 10
       ),
     depthCol = "nFrags",
     varFeatures = 25000,
     dimsToUse = 1:30
   )

### RNA 1:30 ### --> interactive :)

projIBD <- addIterativeLSI(
     ArchRProj = projIBD,
     force = TRUE,
     useMatrix = "GeneExpressionMatrix",
     name = "RNA_IterLSIv1_Dims1to30",
     iterations = 3,
     clusterParams = list(
         resolution = 0.2,
         sampleCells = 10000,
         n.start = 10
       ),
     depthCol = "Gex_nUMI",
     varFeatures = 2500,
     firstSelection = "variable",
     binarize = FALSE,
     saveIterations = FALSE,
     dimsToUse = 1:30
   )

### Combo 1:30 ### --> interactive :)

projIBD <- addCombinedDims(
     projIBD,
     reducedDims = c("ATAC_IterLSIv1_Dims1to30", "RNA_IterLSIv1_Dims1to30"),
     name =  "Combo_IterLSIv1_Dims1to30")

saveArchRProject(projIBD)

##################
### 2.0 LSI v2 ###
##################

### ATAC 1:15 ### --> interactive :)

projIBD <- addIterativeLSI(
     ArchRProj = projIBD,
     force = TRUE,
     useMatrix = "TileMatrix",
     name = "ATAC_IterLSIv2_Dims1to15",
     iterations = 3,
     clusterParams = list(
         resolution = 0.8,
         sampleCells = 10000,
         n.start = 10
       ),
     depthCol = "nFrags",
     varFeatures = 25000,
     dimsToUse = 1:15
   )

### RNA 1:15 ### --> interactive :)

projIBD <- addIterativeLSI(
     ArchRProj = projIBD,
     force = TRUE,
     useMatrix = "GeneExpressionMatrix",
     name = "RNA_IterLSIv2_Dims1to15",
     iterations = 3,
     clusterParams = list(
         resolution = 0.8,
         sampleCells = 10000,
         n.start = 10
       ),
     depthCol = "Gex_nUMI",
     varFeatures = 2500,
     firstSelection = "variable",
     binarize = FALSE,
     saveIterations = FALSE,
     dimsToUse = 1:15
   )

### Combo 1:15 ### --> interactive :)

projIBD <- addCombinedDims(
     projIBD,
     reducedDims = c("ATAC_IterLSIv2_Dims1to15", "RNA_IterLSIv2_Dims1to15"),
     name =  "Combo_IterLSIv2_Dims1to15")

### ATAC 1:20 ### --> interactive :)

projIBD <- addIterativeLSI(
     ArchRProj = projIBD,
     force = TRUE,
     useMatrix = "TileMatrix",
     name = "ATAC_IterLSIv2_Dims1to20",
     iterations = 3,
     clusterParams = list(
         resolution = 0.8,
         sampleCells = 10000,
         n.start = 10
       ),
     depthCol = "nFrags",
     varFeatures = 25000,
     dimsToUse = 1:20
   )

### RNA 1:20 ### --> interactive :)

projIBD <- addIterativeLSI(
     ArchRProj = projIBD,
     force = TRUE,
     useMatrix = "GeneExpressionMatrix",
     name = "RNA_IterLSIv2_Dims1to20",
     iterations = 3,
     clusterParams = list(
         resolution = 0.8,
         sampleCells = 10000,
         n.start = 10
       ),
     depthCol = "Gex_nUMI",
     varFeatures = 2500,
     firstSelection = "variable",
     binarize = FALSE,
     saveIterations = FALSE,
     dimsToUse = 1:20
   )

### Combo 1:20 ### --> interactive :)

projIBD <- addCombinedDims(
     projIBD,
     reducedDims = c("ATAC_IterLSIv2_Dims1to20", "RNA_IterLSIv2_Dims1to20"),
     name =  "Combo_IterLSIv2_Dims1to20")

### ATAC 1:30 ### --> interactive :)

projIBD <- addIterativeLSI(
     ArchRProj = projIBD,
     force = TRUE,
     useMatrix = "TileMatrix",
     name = "ATAC_IterLSIv2_Dims1to30",
     iterations = 3,
     clusterParams = list(
         resolution = 0.8,
         sampleCells = 10000,
         n.start = 10
       ),
     depthCol = "nFrags",
     varFeatures = 25000,
     dimsToUse = 1:30
   )

### RNA 1:30 ### --> interactive :)

projIBD <- addIterativeLSI(
  ArchRProj = projIBD,
  force = TRUE,
  useMatrix = "GeneExpressionMatrix",
  name = "RNA_IterLSIv2_Dims1to30",
  iterations = 3,
  clusterParams = list(
    resolution = 0.8,
    sampleCells = 10000,
    n.start = 10
  ),
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  saveIterations = FALSE,
  dimsToUse = 1:30
)

### Combo 1:30 ###

projIBD <- addCombinedDims(
  projIBD,
  reducedDims = c("ATAC_IterLSIv2_Dims1to30", "RNA_IterLSIv2_Dims1to30"),
  name =  "Combo_IterLSIv2_Dims1to30")

saveArchRProject(projIBD)

##################
### 3.0 LSI v3 ###
##################

### ATAC 1:15 ### --> interactive :)

projIBD <- addIterativeLSI(
  ArchRProj = projIBD,
  force = TRUE,
  useMatrix = "TileMatrix",
  name = "ATAC_IterLSIv3_Dims1to15",
  iterations = 5,
  clusterParams = list(
    resolution = c(0.2, 0.4, 0.5, 0.8),
    sampleCells = 10000,
    n.start = 10
  ),
  depthCol = "nFrags",
  varFeatures = 25000,
  dimsToUse = 1:15
)

### RNA 1:15 ###

projIBD <- addIterativeLSI(
  ArchRProj = projIBD,
  force = TRUE,
  useMatrix = "GeneExpressionMatrix",
  name = "RNA_IterLSIv3_Dims1to15",
  iterations = 5,
  clusterParams = list(
    resolution = c(0.2, 0.4, 0.5, 0.8),
    sampleCells = 10000,
    n.start = 10
  ),
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  saveIterations = FALSE,
  dimsToUse = 1:15
)

### Combo 1:15 ###

projIBD <- addCombinedDims(
  projIBD,
  reducedDims = c("ATAC_IterLSIv3_Dims1to15", "RNA_IterLSIv3_Dims1to15"),
  name =  "Combo_IterLSIv3_Dims1to15")

### ATAC 1:20 ###

projIBD <- addIterativeLSI(
  ArchRProj = projIBD,
  force = TRUE,
  useMatrix = "TileMatrix",
  name = "ATAC_IterLSIv3_Dims1to20",
  iterations = 5,
  clusterParams = list(
    resolution = c(0.2, 0.4, 0.5, 0.8),
    sampleCells = 10000,
    n.start = 10
  ),
  depthCol = "nFrags",
  varFeatures = 25000,
  dimsToUse = 1:20
)

### RNA 1:20 ###

projIBD <- addIterativeLSI(
  ArchRProj = projIBD,
  force = TRUE,
  useMatrix = "GeneExpressionMatrix",
  name = "RNA_IterLSIv3_Dims1to20",
  iterations = 5,
  clusterParams = list(
    resolution = c(0.2, 0.4, 0.5, 0.8),
    sampleCells = 10000,
    n.start = 10
  ),
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  saveIterations = FALSE,
  dimsToUse = 1:20
)

### Combo 1:20 ###

projIBD <- addCombinedDims(
  projIBD,
  reducedDims = c("ATAC_IterLSIv3_Dims1to20", "RNA_IterLSIv3_Dims1to20"),
  name =  "Combo_IterLSIv3_Dims1to20")

### ATAC 1:30 ###

projIBD <- addIterativeLSI(
  ArchRProj = projIBD,
  force = TRUE,
  useMatrix = "TileMatrix",
  name = "ATAC_IterLSIv3_Dims1to30",
  iterations = 5,
  clusterParams = list(
    resolution = c(0.2, 0.4, 0.5, 0.8),
    sampleCells = 10000,
    n.start = 10
  ),
  depthCol = "nFrags",
  varFeatures = 25000,
  dimsToUse = 1:30
)

### RNA 1:30 ###

projIBD <- addIterativeLSI(
  ArchRProj = projIBD,
  force = TRUE,
  useMatrix = "GeneExpressionMatrix",
  name = "RNA_IterLSIv3_Dims1to30",
  iterations = 5,
  clusterParams = list(
    resolution = c(0.2, 0.4, 0.5, 0.8),
    sampleCells = 10000,
    n.start = 10
  ),
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  saveIterations = FALSE,
  dimsToUse = 1:30
)

### Combo 1:30 ###

projIBD <- addCombinedDims(
  projIBD,
  reducedDims = c("ATAC_IterLSIv3_Dims1to30", "RNA_IterLSIv3_Dims1to30"),
  name =  "Combo_IterLSIv3_Dims1to30")

saveArchRProject(projIBD)

######################
### 4.0 Harmony v1 ###
######################

####### LSI v1 ####### 

## ATAC Only 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv1_Dims1to15",
  name = "Harmonyv1_ATAC_IterLSIv1_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## Combo 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv1_Dims1to15",
  name = "Harmonyv1_Combo_IterLSIv1_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## ATAC Only 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv1_Dims1to20",
  name = "Harmonyv1_ATAC_IterLSIv1_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## Combo 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv1_Dims1to20",
  name = "Harmonyv1_Combo_IterLSIv1_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## ATAC Only 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv1_Dims1to30",
  name = "Harmonyv1_ATAC_IterLSIv1_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## Combo 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv1_Dims1to30",
  name = "Harmonyv1_Combo_IterLSIv1_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

saveArchRProject(projIBD)

####### LSI v2 ####### 

## ATAC Only 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv2_Dims1to15",
  name = "Harmonyv1_ATAC_IterLSIv2_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## Combo 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv2_Dims1to15",
  name = "Harmonyv1_Combo_IterLSIv2_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## ATAC Only 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv2_Dims1to20",
  name = "Harmonyv1_ATAC_IterLSIv2_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## Combo 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv2_Dims1to20",
  name = "Harmonyv1_Combo_IterLSIv2_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## ATAC Only 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv2_Dims1to30",
  name = "Harmonyv1_ATAC_IterLSIv2_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## Combo 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv2_Dims1to30",
  name = "Harmonyv1_Combo_IterLSIv2_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

saveArchRProject(projIBD)

####### LSI v3 ####### 

## ATAC Only 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv3_Dims1to15",
  name = "Harmonyv1_ATAC_IterLSIv3_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## Combo 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv3_Dims1to15",
  name = "Harmonyv1_Combo_IterLSIv3_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## ATAC Only 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv3_Dims1to20",
  name = "Harmonyv1_ATAC_IterLSIv3_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## Combo 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv3_Dims1to20",
  name = "Harmonyv1_Combo_IterLSIv3_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## ATAC Only 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv3_Dims1to30",
  name = "Harmonyv1_ATAC_IterLSIv3_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

## Combo 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv3_Dims1to30",
  name = "Harmonyv1_Combo_IterLSIv3_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.25,
  force = TRUE
)

saveArchRProject(projIBD)

######################
### 5.0 Harmony v2 ###
######################

####### LSI v1 ####### 

## ATAC Only 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv1_Dims1to15",
  name = "Harmonyv2_ATAC_IterLSIv1_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## Combo 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv1_Dims1to15",
  name = "Harmonyv2_Combo_IterLSIv1_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## ATAC Only 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv1_Dims1to20",
  name = "Harmonyv2_ATAC_IterLSIv1_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## Combo 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv1_Dims1to20",
  name = "Harmonyv2_Combo_IterLSIv1_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## ATAC Only 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv1_Dims1to30",
  name = "Harmonyv2_ATAC_IterLSIv1_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## Combo 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv1_Dims1to30",
  name = "Harmonyv2_Combo_IterLSIv1_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

saveArchRProject(projIBD)

####### LSI v2 ####### 

## ATAC Only 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv2_Dims1to15",
  name = "Harmonyv2_ATAC_IterLSIv2_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## Combo 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv2_Dims1to15",
  name = "Harmonyv2_Combo_IterLSIv2_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## ATAC Only 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv2_Dims1to20",
  name = "Harmonyv2_ATAC_IterLSIv2_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## Combo 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv2_Dims1to20",
  name = "Harmonyv2_Combo_IterLSIv2_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## ATAC Only 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv2_Dims1to30",
  name = "Harmonyv2_ATAC_IterLSIv2_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## Combo 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv2_Dims1to30",
  name = "Harmonyv2_Combo_IterLSIv2_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

saveArchRProject(projIBD)

####### LSI v3 ####### 

## ATAC Only 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv3_Dims1to15",
  name = "Harmonyv2_ATAC_IterLSIv3_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## Combo 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv3_Dims1to15",
  name = "Harmonyv2_Combo_IterLSIv3_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## ATAC Only 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv3_Dims1to20",
  name = "Harmonyv2_ATAC_IterLSIv3_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## Combo 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv3_Dims1to20",
  name = "Harmonyv2_Combo_IterLSIv3_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## ATAC Only 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv3_Dims1to30",
  name = "Harmonyv2_ATAC_IterLSIv3_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

## Combo 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv3_Dims1to30",
  name = "Harmonyv2_Combo_IterLSIv3_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.5,
  force = TRUE
)

saveArchRProject(projIBD)

######################
### 6.0 Harmony v3 ###
######################

####### LSI v1 ####### 

## ATAC Only 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv1_Dims1to15",
  name = "Harmonyv3_ATAC_IterLSIv1_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## Combo 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv1_Dims1to15",
  name = "Harmonyv3_Combo_IterLSIv1_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## ATAC Only 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv1_Dims1to20",
  name = "Harmonyv3_ATAC_IterLSIv1_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## Combo 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv1_Dims1to20",
  name = "Harmonyv3_Combo_IterLSIv1_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## ATAC Only 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv1_Dims1to30",
  name = "Harmonyv3_ATAC_IterLSIv1_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## Combo 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv1_Dims1to30",
  name = "Harmonyv3_Combo_IterLSIv1_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

saveArchRProject(projIBD)

####### LSI v2 ####### 

## ATAC Only 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv2_Dims1to15",
  name = "Harmonyv3_ATAC_IterLSIv2_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## Combo 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv2_Dims1to15",
  name = "Harmonyv3_Combo_IterLSIv2_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## ATAC Only 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv2_Dims1to20",
  name = "Harmonyv3_ATAC_IterLSIv2_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## Combo 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv2_Dims1to20",
  name = "Harmonyv3_Combo_IterLSIv2_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## ATAC Only 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv2_Dims1to30",
  name = "Harmonyv3_ATAC_IterLSIv2_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## Combo 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv2_Dims1to30",
  name = "Harmonyv3_Combo_IterLSIv2_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

saveArchRProject(projIBD)

####### LSI v3 ####### 

## ATAC Only 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv3_Dims1to15",
  name = "Harmonyv3_ATAC_IterLSIv3_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## Combo 1:15 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv3_Dims1to15",
  name = "Harmonyv3_Combo_IterLSIv3_Dims1to15",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## ATAC Only 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv3_Dims1to20",
  name = "Harmonyv3_ATAC_IterLSIv3_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## Combo 1:20 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv3_Dims1to20",
  name = "Harmonyv3_Combo_IterLSIv3_Dims1to20",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## ATAC Only 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "ATAC_IterLSIv3_Dims1to30",
  name = "Harmonyv3_ATAC_IterLSIv3_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

## Combo 1:30 ##

projIBD <- addHarmony(
  ArchRProj = projIBD,
  reducedDims = "Combo_IterLSIv3_Dims1to30",
  name = "Harmonyv3_Combo_IterLSIv3_Dims1to30",
  groupBy = "Batch",
  corCutOff = 0.75,
  force = TRUE
)

saveArchRProject(projIBD)
