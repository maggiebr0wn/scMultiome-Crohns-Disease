#!/usr/bin/env Rscript

.libPaths("/storage/home/mfisher42/bin/R-4.2.0/lib64/R/library")

library(ArchR)
library(ggplot2)
library(parallel)
library(scales)
library(dplyr)
library(Seurat)
library(stringr)

library(BSgenome.Hsapiens.UCSC.hg38)

### 11-18-2022 Add Pseudobulk Replicates, then do Peak Calling ###

addArchRThreads(threads = 4)
addArchRGenome("hg38")

setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022")
singletons_projIBD <- loadArchRProject("Singletons-Final-projIBD-ATACandGEX-11-2022")

# add pseudobulk reps
singletons_projIBD$Celltypes_Sample <- paste(singletons_projIBD$Celltypes, singletons_projIBD$My_IDs, sep = "_")

# took 2.5 hours with 4 threads
singletons_projIBD <- addGroupCoverages(
                                        ArchRProj = singletons_projIBD,
                                        groupBy = "Celltypes_Sample",
                                        useLabels = FALSE
                                        )


saveArchRProject(singletons_projIBD)

# call peaks with MACS2; 38 min
singletons_projIBD <- addReproduciblePeakSet(
                                             ArchRProj = singletons_projIBD,
                                             groupBy = "Celltypes_Sample",
                                             pathToMacs2 = "/storage/home/mfisher42/.conda/envs/R_env/bin/macs2"
                                             )

saveArchRProject(singletons_projIBD)

# Add peak matrix
singletons_projIBD <- addPeakMatrix(singletons_projIBD)
# check:
getAvailableMatrices(singletons_projIBD)

saveArchRProject(singletons_projIBD)
