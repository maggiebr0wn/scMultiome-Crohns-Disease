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

# 11-18-2022
# Add final labels after finishing
# - freemuxlet (all, filtered) --> base final labels on freemuxFilt
# - demuxlet (filter for matches with FreemuxFilt)
# Also add all metadata labels

addArchRThreads(threads = 4)
addArchRGenome("hg38")

### Load ArchR Object ###
setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022")
projIBD <- loadArchRProject("Filtered-Final-projIBD-ATACandGEX-07-2022")

# load metadata: from Seurat object with labels
metadata_labs <- read.csv("metadata_112022.csv", row.names = 1)
# sort:
sorted_labs <- metadata_labs[order(match(metadata_labs$Barcodes, getCellNames(projIBD))), ]

# Add metadata: Freemux_GEX_all
projIBD <- addCellColData(ArchRProj = projIBD,
                          data = sorted_labs$Freemux_GEX_all,
                          name = "Freemux_GEX_all",
                          cells = sorted_labs$Barcodes,
                          force = TRUE
                          )

# Add metadata: Freemux_GEX_filt
projIBD <- addCellColData(ArchRProj = projIBD,
                          data = sorted_labs$Freemux_GEX_filt,
                          name = "Freemux_GEX_filt",
                          cells = sorted_labs$Barcodes,
                          force = TRUE
                          )

# Add metadata: Demux_GEX_filt
projIBD <- addCellColData(ArchRProj = projIBD,
                          data = sorted_labs$Demux_GEX_filt,
                          name = "Demux_GEX_filt",
                          cells = sorted_labs$Barcodes,
                          force = TRUE
                          )

# Add metadata: My_IDs
projIBD <- addCellColData(ArchRProj = projIBD,
                          data = sorted_labs$My_IDs,
                          name = "My_IDs",
                          cells = sorted_labs$Barcodes,
                          force = TRUE
                          )

# Add metadata: PBMC_IDs
projIBD <- addCellColData(ArchRProj = projIBD,
                          data = sorted_labs$PBMC_IDs,
                          name = "PBMC_IDs",
                          cells = sorted_labs$Barcodes,
                          force = TRUE
                          )

# Add metadata: Genotype_IDs
projIBD <- addCellColData(ArchRProj = projIBD,
                          data = sorted_labs$Genotype_IDs,
                          name = "Genotype_IDs",
                          cells = sorted_labs$Barcodes,
                          force = TRUE
                          )

# Add metadata: Status
projIBD <- addCellColData(ArchRProj = projIBD,
                          data = sorted_labs$Status,
                          name = "Status",
                          cells = sorted_labs$Barcodes,
                          force = TRUE
                          )
# Add metadata: Sex
projIBD <- addCellColData(ArchRProj = projIBD,
                          data = sorted_labs$Sex,
                          name = "Sex",
                          cells = sorted_labs$Barcodes,
                          force = TRUE
                          )

# Add metadata: Status_Donor
projIBD <- addCellColData(ArchRProj = projIBD,
                          data = sorted_labs$Status_Donor,
                          name = "Status_Donor",
                          cells = sorted_labs$Barcodes,
                          force = TRUE
                          )

# Add metadata: Coarse_Celltypes
projIBD <- addCellColData(ArchRProj = projIBD,
                          data = sorted_labs$Coarse_Celltypes,
                          name = "Coarse_Celltypes",
                          cells = sorted_labs$Barcodes,
                          force = TRUE
                          )

# Add metadata: Celltypes
projIBD <- addCellColData(ArchRProj = projIBD,
                          data = sorted_labs$Celltypes,
                          name = "Celltypes",
                          cells = sorted_labs$Barcodes,
                          force = TRUE
                          )

# Subset ArchR Project
idxFiltered_Pass <- which(projIBD$cellNames %in% sorted_labs$Barcodes)
cellsPass <- projIBD$cellNames[idxFiltered_Pass]
length(cellsPass)

singletons_projIBD <- subsetArchRProject(
  ArchRProj = projIBD,
  cells = cellsPass,
  outputDirectory = "Singletons-Final-projIBD-ATACandGEX-11-2022",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

saveArchRProject(singletons_projIBD)






