#!/usr/bin/env Rscript

.libPaths("/storage/home/mfisher42/bin/R-4.2.0/lib64/R/library")

library(ArchR)
library(ggplot2)
library(JASPAR2020)
library(parallel)
library(scales)
library(dplyr)
library(DirichletMultinomial)
library(motifmatchr)
library(Seurat)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)

# 04-11-2023
# ID DARs for T Cells:
# GroupC versus all

addArchRThreads(threads = 6)
addArchRGenome("hg38")

### Load ArchR Object ###
setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022")
projIBD <- loadArchRProject("Singletons-Final-projIBD-ATACandGEX-11-2022")

# Subset T Cells
idxSample <- BiocGenerics::which(projIBD$Coarse_Celltypes %in% "T_Cells")
cellsSample <- projIBD$cellNames[idxSample]
tcells_projIBD <- subsetCells(ArchRProj = projIBD, cellNames = cellsSample)

# only keep cells with more than 25 for pseudobulk
tcells_projIBD$Pseudobulk <- paste(tcells_projIBD$Celltypes, tcells_projIBD$Status_Donor, sep = "_")

pbs <- unique(tcells_projIBD$Pseudobulk)
keep <- NA
for (pb in pbs){
  if (length(which(tcells_projIBD$Pseudobulk == pb)) >= 25){
    print(TRUE)
    keep <- c(keep, pb)
  } else {
    print(FALSE)
  }
}

keep <- keep[complete.cases(keep)]
idxSample <- BiocGenerics::which(tcells_projIBD$Pseudobulk %in% keep)
cellsSample <- tcells_projIBD$cellNames[idxSample]
tcells_projIBD <- subsetCells(ArchRProj = tcells_projIBD, cellNames = cellsSample)

# Remove unknowns
keep <- setdiff(unique(tcells_projIBD$Status_Donor), c("PossibleIBD_1", "Unknown_1"))
idxSample <- BiocGenerics::which(tcells_projIBD$Status_Donor %in% keep)
cellsSample <- tcells_projIBD$cellNames[idxSample]
tcells_projIBD <- subsetCells(ArchRProj = tcells_projIBD, cellNames = cellsSample)

# Add Group Labels
cM <- confusionMatrix(tcells_projIBD$Status_Donor, tcells_projIBD$Status_Donor)
labelOld <- rownames(cM)

remapClust <- c(
        "Healthy_1" = "Group_D", "Crohns_1" = "Group_C",  "Crohns_2" = "Group_B", "Crohns_3" = "Group_D", "Crohns_4" = "Group_A", "Healthy_2" = "Group_D",
        "Healthy_5" = "Group_B", "Healthy_3" = "Group_B", "Healthy_4" = "Group_B", "Crohns_9" = "Group_C", "Crohns_6" = "Group_D", "Crohns_5" = "Group_D",
        "Crohns_8" = "Group_B", "Crohns_7" = "Group_B", "Crohns_10" = "Group_B", "Crohns_11" = "Group_C",  "Crohns_12" = "Group_A", "Crohns_13" = "Group_B",
        "Healthy_6" = "Group_B",  "Healthy_7" = "Group_C", "Crohns_14" = "Group_C", "Crohns_15" = "Group_D", "Healthy_8" = "Group_B", "Healthy_9" = "Group_C",
        "Crohns_16" = "Group_C",  "Crohns_17" = "Group_D", "Crohns_18" = "Group_C", "Crohns_19" = "Group_C"
)

remapClust <- remapClust[names(remapClust) %in% labelOld]

tcells_projIBD$Groups <- mapLabels(tcells_projIBD$Status_Donor, newLabels = remapClust, oldLabels = labelOld)

markerTest <- getMarkerFeatures(
        ArchRProj = tcells_projIBD,
        useMatrix = "PeakMatrix",
        groupBy = "Groups",
        testMethod = "binomial",
        binarize = TRUE,
        bias = c("TSSEnrichment", "log10(nFrags)"),
        useGroups = "Group_C",
        bgdGroups = c("Group_A", "Group_B", "Group_D")
)

sigMarkers <- getMarkers(markerTest, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.25") # get signif

filename <- "/storage/home/mfisher42/scProjects/CD_Subra/Subset_Specific_DEGs_DARs/All_DARs_04112023/T_Cells/GroupC_DARs.csv"
write.csv(sigMarkers$Group_C, file = filename)
