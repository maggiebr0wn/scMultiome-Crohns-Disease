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
# ID DARs for B Cells:
# Group2 versus all

addArchRThreads(threads = 6)
addArchRGenome("hg38")

### Load ArchR Object ###
setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022")
projIBD <- loadArchRProject("Singletons-Final-projIBD-ATACandGEX-11-2022")

# Subset  B Cells
idxSample <- BiocGenerics::which(projIBD$Coarse_Celltypes %in% "B_Cells")
cellsSample <- projIBD$cellNames[idxSample]
bcells_projIBD <- subsetCells(ArchRProj = projIBD, cellNames = cellsSample)

# only keep cells with more than 25 for pseudobulk
bcells_projIBD$Pseudobulk <- paste(bcells_projIBD$Celltypes, bcells_projIBD$Status_Donor, sep = "_")

pbs <- unique(bcells_projIBD$Pseudobulk)
keep <- NA
for (pb in pbs){
  if (length(which(bcells_projIBD$Pseudobulk == pb)) >= 25){
    print(TRUE)
    keep <- c(keep, pb)
  } else {
    print(FALSE)
  }
}

keep <- keep[complete.cases(keep)]
idxSample <- BiocGenerics::which(bcells_projIBD$Pseudobulk %in% keep)
cellsSample <- bcells_projIBD$cellNames[idxSample]
bcells_projIBD <- subsetCells(ArchRProj = bcells_projIBD, cellNames = cellsSample)

# Remove Healthy_6, Crohns_5 and unknowns
keep <- setdiff(unique(bcells_projIBD$Status_Donor), c("Healthy_6", "Crohns_5", "PossibleIBD_1", "Unknown_1"))
idxSample <- BiocGenerics::which(bcells_projIBD$Status_Donor %in% keep)
cellsSample <- bcells_projIBD$cellNames[idxSample]
bcells_projIBD <- subsetCells(ArchRProj = bcells_projIBD, cellNames = cellsSample)

ident1 <- c("Crohns_4", "Crohns_12")

ident2 <- setdiff(unique(bcells_projIBD$Status_Donor), ident1)

# Add Group Labels
cM <- confusionMatrix(bcells_projIBD$Status_Donor, bcells_projIBD$Status_Donor)
labelOld <- rownames(cM)

remapClust <- c(
        "Healthy_1" = "Group_4", "Crohns_1" = "Group_4",  "Crohns_2" = "Group_2", "Crohns_3" = "Group_3", "Healthy_2" = "Group_3", "Crohns_4" = "Group_1",
        "Healthy_4" = "Group_2", "Healthy_5" = "Group_2", "Healthy_3" = "Group_2", "Crohns_6" = "Group_4", "Crohns_9" = "Group_4", "Crohns_8" = "Group_2",
        "Crohns_7" = "Group_2", "Crohns_10" = "Group_3", "Crohns_11" = "Group_4", "Crohns_13" = "Group_3", "Crohns_12" = "Group_1", "Healthy_7" = "Group_4",
        "Crohns_14" = "Group_4", "Crohns_15" = "Group_3", "Healthy_8" = "Group_3", "Crohns_16" = "Group_4", "Healthy_9" = "Group_4", "Crohns_17" = "Group_3",
        "Crohns_18" = "Group_4", "Crohns_19" = "Group_3"
)

remapClust <- remapClust[names(remapClust) %in% labelOld]

bcells_projIBD$Groups <- mapLabels(bcells_projIBD$Status_Donor, newLabels = remapClust, oldLabels = labelOld)

markerTest <- getMarkerFeatures(
        ArchRProj = bcells_projIBD,
        useMatrix = "PeakMatrix",
        groupBy = "Groups",
        testMethod = "binomial",
        binarize = TRUE,
        bias = c("TSSEnrichment", "log10(nFrags)"),
        useGroups = "Group_2",
        bgdGroups = c("Group_1", "Group_3", "Group_4")
)

sigMarkers <- getMarkers(markerTest, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 0.25") # get signif

filename <- "/storage/home/mfisher42/scProjects/CD_Subra/Subset_Specific_DEGs_DARs/All_DARs_04112023/B_Cells/Group2_DARs.csv"
write.csv(sigMarkers$Group_2, file = filename)
