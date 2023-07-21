# 07-06-2022

library(Seurat)
library(ArchR)
library(ggplot2)
library(parallel)

addArchRThreads(threads = 10)
addArchRGenome("hg38")

setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/")

############################
### 1.0 load arrow files ###
############################

ArrowFiles <- c("Sample_0.arrow", "Sample_1.arrow", "Sample_2.arrow", "Pool_2.arrow",
        "Pool_3.arrow", "Pool_4.arrow", "Pool_5.arrow", "Pool_6.arrow",
        "Pool_7.arrow", "Pool_8.arrow", "Pool_9.arrow", "Pool_10.arrow", "Pool_11.arrow")

##############################
### 2.0 Make ArchR Project ###
##############################

projIBD <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "projIBD",
  copyArrows = TRUE
)

saveArchRProject(ArchRProj = projIBD, outputDirectory = "Unfiltered-projIBD-07-2022", load = FALSE)
