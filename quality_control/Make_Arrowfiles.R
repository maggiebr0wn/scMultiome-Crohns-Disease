# 07-06-2022

library(ArchR)
library(Seurat)
library(parallel)

# 1.0 Generate ArrowFiles

addArchRThreads(threads = 10)
addArchRGenome("hg38")

setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/")

# get fragment files:

frag_files <-c(
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Sample_0/outs/atac_fragments.tsv.gz", # Sample_0
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Sample_1/outs/atac_fragments.tsv.gz", # Sample_1
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Sample_2/outs/atac_fragments.tsv.gz", # Sample_2
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_2/outs/atac_fragments.tsv.gz", # Pool_2
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_3/outs/atac_fragments.tsv.gz", # Pool_3
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_4/outs/atac_fragments.tsv.gz", # Pool_4
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_5/outs/atac_fragments.tsv.gz", # Pool_5
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_6/outs/atac_fragments.tsv.gz", # Pool_6
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_7/outs/atac_fragments.tsv.gz", # Pool_7
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_8/outs/atac_fragments.tsv.gz", # Pool_8
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_9/outs/atac_fragments.tsv.gz", # Pool_9
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_10/outs/atac_fragments.tsv.gz", # Pool_10
               "/storage/home/mfisher42/scProjects/CD_Subra/rawdata/Maggie_Organized_Data/cellranger_output/Pool_11/outs/atac_fragments.tsv.gz" # Pool_11
               )

names(frag_files) <- c("Sample_0", "Sample_1", "Sample_2", "Pool_2", "Pool_3", "Pool_4", "Pool_5", "Pool_6",
                       "Pool_7", "Pool_8", "Pool_9", "Pool_10", "Pool_11"
                       )

ArrowFiles <- createArrowFiles(
        inputFiles = frag_files,
        sampleNames = names(frag_files),
        minTSS = 0,
        minFrags = 1000,
        maxFrags = 1e+06,
        addTileMat = TRUE,
        addGeneScoreMat = TRUE
        )

ArrowFiles <-c("Sample_0.arrow", "Sample_1.arrow", "Sample_2.arrow", "Pool_2.arrow",
        "Pool_3.arrow", "Pool_4.arrow", "Pool_5.arrow", "Pool_6.arrow",
        "Pool_7.arrow", "Pool_8.arrow", "Pool_9.arrow", "Pool_10.arrow", "Pool_11.arrow")

