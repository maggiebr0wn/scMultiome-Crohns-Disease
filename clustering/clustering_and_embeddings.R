# 07-12-2022

library(Seurat)
library(ArchR)
library(ggplot2)
library(parallel)

# Clustering and Embeddings
# Cluster on various dimensionality reductions; loop through

addArchRThreads(threads = 6)
addArchRGenome("hg38")

getwd()

setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022")
projIBD <- loadArchRProject("Filtered-Final-projIBD-ATACandGEX-07-2022")
# numberOfCells(1): 81280
# medianTSS(1): 17.689
# medianFrags(1): 8369.5

#table(projIBD$Sample)
# Pool_10  Pool_11   Pool_2   Pool_3   Pool_4   Pool_5   Pool_6   Pool_7
#    6178     4669     5636     5575     6493     5489     5835     6310
#  Pool_8   Pool_9 Sample_0 Sample_1 Sample_2
#    9845     5560     4589     4104     7588


### Cluster for each Dimensionality Reduction ###
# Resolutions: 0.2, 0.5, 0.8, 1.0

# get all DimReds:
dimred_list <- names(projIBD@reducedDims)
# remove the RNA only ones
rna <- c("RNA_IterLSIv1_Dims1to15", "RNA_IterLSIv1_Dims1to20", "RNA_IterLSIv1_Dims1to30",
         "RNA_IterLSIv2_Dims1to15", "RNA_IterLSIv2_Dims1to20", "RNA_IterLSIv2_Dims1to30",
         "RNA_IterLSIv3_Dims1to15", "RNA_IterLSIv3_Dims1to20", "RNA_IterLSIv3_Dims1to30"
         )

final_dimred_list <- dimred_list[!dimred_list %in% rna]

resolution_list <- c(0.2, 0.5, 0.8, 1.0)

outpath <- "/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/Filtered-Final-projIBD-ATACandGEX-07-2022/CLUSTS_UMAPS/"

for (dimred in final_dimred_list){
  print(dimred)
  # Cluster and generate embeddings for resolutions
  for (res in resolution_list){
    print(res)
    print("Clustering! :)")
    # add Clusters
    char_res <- as.character(res)
    char_res <- gsub("\\.", "", char_res)
    clust_name <- paste("Clusts_Res", char_res, "_", dimred, sep = "")
    projIBD <- addClusters(
      input = projIBD,
      reducedDims = dimred,
      method = "Seurat",
      nam = clust_name,
      resolution = res,
      force = TRUE
    )
    write.csv(eval(parse(text=(paste("projIBD$",clust_name)))), file = paste(outpath, clust_name, ".csv", sep = ""))
    # add UMAP
    print("Adding UMAP! :)")
    umap_name <- paste(clust_name, "_UMAP", sep = "")
    projIBD <- addUMAP(
      ArchRProj = projIBD,
      reducedDims = dimred,
      name = umap_name,
      nNeighbors = 30,
      minDist = 0.5,
      metric = "cosine",
      force = TRUE
    )
    # add Embedding
    embedding <- getEmbedding(
      ArchRProj = projIBD,
      embedding = umap_name,
      returnDF = TRUE)
    colnames(embedding) <- c("UMAP_1", "UMAP_2")
    write.csv(embedding, file = paste(outpath, umap_name, ".csv", sep = ""))
  }
  saveArchRProject(projIBD)
}
