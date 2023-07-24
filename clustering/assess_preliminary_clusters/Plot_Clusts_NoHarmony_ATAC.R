# 07-14-2022

# ATAC Only, No Harmony
# Generate plots for each iteration:
# 1.) UMAPS by cluster, batch, pool
# 2.) Stacked barplots:
#   - clusters per batch, batches per cluster
#   - clusters per pool, pools per cluster
# 3.) Gene markers:
#   - On UMAPS
#   - On DotPlots

library(ArchR)
library(ggplot2)
library(parallel)
library(scales)
library(dplyr)
library(Seurat)
library(stringr)

addArchRThreads(threads = 4)
addArchRGenome("hg38")

### Load ArchR Object ###
setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022")
projIBD <- loadArchRProject("Filtered-Final-projIBD-ATACandGEX-07-2022")
# numberOfCells(1): 77871

### Get Cluster names and UMAP labels for this set ###
cluster_files <- list.files(path = "/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/Filtered-Final-projIBD-ATACandGEX-07-2022/CLUSTS_UMAPS/NoHarmony_ATAC",
                    pattern = ".csv")

# only keep cluster_files
cluster_files <- str_subset(cluster_files, "Clusts", negate = FALSE)
clusters <- gsub('.{4}$', '', cluster_files)

umap_files <- list.files(path = "/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/Filtered-Final-projIBD-ATACandGEX-07-2022/CLUSTS_UMAPS/UMAP_Coords/NoHarmony_ATAC",
                            pattern = ".csv")
umaps <- gsub('.{4}$', '', umap_files)

names(umaps) <- clusters

### Loop through each UMAP/Cluster set for plots: ###
setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/Filtered-Final-projIBD-ATACandGEX-07-2022/CLUSTS_UMAPS/NoHarmony_ATAC/Plots")

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

for (umap_set in umaps){
  print(umap_set)
  ### 1.0 UMAPs: groupby cluster, Sample, Batch ###
  print("Making UMAPS")
  # color by cluster
  p1 <- plotEmbedding(
    ArchRProj = projIBD,
    colorBy = "cellColData",
    name = names(umaps)[umaps == umap_set],
    embedding = umap_set
  ) + theme(legend.text=element_text(size=8))
  
  # plot by sample
  p2 <- plotEmbedding(
    ArchRProj = projIBD,
    colorBy = "cellColData",
    name = "Sample",
    embedding = umap_set
  ) + theme(legend.text=element_text(size=8))
  
  # color by batch
  p3 <- plotEmbedding(
    ArchRProj = projIBD,
    colorBy = "cellColData",
    name = "Batch",
    embedding = umap_set
  ) + theme(legend.text=element_text(size=8))
  
  # save UMAP plots:
  filename <- paste(umap_set, ".pdf", sep = "")
  pdf(file = filename)
  print(p1)
  print(p2)
  print(p3)
  dev.off()
  
  
  ### 2.0 Stacked Barplots: assess batch/cluster, sample/cluster proportions ###
  print("Making barplots: Sample per Cluster")
  # sample per cluster
  temp <- getCellColData(projIBD, select = names(umaps)[umaps == umap_set])
  colnames(temp) <- "x"
  barplot_data <- data.frame(projIBD$Sample, temp$x, 1)
  colnames(barplot_data) <- c("Sample", "Cluster", "X1")
  # Stacked + percent
  p1 <- ggplot(barplot_data, aes(fill=Sample, y=X1, x=Cluster)) +
    ggtitle("Proportion of Samples per Cluster") +
    xlab("Cluster Label") +
    ylab("Proportion of Cells") +
    geom_bar(position="fill", stat="identity") +
    theme(text = element_text(size=14)) +
    scale_x_discrete(limits = str_sort(unique(barplot_data$Cluster), numeric = TRUE)) +
    RotatedAxis()
  p1 <- p1 + guides(fill=guide_legend(title="Sample"))
  
  print("Making barplots: Batch per Cluster")
  # batch per cluster
  temp <- getCellColData(projIBD, select = names(umaps)[umaps == umap_set])
  colnames(temp) <- "x"
  barplot_data <- data.frame(projIBD$Batch, temp$x, 1)
  colnames(barplot_data) <- c("Batch", "Cluster", "X1")
  # Stacked + percent
  p2 <- ggplot(barplot_data, aes(fill=Batch, y=X1, x=Cluster)) +
    ggtitle("Proportion of Batches per Cluster") +
    xlab("Cluster Label") +
    ylab("Proportion of Cells") +
    geom_bar(position="fill", stat="identity") +
    theme(text = element_text(size=14)) +
    scale_x_discrete(limits = str_sort(unique(barplot_data$Cluster), numeric = TRUE)) +
    RotatedAxis()
  p2 <- p2 + guides(fill=guide_legend(title="Batch"))
  
  
  filename <- paste(names(umaps)[umaps == umap_set], "_per_cluster_barplots.pdf", sep = "")
  pdf(filename)
  print(p1)
  print(p2)
  dev.off()
  
  print("Making barplots: Cluster per Sample")
  # cluster per sample
  temp <- getCellColData(projIBD, select = names(umaps)[umaps == umap_set])
  colnames(temp) <- "x"
  barplot_data <- data.frame(projIBD$Sample, temp$x, 1)
  colnames(barplot_data) <- c("Sample", "Cluster", "X1")
  # Stacked + percent
  p1 <- ggplot(barplot_data, aes(fill=Cluster, y=X1, x=Sample)) +
    ggtitle("Proportion of Clusters per Sample") +
    xlab("Sample Label") +
    ylab("Proportion of Cells") +
    geom_bar(position="fill", stat="identity") +
    theme(text = element_text(size=14)) +
    scale_x_discrete(limits = str_sort(unique(barplot_data$Sample), numeric = TRUE)) +
    RotatedAxis()
  p1 <- p1 + guides(fill=guide_legend(title="Cluster"))
  
  print("Making barplots: Cluster per Batch")
  # cluster per batch
  temp <- getCellColData(projIBD, select = names(umaps)[umaps == umap_set])
  colnames(temp) <- "x"
  barplot_data <- data.frame(projIBD$Batch, temp$x, 1)
  colnames(barplot_data) <- c("Batch", "Cluster", "X1")
  
  # Stacked + percent
  p2 <- ggplot(barplot_data, aes(fill=Cluster, y=X1, x=Batch)) +
    ggtitle("Proportion of Clusters per Batch") +
    xlab("Batch Label") +
    ylab("Proportion of Cells") +
    geom_bar(position="fill", stat="identity") +
    theme(text = element_text(size=14)) +
    scale_x_discrete(limits = str_sort(unique(barplot_data$Batch), numeric = TRUE)) +
    RotatedAxis()
  p2 <- p2 + guides(fill=guide_legend(title="Cluster"))
  
  filename <- paste(names(umaps)[umaps == umap_set], "_per_sample_batch_barplots.pdf", sep = "")
  pdf(filename)
  print(p1)
  print(p2)
  dev.off()
  
  ### 3.0 Gene markers: ATAC/GEX UMAPs, ATAC/GEX Dotplots (create Seurat Object?) ###
  print("Making Marker Gene UMAPS")
  #addArchRThreads(threads = 1)

  lsi <- paste(
	       str_split_fixed(umap_set, "_", 6)[,3],
	       str_split_fixed(umap_set, "_", 6)[,4],
	       str_split_fixed(umap_set, "_", 6)[,5],
	       sep = "_"
	       )

  projIBD <- addImputeWeights(projIBD, reducedDims = lsi)
  #sub_projIBD <-  addImputeWeights(sub_projIBD, reducedDims = names(lsi_umaps)[lsi_umaps == umap_set])
  # umap: ATAC
  filename <- paste(names(umaps)[umaps == umap_set], "_ATAC_UMAP_GENE_MARKERS.pdf", sep = "")
  pdf(file = filename)
  for (gene in gene_markers){
    p <- plotEmbedding(
      ArchRProj = projIBD,
      colorBy = "GeneScoreMatrix",
      name = gene,
      embedding = umap_set,
      imputeWeights = getImputeWeights(projIBD),
      plotAs = "points"
    )
    plot(p)
  }
  dev.off()
  # umap: GEX
  filename <- paste(names(umaps)[umaps == umap_set], "_GEX_UMAP_GENE_MARKERS.pdf", sep = "")
  pdf(file = filename)
  for (gene in gene_markers){
    p <- plotEmbedding(
      ArchRProj = projIBD,
      colorBy = "GeneExpressionMatrix",
      name = gene,
      embedding = umap_set,
      imputeWeights = getImputeWeights(projIBD),
      plotAs = "points"
    )
    print(p)
    
  }
  dev.off()
  
  # trackplots: ATAC
  print("Making track plots")
  p <- plotBrowserTrack(
    ArchRProj = projIBD,
    groupBy = names(umaps)[umaps == umap_set],
    geneSymbol = gene_markers,
    upstream = 50000,
    downstream = 50000
  )
  
  pdf(file = paste(names(umaps)[umaps == umap_set], "_TRACK_Plots_GENE_MARKERS.pdf", sep = ""),
      height = 6, width = 4)
  for (gene in gene_markers){
    grid::grid.newpage()
    eval(parse(text=(paste("grid::grid.draw(p$",gene,")"))))
  }
  dev.off()
}
  
