#!/usr/bin/env Rscript

library(ggplot2)
library(parallel)
library(scales)
library(dplyr)
library(pheatmap)
library(Seurat)
library(stringr)
library(RColorBrewer)

library(CellChat)
library(patchwork)

### 12-07-22 ###
# CellChat for CD4 and CD12 vs Healthy

# Load Seurat object
ibd_seurat <- readRDS("/Users/maggiebrown/scProject/CD_Multiome/Multiome_07_2022/seurat_obj_FINAL_LABELS_112022.rds")
ibd_seurat <- SetIdent(ibd_seurat, value = "Status_Donor")

crohns_ibd_seurat <- subset(ibd_seurat, idents = c("Crohns_4","Crohns_12"))
healthy_ibd_seurat <- subset(ibd_seurat, idents =  c( "Healthy_1", "Healthy_2", "Healthy_3",
                                                     "Healthy_4", "Healthy_5", "Healthy_6",
                                                     "Healthy_7", "Healthy_8", "Healthy_9"))

setwd("/Users/maggiebrown/scProject/CD_Multiome/Multiome_07_2022/CellChat/CD4_CD12")

#### Make CellChat objects; Crohns #### 
crohns.input <- GetAssayData(crohns_ibd_seurat, assay = "RNA", slot = "data") # normalized data matrix
crohns_ibd_seurat <- SetIdent(crohns_ibd_seurat, value = "Celltypes")
labels <- Idents(crohns_ibd_seurat)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
crohns_cellchat <- createCellChat(object = crohns.input, meta = meta, group.by = "group")
crohns_cellchat <- addMeta(crohns_cellchat, meta = meta, meta.name = "group")
crohns_cellchat <- setIdent(crohns_cellchat, ident.use = "group") # set "labels" as default cell identity
levels(crohns_cellchat@idents) # show factor levels of the cell labels
# cellchat DB
CellChatDB <- CellChatDB.human 
crohns_cellchat@DB <- CellChatDB
# subset the expression data of signaling genes for saving computation cost
crohns_cellchat <- subsetData(crohns_cellchat) 
crohns_cellchat <- identifyOverExpressedGenes(crohns_cellchat)
crohns_cellchat <- identifyOverExpressedInteractions(crohns_cellchat)
# compute communication probability
crohns_cellchat <- computeCommunProb(crohns_cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
crohns_cellchat <- filterCommunication(crohns_cellchat, min.cells = 15)
crohns_cellchat <- computeCommunProbPathway(crohns_cellchat)
# aggregate cell-cell network
crohns_cellchat <- aggregateNet(crohns_cellchat)
# viz all cell-cell communications
groupSize <- as.numeric(table(crohns_cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf(file = "cd4_cd12_aggregate_interactions.pdf")
netVisual_circle(crohns_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(crohns_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
# see pathways:
crohns_cellchat@netP$pathways
# [1] "TGFb"     "APP"      "MHC-I"    "CLEC"     "VISFATIN" "ADGRE5"   "MIF"      "ALCAM"    "CD6"     
# [10] "CD86"     "ITGB2"    "SEMA4"    "THBS"     "CD99"     "PECAM1"   "PTPRM"    "NCAM" 
pdf(file = "cd4_cd12_each_pathway.pdf")
for (i in crohns_cellchat@netP$pathways){
  print(netVisual_aggregate(crohns_cellchat, signaling = i, layout = "circle"))
}
dev.off()

plotGeneExpression(crohns_cellchat, signaling = "TGFb")
# TGFB1, TGFBR1, TGFBR2
plotGeneExpression(crohns_cellchat, signaling = "CLEC")
# CLEC2D, CLEC2B, CD69, KLRB1


#### Make CellChat objects; Crohns #### 
healthy.input <- GetAssayData(healthy_ibd_seurat, assay = "RNA", slot = "data") # normalized data matrix
healthy_ibd_seurat <- SetIdent(healthy_ibd_seurat, value = "Celltypes")
labels <- Idents(healthy_ibd_seurat)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
healthy_cellchat <- createCellChat(object = healthy.input, meta = meta, group.by = "group")
healthy_cellchat <- addMeta(healthy_cellchat, meta = meta, meta.name = "group")
healthy_cellchat <- setIdent(healthy_cellchat, ident.use = "group") # set "labels" as default cell identity
levels(healthy_cellchat@idents) # show factor levels of the cell labels
# cellchat DB
CellChatDB <- CellChatDB.human 
healthy_cellchat@DB <- CellChatDB
# subset the expression data of signaling genes for saving computation cost
healthy_cellchat <- subsetData(healthy_cellchat) 
healthy_cellchat <- identifyOverExpressedGenes(healthy_cellchat)
healthy_cellchat <- identifyOverExpressedInteractions(healthy_cellchat)
# compute communication probability
healthy_cellchat <- computeCommunProb(healthy_cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
healthy_cellchat <- filterCommunication(healthy_cellchat, min.cells = 15)
healthy_cellchat <- computeCommunProbPathway(healthy_cellchat)
# aggregate cell-cell network
healthy_cellchat <- aggregateNet(healthy_cellchat)
# viz all cell-cell communications
groupSize <- as.numeric(table(healthy_cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf(file = "/Users/maggiebrown/scProject/CD_Multiome/Multiome_07_2022/CellChat/healthy/healthy_aggregate_interactions.pdf")
netVisual_circle(healthy_cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(healthy_cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
# see pathways:
healthy_cellchat@netP$pathways
# [1] "MHC-I"    "MHC-II"   "APP"      "MIF"      "CD22"     "CD45"     "ADGRE5"   "CLEC"    
# [9] "GALECTIN" "COLLAGEN" "CD99"     "ITGB2"    "PECAM1"   "ALCAM"    "CD6"      "TGFb"    
# [17] "NEGR"     "VISFATIN" "SEMA4"    "CD86"     "SELPLG"   "BAFF"     "CCL"      "CD23"    
# [25] "ANNEXIN"  "LAMININ"  "GRN"      "SEMA3"    "IL16"     "BTLA"     "NCAM"     "ANGPTL"  
# [33] "CADM"     "MPZ"     
pdf(file = "/Users/maggiebrown/scProject/CD_Multiome/Multiome_07_2022/CellChat/healthy/healthy_each_pathway.pdf")
for (i in healthy_cellchat@netP$pathways){
  print(netVisual_aggregate(healthy_cellchat, signaling = i, layout = "circle"))
}
dev.off()

### Comparison ###
object.list <- list(Crohns = crohns_cellchat, Healthy = healthy_cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("Crohns", "Healthy")) # set factor level
pdf(file = "tgfb_cd_h_violin.pdf", width = "20", height = "20")
plotGeneExpression(cellchat, signaling = "TGFb", split.by = "datasets", colors.ggplot = T)
dev.off()
# Chord diagram
pathways.show <- c("TGFb") 
par(mfrow = c(1,2), xpd=TRUE)
pdf(file = "tgfb_cd_h.pdf", width = "20", height = "20")
for (i in 1:length(object.list)) {
  print(netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i])))
}
dev.off()

# perform differential expression analysis
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "Crohns"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", 
                                       pos.dataset = pos.dataset, 
                                       features.name = features.name, 
                                       only.pos = FALSE, 
                                       thresh.pc = 0.1, 
                                       thresh.fc = 0.1, 
                                       thresh.p = 1)

net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in Crohns
net.up <- subsetCommunication(cellchat, net = net, datasets = "Crohns",ligand.logFC = 0.2, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

# violin plot
ibd_seurat <- SetIdent(ibd_seurat, value = "Status_Donor")

group1_ibd_seurat <- subset(ibd_seurat, idents = c("Crohns_4","Crohns_12",
                                                   "Healthy_1", "Healthy_2", "Healthy_3",
                                                      "Healthy_4", "Healthy_5", "Healthy_6",
                                                      "Healthy_7", "Healthy_8", "Healthy_9"))

new.cluster.ids <- c("Healthy", "Healthy", "Healthy", "Healthy", 
                     "Crohns_Group1", "Healthy", "Crohns_Group1",
                     "Healthy", "Healthy", "Healthy", "Healthy")
names(new.cluster.ids) <- levels(group1_ibd_seurat)
group1_ibd_seurat <- RenameIdents(group1_ibd_seurat, new.cluster.ids)
group1_ibd_seurat <- SetIdent(group1_ibd_seurat, value = "Celltypes")

pdf(file = "cd4_cd12_cellchat_vlnplots.pdf", height = 6, width = 18)
VlnPlot(group1_ibd_seurat, features = "TGFB1", split.by = "Status")
VlnPlot(group1_ibd_seurat, features = "TGFBR1", split.plot = TRUE, split.by = "Status")
VlnPlot(group1_ibd_seurat, features = "TGFBR2", split.plot = TRUE, split.by = "Status")
VlnPlot(group1_ibd_seurat, features = "CLEC2D", split.plot = TRUE, split.by = "Status")
VlnPlot(group1_ibd_seurat, features = "CLEC2B", split.plot = TRUE, split.by = "Status")
VlnPlot(group1_ibd_seurat, features = "CD69", split.plot = TRUE, split.by = "Status")
VlnPlot(group1_ibd_seurat, features = "KLRB1", split.plot = TRUE, split.by = "Status")
dev.off()



