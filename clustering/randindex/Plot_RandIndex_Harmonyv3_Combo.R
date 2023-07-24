#!/usr/bin/env Rscript

library(data.table)  
library(dplyr) 
library(Dune)
library(ggplot2)
library(Rmagic)
library(Seurat)
library(stringr)

# 07-13-2022
# Plot RandIndex for:
# 1.) ATAC only, with Harmony (n = 108) 
# 2.) ATAC only, without Harmony (n = 36)
# 3.) ATAC+RNA, with Harmony (n = 108) -- Only this for this script
# 4.) ATAC+RNA, without Harmony (n = 36)  
# Plots are generate separately

########################################
### 1.0 ATAC+RNA, Harmonyv3 (n = 38) ###
########################################
setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/Filtered-Final-projIBD-ATACandGEX-07-2022/CLUSTS_UMAPS/Harmony_Combo/RandIndex/Harmonyv3")

merger <- read.csv("merger_currentMat.csv", row.names=1)

### ALL Plotted Together ###
# plot ARIs as barplot
getARIs <- ARIs(merger)

sumdata=data.frame(value=apply(getARIs,2,mean))
sumdata$key=rownames(sumdata)

# save plots
pdf(file = "RandIndex_heatmap_all.pdf", height = 15, width = 15)
plotARIs(clusMat = merger) + RotatedAxis()
dev.off()

pdf(file = "RandIndex_barplot_all.pdf", height = 8, width = 8)
ggplot(data=sumdata, aes(x=key, y=value, fill=key)) +
  geom_bar(colour="black", stat="identity") + coord_flip() + 
  theme(legend.position="none")
dev.off()

### ALL RES 0.2 Plotted Together ###
merger02 <- merger[ , grepl( "Res02" , names(merger) ) ]
# plot ARIs as barplot
getARIs <- ARIs(merger02)

sumdata=data.frame(value=apply(getARIs,2,mean))
sumdata$key=rownames(sumdata)

# save plots
pdf(file = "RandIndex_heatmap_Res02.pdf", height = 6, width = 7)
plotARIs(clusMat = merger02) + RotatedAxis()
dev.off()

pdf(file = "RandIndex_barplot_Res02.pdf", height = 5, width = 5)
ggplot(data=sumdata, aes(x=key, y=value, fill=key)) +
  geom_bar(colour="black", stat="identity") + coord_flip() + 
  theme(legend.position="none")
dev.off()

### ALL RES 0.5 Plotted Together ###
merger05 <- merger[ , grepl( "Res05" , names(merger) ) ]
# plot ARIs as barplot
getARIs <- ARIs(merger05)

sumdata=data.frame(value=apply(getARIs,2,mean))
sumdata$key=rownames(sumdata)

# save plots
pdf(file = "RandIndex_heatmap_Res05.pdf", height = 6, width = 7)
plotARIs(clusMat = merger05) + RotatedAxis()
dev.off()

pdf(file = "RandIndex_barplot_Res05.pdf", height = 5, width = 5)
ggplot(data=sumdata, aes(x=key, y=value, fill=key)) +
  geom_bar(colour="black", stat="identity") + coord_flip() + 
  theme(legend.position="none")
dev.off()

### ALL RES 0.8 Plotted Together ###
merger08 <- merger[ , grepl( "Res08" , names(merger) ) ]
# plot ARIs as barplot
getARIs <- ARIs(merger08)

sumdata=data.frame(value=apply(getARIs,2,mean))
sumdata$key=rownames(sumdata)

# save plots
pdf(file = "RandIndex_heatmap_Res08.pdf", height = 6, width = 7)
plotARIs(clusMat = merger08) + RotatedAxis()
dev.off()

pdf(file = "RandIndex_barplot_Res08.pdf", height = 5, width = 5)
ggplot(data=sumdata, aes(x=key, y=value, fill=key)) +
  geom_bar(colour="black", stat="identity") + coord_flip() + 
  theme(legend.position="none")
dev.off()

### ALL RES 1.0 Plotted Together ###
merger1 <- merger[ , grepl( "Res1" , names(merger) ) ]
# plot ARIs as barplot
getARIs <- ARIs(merger1)

sumdata=data.frame(value=apply(getARIs,2,mean))
sumdata$key=rownames(sumdata)

# save plots
pdf(file = "RandIndex_heatmap_Res10.pdf", height = 6, width = 7)
plotARIs(clusMat = merger1) + RotatedAxis()
dev.off()

pdf(file = "RandIndex_barplot_Res10.pdf", height = 5, width = 5)
ggplot(data=sumdata, aes(x=key, y=value, fill=key)) +
  geom_bar(colour="black", stat="identity") + coord_flip() + 
  theme(legend.position="none")
dev.off()

