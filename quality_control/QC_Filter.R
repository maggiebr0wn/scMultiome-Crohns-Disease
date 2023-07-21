# 07-06-2022

library(Seurat)
library(ArchR)
library(ggplot2)
library(parallel)

addArchRThreads(threads = 10)
addArchRGenome("hg38")

# This script assesses the QC for the prefiltered projIBD
# Filters out low quality cells

##############################
### 1.0 Load ArchR Project ###
##############################

setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/")

projIBD <- loadArchRProject("Unfiltered-projIBD-07-2022")
# numberOfCells(1): 230791
# medianTSS(1): 4.28
# medianFrags(1): 2578

table(projIBD$Sample)
# Pool_10  Pool_11   Pool_2   Pool_3   Pool_4   Pool_5   Pool_6   Pool_7
# 10580    32795     6776    11102     8690    51976    10218    45546
# Pool_8   Pool_9 Sample_0 Sample_1 Sample_2
# 11486    13288     5086     6215    17033

#########################
### 2.0 Make QC Plots ###
#########################

### TSS Enrichment Ridge Plots

p1 <- plotGroups(
  ArchRProj = projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "ridges"
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

### TSS Enrichment Violin Plots
p2 <- plotGroups(
  ArchRProj = projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

### log10(fragments) Ridge Plots
p3 <- plotGroups(
  ArchRProj = projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "ridges"
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

### log10(fragments) Violin Plots
p4 <- plotGroups(
  ArchRProj = projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics.pdf", ArchRProj = projIBD, addDOC = FALSE, width = 7, height = 5)


### Plot Sample Size Frag Distr and TSS Enrichment Profiles
p1 <- plotFragmentSizes(ArchRProj = projIBD)
p2 <- plotTSSEnrichment(ArchRProj = projIBD)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projIBD, addDOC = FALSE, width = 7, height = 5)


####################################
### 3.0 Filter Low Quality Cells ###
####################################

idxFiltered_Pass <-  which(projIBD$TSSEnrichment >= 8 & # TSS min
                             projIBD$TSSEnrichment <= 30 & # TSS max
                             log10(projIBD$nFrags) >= 3 & # log10(nfrags) min
                             log10(projIBD$nFrags) <= 4.75) # log10(nfrags) max

cellsPass <- projIBD$cellNames[idxFiltered_Pass]

filtered_projIBD <- subsetArchRProject(
  ArchRProj = projIBD,
  cells = cellsPass,
  outputDirectory = "Filtered-projIBD-ATAC-07-2022",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

saveArchRProject(ArchRProj = filtered_projIBD, outputDirectory = "Filtered-projIBD-ATAC-07-2022", load = FALSE)
filtered_projIBD <- loadArchRProject("Filtered-projIBD-ATAC-07-2022")

table(filtered_projIBD$Sample)


############################################
### 4.0 Make QC Plots for Filtered Cells ###
############################################

### TSS Enrichment Ridge Plots

p1 <- plotGroups(
  ArchRProj = filtered_projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "ridges"
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

### TSS Enrichment Violin Plots
p2 <- plotGroups(
  ArchRProj = filtered_projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

### log10(fragments) Ridge Plots
p3 <- plotGroups(
  ArchRProj = filtered_projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "ridges"
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

### log10(fragments) Violin Plots
p4 <- plotGroups(
  ArchRProj = filtered_projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics_FILTERED.pdf", ArchRProj = filtered_projIBD, addDOC = FALSE, width = 7, height = 5)

### Plot Sample Size Frag Distr and TSS Enrichment Profiles
p1 <- plotFragmentSizes(ArchRProj = filtered_projIBD)
p2 <- plotTSSEnrichment(ArchRProj = filtered_projIBD)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile_FILTERED.pdf", ArchRProj = filtered_projIBD, addDOC = FALSE, width = 7, height = 5)

####################################
### 3.0 Filter Low Quality Cells ###
####################################

idxFiltered_Pass <-  which(projIBD$TSSEnrichment >= 8 & # TSS min
                             projIBD$TSSEnrichment <= 30 & # TSS max
                             log10(projIBD$nFrags) >= 3 & # log10(nfrags) min
                             log10(projIBD$nFrags) <= 4.75) # log10(nfrags) max

cellsPass <- projIBD$cellNames[idxFiltered_Pass]

filtered_projIBD <- subsetArchRProject(
  ArchRProj = projIBD,
  cells = cellsPass,
  outputDirectory = "Filtered-projIBD-ATAC-07-2022",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = FALSE
)

saveArchRProject(ArchRProj = filtered_projIBD, outputDirectory = "Filtered-projIBD-ATAC-07-2022", load = FALSE)
filtered_projIBD <- loadArchRProject("Filtered-projIBD-ATAC-07-2022")
# numberOfCells(1): 97593
# medianTSS(1): 17.184
# medianFrags(1): 6882

table(filtered_projIBD$Sample)
# Pool_10  Pool_11   Pool_2   Pool_3   Pool_4   Pool_5   Pool_6   Pool_7
# 7484     5585     6387     7812     7823     6281     6670     7340
# Pool_8   Pool_9 Sample_0 Sample_1 Sample_2
# 11265     6747     4983     5768    13448

############################################
### 4.0 Make QC Plots for Filtered Cells ###
############################################

### TSS Enrichment Ridge Plots

p1 <- plotGroups(
  ArchRProj = filtered_projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "ridges"
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

### TSS Enrichment Violin Plots
p2 <- plotGroups(
  ArchRProj = filtered_projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

### log10(fragments) Ridge Plots
p3 <- plotGroups(
  ArchRProj = filtered_projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "ridges"
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))


### log10(fragments) Violin Plots
p4 <- plotGroups(
  ArchRProj = filtered_projIBD,
  groupBy = "Sample",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
) + theme(axis.text = element_text(size = 10), axis.title = element_text(size = 10))

plotPDF(p1,p2,p3,p4, name = "QC-Sample-Statistics_FILTERED.pdf", ArchRProj = filtered_projIBD, addDOC = FALSE, width = 7, height = 5)

### Plot Sample Size Frag Distr and TSS Enrichment Profiles
p1 <- plotFragmentSizes(ArchRProj = filtered_projIBD)
p2 <- plotTSSEnrichment(ArchRProj = filtered_projIBD)
plotPDF(p1,p2, name = "QC-Sample-FragSizes-TSSProfile_FILTERED.pdf", ArchRProj = filtered_projIBD, addDOC = FALSE, width = 7, height = 5)

