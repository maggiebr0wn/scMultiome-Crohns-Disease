# Load the required library

library(Biobase)
library(ggplot2)
library(parallel)
library(scales)
library(dplyr)
library(Matrix)
library(pheatmap)
library(readxl)
library(Seurat)
library(snm)
library(stringr)
library(sva)
library(RColorBrewer)

library(ExpressionNormalizationWorkflow)

### 01-26-22 ###
# Pseudobulk/cluster B Cells with Ward's method

# Load Seurat object
ibd_seurat <- readRDS("/Users/maggiebrown/scProject/CD_Multiome/Multiome_07_2022/seurat_obj_FINAL_LABELS_112022.rds")
ibd_seurat <- SetIdent(ibd_seurat, value = "Coarse_Celltypes")

# only keep cells with more than 25 for pseudobulk
ibd_seurat$Pseudobulk <- paste(ibd_seurat$Celltypes, ibd_seurat$Status_Donor, sep = "_")

pbs <- unique(ibd_seurat$Pseudobulk)
keep <- NA
for (pb in pbs){
  if (length(which(ibd_seurat$Pseudobulk == pb)) >= 25){
    print(TRUE)
    keep <- c(keep, pb)
  } else {
    print(FALSE)
  }
}

keep <- keep[complete.cases(keep)]
ibd_seurat <- SetIdent(ibd_seurat, value = "Pseudobulk")
sub_ibd_seurat <- subset(ibd_seurat, idents = keep)
sub_ibd_seurat <- SetIdent(sub_ibd_seurat, value = "Coarse_Celltypes")

### Sum all counts for each gene, express as cpm ###

## B cells ##
bcells <- subset(sub_ibd_seurat, idents = c("B_Cells"))
# remove H6, CD5, and unknowns
bcells <- SetIdent(bcells, value = "Status_Donor")
keep <- setdiff(unique(bcells$Status_Donor), c("Healthy_6", "Crohns_5", "PossibleIBD_1", "Unknown_1"))
bcells <- subset(bcells, idents = keep)

# get raw counts
gex <- bcells@assays$OriginalRNA@counts
gex_df <- as.data.frame(gex)

#Add pseudbulk name
colnames(gex_df) <- paste(bcells$Coarse_Celltypes, bcells$Status_Donor, sep = "_")

# Sum the total counts for each gene across all cells within the same "pseudobulk"
gene_counts <- sapply(split.default(gex_df, names(gex_df)), rowSums, na.rm = TRUE)

# Divide the count for each gene by the total count for all genes
gene_counts <- data.frame(gene_counts)
pb_totals <- as.data.frame(colSums(gene_counts))
colnames(pb_totals) <- "Total_Counts"

gene_counts_new <- gene_counts
gene_counts_new$B_Cells_Crohns_1 <- gene_counts$B_Cells_Crohns_1/pb_totals["B_Cells_Crohns_1",]
gene_counts_new$B_Cells_Crohns_2 <- gene_counts$B_Cells_Crohns_2/pb_totals["B_Cells_Crohns_2",]
gene_counts_new$B_Cells_Crohns_3 <- gene_counts$B_Cells_Crohns_3/pb_totals["B_Cells_Crohns_3",]
gene_counts_new$B_Cells_Crohns_4 <- gene_counts$B_Cells_Crohns_4/pb_totals["B_Cells_Crohns_4",]
gene_counts_new$B_Cells_Crohns_6 <- gene_counts$B_Cells_Crohns_6/pb_totals["B_Cells_Crohns_6",]
gene_counts_new$B_Cells_Crohns_7 <- gene_counts$B_Cells_Crohns_7/pb_totals["B_Cells_Crohns_7",]
gene_counts_new$B_Cells_Crohns_8 <- gene_counts$B_Cells_Crohns_8/pb_totals["B_Cells_Crohns_8",]
gene_counts_new$B_Cells_Crohns_9 <- gene_counts$B_Cells_Crohns_9/pb_totals["B_Cells_Crohns_9",]
gene_counts_new$B_Cells_Crohns_10 <- gene_counts$B_Cells_Crohns_10/pb_totals["B_Cells_Crohns_10",]
gene_counts_new$B_Cells_Crohns_11 <- gene_counts$B_Cells_Crohns_11/pb_totals["B_Cells_Crohns_11",]
gene_counts_new$B_Cells_Crohns_12 <- gene_counts$B_Cells_Crohns_12/pb_totals["B_Cells_Crohns_12",]
gene_counts_new$B_Cells_Crohns_13 <- gene_counts$B_Cells_Crohns_13/pb_totals["B_Cells_Crohns_13",]
gene_counts_new$B_Cells_Crohns_14 <- gene_counts$B_Cells_Crohns_14/pb_totals["B_Cells_Crohns_14",]
gene_counts_new$B_Cells_Crohns_15 <- gene_counts$B_Cells_Crohns_15/pb_totals["B_Cells_Crohns_15",]
gene_counts_new$B_Cells_Crohns_16 <- gene_counts$B_Cells_Crohns_16/pb_totals["B_Cells_Crohns_16",]
gene_counts_new$B_Cells_Crohns_17 <- gene_counts$B_Cells_Crohns_17/pb_totals["B_Cells_Crohns_17",]
gene_counts_new$B_Cells_Crohns_18 <- gene_counts$B_Cells_Crohns_18/pb_totals["B_Cells_Crohns_18",]
gene_counts_new$B_Cells_Crohns_19 <- gene_counts$B_Cells_Crohns_19/pb_totals["B_Cells_Crohns_19",]

gene_counts_new$B_Cells_Healthy_1 <- gene_counts$B_Cells_Healthy_1/pb_totals["B_Cells_Healthy_1",]
gene_counts_new$B_Cells_Healthy_2 <- gene_counts$B_Cells_Healthy_2/pb_totals["B_Cells_Healthy_2",]
gene_counts_new$B_Cells_Healthy_3 <- gene_counts$B_Cells_Healthy_3/pb_totals["B_Cells_Healthy_3",]
gene_counts_new$B_Cells_Healthy_4 <- gene_counts$B_Cells_Healthy_4/pb_totals["B_Cells_Healthy_4",]
gene_counts_new$B_Cells_Healthy_5 <- gene_counts$B_Cells_Healthy_5/pb_totals["B_Cells_Healthy_5",]
gene_counts_new$B_Cells_Healthy_7 <- gene_counts$B_Cells_Healthy_7/pb_totals["B_Cells_Healthy_7",]
gene_counts_new$B_Cells_Healthy_8 <- gene_counts$B_Cells_Healthy_8/pb_totals["B_Cells_Healthy_8",]
gene_counts_new$B_Cells_Healthy_9 <- gene_counts$B_Cells_Healthy_9/pb_totals["B_Cells_Healthy_9",]

cpm <- gene_counts_new * 1000000

final_counts <- log2(cpm + 1)

### Hierarchical Clustering: Top 5k genes ###

sums <- rowSums(final_counts)
top5k <- names(head(sort(sums, decreasing = TRUE), 5000))

filt_final_counts <- final_counts[rownames(final_counts) %in% top5k,]
filt_final_counts_mat <- as.matrix(filt_final_counts)

mat <- scale(t(as.matrix(filt_final_counts_mat)))
#mat <- t(as.matrix(filt_final_counts_mat))

new_mat <- mat[ , colSums(is.na(mat)) == 0]

p1 <- pheatmap(new_mat, color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100),
               cluster_rows = TRUE, 
               cluster_cols = TRUE, 
               cutree_rows = 4,
               show_colnames = FALSE,
               clustering_method = "ward.D2")

pdf(file = "/Users/maggiebrown/scProject/CD_Multiome/Multiome_07_2022/Pseudobulk/hierarchical_clustering_02062023/b_cells_WARDs_pseudobulk_withHealthy_5k_genes.pdf")
p1
dev.off()
