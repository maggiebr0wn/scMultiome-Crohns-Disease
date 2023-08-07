#!/usr/bin/env Rscript

.libPaths("/storage/home/mfisher42/bin/R-4.2.2/library")

library(pheatmap)
library(RColorBrewer)
library(Seurat)

# 05-09-2023
# Perform hierarchical clustering on B Cells GEX

# Load Seurat object; use rebuilt seurat object to ensure raw counts are used
onek1k_seurat <- readRDS("/storage/home/mfisher42/scProjects/Onek1k/GEO_data/rawdata_seurat_filtered.rds")
onek1k_seurat <- SetIdent(onek1k_seurat, value = "cell_type")

# Pull out B cells:

bcell_labels <- c(
           "memory B cell", "naive B cell", "transitional stage B cell",
           "plasmablast"
           )

bcells_onek1k <- subset(onek1k_seurat, idents = bcell_labels)

# only keep cells with more than 25 for pseudobulk
bcells_onek1k$Pseudobulk <- paste(bcells_onek1k$donor_id, bcells_onek1k$cell_type, sep = "_")

pbs <- unique(bcells_onek1k$Pseudobulk)
keep <- NA
for (pb in pbs){
  if (length(which(bcells_onek1k$Pseudobulk == pb)) >= 25){
    print(TRUE)
    keep <- c(keep, pb)
  } else {
    print(FALSE)
  }
}

keep <- keep[complete.cases(keep)]
bcells_onek1k <- SetIdent(bcells_onek1k, value = "Pseudobulk")
onek1k_seurat <- readRDS("/storage/home/mfisher42/scProjects/Onek1k/GEO_data/rawdata_seurat_filtered.rds")
onek1k_seurat <- SetIdent(onek1k_seurat, value = "cell_type")
sub_bcells_onek1k <- subset(bcells_onek1k, idents = keep)
sub_bcells_onek1k <- SetIdent(sub_bcells_onek1k, value = "donor_id")

# downsample 500 cells per donor (currently, too many cells, data too big)

small_bcells <- subset(sub_bcells_onek1k, downsample = 500)

# get raw counts
gex <- small_bcells@assays$RNA@counts

sums <- rowSums(gex)
top15k <- names(head(sort(sums, decreasing = TRUE), 15000))

sub_gex <- gex[rownames(gex) %in% top15k,]

gex_df <- as.data.frame(sub_gex)

# Add pseudobulk name
colnames(gex_df) <- paste("B_Cells", small_bcells$donor_id, sep = "_")

# Sum the total counts for each gene across all cells within the same "pseudobulk"
gene_counts <- sapply(split.default(gex_df, names(gex_df)), rowSums, na.rm = TRUE)

# Divide the count for each gene by the total count for all genes
gene_counts <- data.frame(gene_counts)
pb_totals <- as.data.frame(colSums(gene_counts))
colnames(pb_totals) <- "Total_Counts"

gene_counts_new <- gene_counts

for (col in colnames(gene_counts_new)){
        print(col)
        eval(parse(text = paste("gene_counts_new$",col," <- gene_counts$",col,"/pb_totals['",col,"',]", sep = "")))
}

cpm <- gene_counts_new * 1000000

final_counts <- log2(cpm + 1)

### Hierarchical Clustering: Top 5k genes ###

sums <- rowSums(final_counts)
top5k <- names(head(sort(sums, decreasing = TRUE), 5000))

filt_final_counts <- final_counts[rownames(final_counts) %in% top5k,]
filt_final_counts_mat <- as.matrix(filt_final_counts)

mat <- scale(t(as.matrix(filt_final_counts_mat)))

new_mat <- mat[ , colSums(is.na(mat)) == 0]

p1 <- pheatmap(new_mat, color = colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100),
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               cutree_rows = 3,
               show_colnames = FALSE,
               show_rownames = FALSE,
               clustering_method = "ward.D2")

pdf(file = "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/b_cells_WARDs_pseudobulk_5k_genes.pdf")
p1
dev.off()

# write list of donors per subset:

donor_cluster <- data.frame(cbind(new_mat, cluster = cutree(p1$tree_row, k = 3)))
bcells_subsets <- data.frame(rownames(donor_cluster), donor_cluster$cluster)

write.csv(bcells_subsets, file = "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/bcell_3subsets.csv")
