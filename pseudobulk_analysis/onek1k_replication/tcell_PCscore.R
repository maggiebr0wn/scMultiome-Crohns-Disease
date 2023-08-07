#!/usr/bin/env Rscript

.libPaths("/storage/home/mfisher42/bin/R-4.2.2/library")

library(ggplot2)
library(ggfortify)
library(pheatmap)
library(RColorBrewer)
library(readr)
library(Seurat)

# 05-24-2023
# PC scores from CD multiome subset specific DEGs

# Load Seurat object; use rebuilt seurat object to ensure raw counts are used
onek1k_seurat <- readRDS("/storage/home/mfisher42/scProjects/Onek1k/GEO_data/rawdata_seurat_filtered.rds")
onek1k_seurat <- SetIdent(onek1k_seurat, value = "cell_type")

# Pull out T cells:
tcell_labels <- c(
           "CD4-positive, alpha-beta T cell",
           "CD8-positive, alpha-beta T cell",
           "gamma-delta T cell",
           "regulatory T cell",
           "naive thymus-derived CD4-positive, alpha-beta T cell",
           "naive thymus-derived CD8-positive, alpha-beta T cell",
           "central memory CD4-positive, alpha-beta T cell",
           "effector memory CD4-positive, alpha-beta T cell",
           "central memory CD8-positive, alpha-beta T cell",
           "effector memory CD8-positive, alpha-beta T cell",
           "CD4-positive, alpha-beta cytotoxic T cell",
           "mucosal invariant T cell"
           )

tcells_onek1k <- subset(onek1k_seurat, idents = tcell_labels)

# Load hierarchical clustering assignments for t cells
subset_labels <- read.csv("/storage/home/mfisher42/scProjects/Onek1k/3_subsets/tcell_3subsets.csv")
subset_labels$rownames.donor_cluster. <- sub("T_Cells_", "", subset_labels$rownames.donor_cluster.) # fix donor ids

# subset tcells for donor IDs used in hierarchical clustering
tcells_onek1k <- SetIdent(tcells_onek1k, value = "donor_id")
sub_tcells_onek1k <- subset(tcells_onek1k, idents = subset_labels$rownames.donor_cluster.)

# add subset labels to seurat object
labels <- subset_labels$donor_cluster.cluster
names(labels) <- subset_labels$rownames.donor_cluster.
sub_tcells_onek1k <- RenameIdents(sub_tcells_onek1k, labels)
sub_tcells_onek1k$Subset <- Idents(sub_tcells_onek1k)

# only keep cells with more than 25 for pseudobulk
sub_tcells_onek1k$Pseudobulk <- paste(sub_tcells_onek1k$donor_id, sub_tcells_onek1k$cell_type, sep = "_")

pbs <- unique(sub_tcells_onek1k$Pseudobulk)
keep <- NA
for (pb in pbs){
  if (length(which(sub_tcells_onek1k$Pseudobulk == pb)) >= 25){
    print(TRUE)
    keep <- c(keep, pb)
  } else {
    print(FALSE)
  }
}

keep <- keep[complete.cases(keep)]
sub_tcells_onek1k <- SetIdent(sub_tcells_onek1k, value = "Pseudobulk")
sub_tcells_onek1k <- subset(sub_tcells_onek1k, idents = keep)
sub_tcells_onek1k <- SetIdent(sub_tcells_onek1k, value = "donor_id")

small_tcells <- subset(sub_tcells_onek1k, downsample = 500)

# get raw counts
gex <- small_tcells@assays$RNA@counts

sums <- rowSums(gex)
top15k <- names(head(sort(sums, decreasing = TRUE), 15000))

sub_gex <- gex[rownames(gex) %in% top15k,]

gex_df <- as.data.frame(sub_gex)

# Add pseudobulk name
colnames(gex_df) <- paste("T_Cells", small_tcells$donor_id, sep = "_")

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

#####################
### Subset DEGs ###
#####################

# use DEGs from CD Multiome Subsets

Group1_DEGs <- read_csv("/storage/home/mfisher42/scProjects/CD_Subra/Subset_Specific_DEGs_DARs/All_DEGs_03312023/T_Cells/GroupA_DEGs.csv")
Group1_DEGs <- data.frame(Group1_DEGs)
rownames(Group1_DEGs) <- Group1_DEGs$...1

Group2_DEGs <- read_csv("/storage/home/mfisher42/scProjects/CD_Subra/Subset_Specific_DEGs_DARs/All_DEGs_03312023/T_Cells/GroupB_DEGs.csv")
Group2_DEGs <- data.frame(Group2_DEGs)
rownames(Group2_DEGs) <- Group2_DEGs$...1

Group3_DEGs <- read_csv("/storage/home/mfisher42/scProjects/CD_Subra/Subset_Specific_DEGs_DARs/All_DEGs_03312023/T_Cells/GroupC_DEGs.csv")
Group3_DEGs <- data.frame(Group3_DEGs)
rownames(Group3_DEGs) <- Group3_DEGs$...1

Group4_DEGs <- read_csv("/storage/home/mfisher42/scProjects/CD_Subra/Subset_Specific_DEGs_DARs/All_DEGs_03312023/T_Cells/GroupD_DEGs.csv")
Group4_DEGs <- data.frame(Group4_DEGs)
rownames(Group4_DEGs) <- Group4_DEGs$...1


### Check top 50 significant DEGs for Subset 1 ###
upreg <- Group1_DEGs[Group1_DEGs$avg_log2FC > 0,]
top50 <- rownames(head(upreg[order(upreg$p_val_adj),], 50))
filt_final_counts <- final_counts[rownames(final_counts) %in% top50,]

counts_mat <- t(filt_final_counts)

# PCA
pca_res <- prcomp(t(filt_final_counts), scale. = TRUE)
counts_df <- data.frame(counts_mat)

# add group ID
group_ids <- data.frame(paste("T_Cells", small_tcells$donor_id, sep = "_"), small_tcells$Subset)
colnames(group_ids) <- c("Donor", "Subset")
group_ids <- unique(group_ids[order(group_ids$Donor),]) # make list unique
sorted_group_ids <- group_ids[order(match(group_ids$Donor, rownames(counts_df))),] # order subset ids to match counts_df
counts_df$Group <- sorted_group_ids$Subset

# make plots: PC plot, violin plot
p1 <- autoplot(pca_res, data = counts_df, colour = "Group")

pc_df <- data.frame(pca_res$x, counts_df$Group)
pc_df$counts_df.Group <- as.factor(pc_df$counts_df.Group)

p2 <- ggplot(pc_df, aes(x=counts_df.Group, y=PC1, fill=counts_df.Group)) + labs(x = "Group") +
  theme(legend.position = "none") + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1))

pdf(file = "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset1.pdf", width = 4, height = 3)
p1
p2
dev.off()

# ANOVA
result <- aov(PC1 ~ counts_df.Group, data = pc_df)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset1_anova.txt"
sink(output_file)
summary(result)
sink()

# Tukey Test
posthoc <- TukeyHSD(result)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset1_tukey.txt"
sink(output_file)
posthoc
sink()

### Check top 50 significant DEGs for Subset 2 ###
upreg <- Group2_DEGs[Group2_DEGs$avg_log2FC > 0,]
top50 <- rownames(head(upreg[order(upreg$p_val_adj),], 50))
filt_final_counts <- final_counts[rownames(final_counts) %in% top50,]

counts_mat <- t(filt_final_counts)

# PCA
pca_res <- prcomp(t(filt_final_counts), scale. = TRUE)
counts_df <- data.frame(counts_mat)

# add group ID
group_ids <- data.frame(paste("T_Cells", small_tcells$donor_id, sep = "_"), small_tcells$Subset)
colnames(group_ids) <- c("Donor", "Subset")
group_ids <- unique(group_ids[order(group_ids$Donor),]) # make list unique
sorted_group_ids <- group_ids[order(match(group_ids$Donor, rownames(counts_df))),] # order subset ids to match counts_df
counts_df$Group <- sorted_group_ids$Subset

# make plots: PC plot, violin plot
p1 <- autoplot(pca_res, data = counts_df, colour = "Group")

pc_df <- data.frame(pca_res$x, counts_df$Group)
pc_df$counts_df.Group <- as.factor(pc_df$counts_df.Group)

p2 <- ggplot(pc_df, aes(x=counts_df.Group, y=PC1, fill=counts_df.Group)) + labs(x = "Group") +
  theme(legend.position = "none") + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1))

pdf(file = "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset2.pdf", width = 4, height = 3)
p1
p2
dev.off()

# ANOVA
result <- aov(PC1 ~ counts_df.Group, data = pc_df)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset2_anova.txt"
sink(output_file)
summary(result)
sink()

# Tukey Test
posthoc <- TukeyHSD(result)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset2_tukey.txt"
sink(output_file)
posthoc
sink()

### Check top 50 significant DEGs for Subset 3 ###
upreg <- Group3_DEGs[Group3_DEGs$avg_log2FC > 0,]
top50 <- rownames(head(upreg[order(upreg$p_val_adj),], 50))
filt_final_counts <- final_counts[rownames(final_counts) %in% top50,]

counts_mat <- t(filt_final_counts)

# PCA
pca_res <- prcomp(t(filt_final_counts), scale. = TRUE)
counts_df <- data.frame(counts_mat)

# add group ID
group_ids <- data.frame(paste("T_Cells", small_tcells$donor_id, sep = "_"), small_tcells$Subset)
colnames(group_ids) <- c("Donor", "Subset")
group_ids <- unique(group_ids[order(group_ids$Donor),]) # make list unique
sorted_group_ids <- group_ids[order(match(group_ids$Donor, rownames(counts_df))),] # order subset ids to match counts_df
counts_df$Group <- sorted_group_ids$Subset

# make plots: PC plot, violin plot
p1 <- autoplot(pca_res, data = counts_df, colour = "Group")

pc_df <- data.frame(pca_res$x, counts_df$Group)
pc_df$counts_df.Group <- as.factor(pc_df$counts_df.Group)

p2 <- ggplot(pc_df, aes(x=counts_df.Group, y=PC1, fill=counts_df.Group)) + labs(x = "Group") +
  theme(legend.position = "none") + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1))

pdf(file = "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset3.pdf", width = 4, height = 3)
p1
p2
dev.off()

# ANOVA
result <- aov(PC1 ~ counts_df.Group, data = pc_df)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset3_anova.txt"
sink(output_file)
summary(result)
sink()

# Tukey Test
posthoc <- TukeyHSD(result)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset3_tukey.txt"
sink(output_file)
posthoc
sink()


### Check top 50 significant DEGs for Subset 4 ###
upreg <- Group4_DEGs[Group4_DEGs$avg_log2FC > 0,]
top50 <- rownames(head(upreg[order(upreg$p_val_adj),], 50))
filt_final_counts <- final_counts[rownames(final_counts) %in% top50,]

counts_mat <- t(filt_final_counts)

# PCA
pca_res <- prcomp(t(filt_final_counts), scale. = TRUE)
counts_df <- data.frame(counts_mat)

# add group ID
group_ids <- data.frame(paste("T_Cells", small_tcells$donor_id, sep = "_"), small_tcells$Subset)
colnames(group_ids) <- c("Donor", "Subset")
group_ids <- unique(group_ids[order(group_ids$Donor),]) # make list unique
sorted_group_ids <- group_ids[order(match(group_ids$Donor, rownames(counts_df))),] # order subset ids to match counts_df
counts_df$Group <- sorted_group_ids$Subset

# make plots: PC plot, violin plot
p1 <- autoplot(pca_res, data = counts_df, colour = "Group")

pc_df <- data.frame(pca_res$x, counts_df$Group)
pc_df$counts_df.Group <- as.factor(pc_df$counts_df.Group)

p2 <- ggplot(pc_df, aes(x=counts_df.Group, y=PC1, fill=counts_df.Group)) + labs(x = "Group") +
  theme(legend.position = "none") + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1))

pdf(file = "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset4.pdf", width = 4, height = 3)
p1
p2
dev.off()

# ANOVA
result <- aov(PC1 ~ counts_df.Group, data = pc_df)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset4_anova.txt"
sink(output_file)
summary(result)
sink()

# Tukey Test
posthoc <- TukeyHSD(result)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/pc_scores_tcells_subset4_tukey.txt"
sink(output_file)
posthoc
sink()

#### PERMUTE GROUP ASSIGNMENTS: Negative Control ####

set.seed(123)

### Check top 50 significant DEGs for Subset 1 ###
upreg <- Group1_DEGs[Group1_DEGs$avg_log2FC > 0,]
top50 <- rownames(head(upreg[order(upreg$p_val_adj),], 50))
filt_final_counts <- final_counts[rownames(final_counts) %in% top50,]

counts_mat <- t(filt_final_counts)

# PCA
pca_res <- prcomp(t(filt_final_counts), scale. = TRUE)
counts_df <- data.frame(counts_mat)

# add group ID
group_ids <- data.frame(paste("T_Cells", small_tcells$donor_id, sep = "_"), small_tcells$Subset)
colnames(group_ids) <- c("Donor", "Subset")
group_ids <- unique(group_ids[order(group_ids$Donor),]) # make list unique
sorted_group_ids <- group_ids[order(match(group_ids$Donor, rownames(counts_df))),] # order subset ids to match counts_df
counts_df$Group <- sorted_group_ids$Subset


# make plots: PC plot, violin plot
pc_df <- data.frame(pca_res$x, counts_df$Group)
pc_df$Permute_Group <- sample(pc_df$counts_df.Group)
pc_df$Permute_Group <- as.factor(pc_df$Permute_Group)

p2 <- ggplot(pc_df, aes(x=Permute_Group, y=PC1, fill=Permute_Group)) + labs(x = "Group") +
  theme(legend.position = "none") + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1))

pdf(file = "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset1.pdf", width = 4, height = 3)
p2
dev.off()

# ANOVA
result <- aov(PC1 ~ Permute_Group, data = pc_df)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset1_anova.txt"
sink(output_file)
summary(result)
sink()

# Tukey Test
posthoc <- TukeyHSD(result)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset1_tukey.txt"
sink(output_file)
posthoc
sink()

### Check top 50 significant DEGs for Subset 2 ###
upreg <- Group2_DEGs[Group2_DEGs$avg_log2FC > 0,]
top50 <- rownames(head(upreg[order(upreg$p_val_adj),], 50))
filt_final_counts <- final_counts[rownames(final_counts) %in% top50,]

counts_mat <- t(filt_final_counts)

# PCA
pca_res <- prcomp(t(filt_final_counts), scale. = TRUE)
counts_df <- data.frame(counts_mat)

# add group ID
group_ids <- data.frame(paste("T_Cells", small_tcells$donor_id, sep = "_"), small_tcells$Subset)
colnames(group_ids) <- c("Donor", "Subset")
group_ids <- unique(group_ids[order(group_ids$Donor),]) # make list unique
sorted_group_ids <- group_ids[order(match(group_ids$Donor, rownames(counts_df))),] # order subset ids to match counts_df
counts_df$Group <- sorted_group_ids$Subset

# make plots: PC plot, violin plot

pc_df <- data.frame(pca_res$x, counts_df$Group)
pc_df$Permute_Group <- sample(pc_df$counts_df.Group)
pc_df$Permute_Group <- as.factor(pc_df$Permute_Group)

p2 <- ggplot(pc_df, aes(x=Permute_Group, y=PC1, fill=Permute_Group)) + labs(x = "Group") +
  theme(legend.position = "none") + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1))

pdf(file = "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset2.pdf", width = 4, height = 3)
p2
dev.off()

# ANOVA
result <- aov(PC1 ~ Permute_Group, data = pc_df)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset2_anova.txt"
sink(output_file)
summary(result)
sink()

# Tukey Test
posthoc <- TukeyHSD(result)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset2_tukey.txt"
sink(output_file)
posthoc
sink()

### Check top 50 significant DEGs for Subset 3 ###
upreg <- Group3_DEGs[Group3_DEGs$avg_log2FC > 0,]
top50 <- rownames(head(upreg[order(upreg$p_val_adj),], 50))
filt_final_counts <- final_counts[rownames(final_counts) %in% top50,]

counts_mat <- t(filt_final_counts)

# PCA
pca_res <- prcomp(t(filt_final_counts), scale. = TRUE)
counts_df <- data.frame(counts_mat)

# add group ID
group_ids <- data.frame(paste("T_Cells", small_tcells$donor_id, sep = "_"), small_tcells$Subset)
colnames(group_ids) <- c("Donor", "Subset")
group_ids <- unique(group_ids[order(group_ids$Donor),]) # make list unique
sorted_group_ids <- group_ids[order(match(group_ids$Donor, rownames(counts_df))),] # order subset ids to match counts_df
counts_df$Group <- sorted_group_ids$Subset

# make plots: PC plot, violin plot
pc_df <- data.frame(pca_res$x, counts_df$Group)
pc_df$Permute_Group <- sample(pc_df$counts_df.Group)
pc_df$Permute_Group <- as.factor(pc_df$Permute_Group)

p2 <- ggplot(pc_df, aes(x=Permute_Group, y=PC1, fill=Permute_Group)) + labs(x = "Group") +
  theme(legend.position = "none") + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1))

pdf(file = "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset3.pdf", width = 4, height = 3)
p2
dev.off()

# ANOVA
result <- aov(PC1 ~ Permute_Group, data = pc_df)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset3_anova.txt"
sink(output_file)
summary(result)
sink()

# Tukey Test
posthoc <- TukeyHSD(result)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset3_tukey.txt"
sink(output_file)
posthoc
sink()


### Check top 50 significant DEGs for Subset 4 ###
upreg <- Group4_DEGs[Group4_DEGs$avg_log2FC > 0,]
top50 <- rownames(head(upreg[order(upreg$p_val_adj),], 50))
filt_final_counts <- final_counts[rownames(final_counts) %in% top50,]

counts_mat <- t(filt_final_counts)

# PCA
pca_res <- prcomp(t(filt_final_counts), scale. = TRUE)
counts_df <- data.frame(counts_mat)

# add group ID
group_ids <- data.frame(paste("T_Cells", small_tcells$donor_id, sep = "_"), small_tcells$Subset)
colnames(group_ids) <- c("Donor", "Subset")
group_ids <- unique(group_ids[order(group_ids$Donor),]) # make list unique
sorted_group_ids <- group_ids[order(match(group_ids$Donor, rownames(counts_df))),] # order subset ids to match counts_df
counts_df$Group <- sorted_group_ids$Subset

# make plots: PC plot, violin plot
pc_df <- data.frame(pca_res$x, counts_df$Group)
pc_df$Permute_Group <- sample(pc_df$counts_df.Group)
pc_df$Permute_Group <- as.factor(pc_df$Permute_Group)

p2 <- ggplot(pc_df, aes(x=Permute_Group, y=PC1, fill=Permute_Group)) + labs(x = "Group") +
  theme(legend.position = "none") + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.1))

pdf(file = "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset4.pdf", width = 4, height = 3)
p2
dev.off()

# ANOVA
result <- aov(PC1 ~ Permute_Group, data = pc_df)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset4_anova.txt"
sink(output_file)
summary(result)
sink()

# Tukey Test
posthoc <- TukeyHSD(result)

output_file <- "/storage/home/mfisher42/scProjects/Onek1k/3_subsets/neg_controls_pc_score/pc_scores_tcells_subset4_tukey.txt"
sink(output_file)
posthoc
sink()


