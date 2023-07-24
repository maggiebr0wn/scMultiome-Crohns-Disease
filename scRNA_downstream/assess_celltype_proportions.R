# Load the required library

library(Biobase)
library(ggplot2)
library(parallel)
library(scales)
library(dplyr)
library(Matrix)
library(readxl)
library(Seurat)
library(speckle)
library(stringr)
library(RColorBrewer)


### 05-09-23 ###
# Assess proportion of cell populations per subset:
# 1.) B cells
# 2.) T cells

# Load Seurat object
ibd_seurat <- readRDS("/Users/maggiebrown/scProject/CD_Multiome/Multiome_07_2022/seurat_obj_FINAL_LABELS_112022.rds")
ibd_seurat <- SetIdent(ibd_seurat, value = "Coarse_Celltypes")

### Proportions of B cell subpopulations per donor subset ###

# add B cell subset labels:
ibd_seurat <- SetIdent(ibd_seurat, value = "Status_Donor")

bcell_subsets <- c(
  "Group_4", "Group_4", "Group_4", "Group_4", "PossibleIBD_1", "Group_2", "Group_3",
  "Group_4", "Group_4", "Crohns_5", "Group_4", "Group_4", "Group_3", "Group_3",
  "Healthy_6", "Group_3", "Group_1", "Unknown_1", "Group_3", "Group_3", "Group_2",
  "Group_2", "Group_4", "Group_2", "Group_3", "Group_1", "Group_3", "Group_4",
  "Group_2", "Group_2" 
)

names(bcell_subsets) <- levels(ibd_seurat)
ibd_seurat <- RenameIdents(ibd_seurat, bcell_subsets)

ibd_seurat$B_Cell_Subsets <- Idents(ibd_seurat)

# Proportion of B cells per Subset
ibd_seurat <- SetIdent(ibd_seurat, value = "B_Cell_Subsets")
ibd_seurat <- subset(ibd_seurat, idents = c("Group_1", "Group_2", "Group_3", "Group_4"))

ibd_seurat <- SetIdent(ibd_seurat, value = "Celltypes")
bcells_seurat <- subset(ibd_seurat, idents = c("Transitional B Cells", "Resting Naive B Cells", 
                                                   "Activated B Cells", "Plasma B Cells"))

barplot_data <- data.frame(bcells_seurat$B_Cell_Subsets, bcells_seurat$Celltypes, 1)
colnames(barplot_data) <- c("Subset", "Celltypes", "X1")

barplot_data$Sample <- factor(barplot_data$Subset, 
                              levels = c("Group_1", "Group_2", "Group_3", "Group_4"))

# Stacked + percent
p1 <- ggplot(barplot_data, aes(fill=Celltypes, y=X1, x=Sample)) +
  ggtitle("Proportion of B Cell Types per Subset") +
  xlab("Subset Label") +
  ylab("Proportion of Cells") +
  geom_bar(position="fill", stat="identity") +
  theme(text = element_text(size=14)) +
  scale_fill_manual(values = c(
    "dodgerblue", "royalblue3", "blue", "navy"
  )) +
  scale_x_discrete(limits = rev(c(
    "Group_1", "Group_2", "Group_3", "Group_4"
  ))) +
  RotatedAxis()
p1 <- p1 + guides(fill=guide_legend(title = "Cell Types")) + coord_flip()

pdf(file = "/Users/maggiebrown/scProject/CD_Multiome/Multiome_07_2022/Pseudobulk/proportion_bcells_per_subset.pdf", 
    width = 5, height = 3)
p1
dev.off()

### Proportions of T cell subpopulations per donor subset ###
ibd_seurat <- readRDS("/Users/maggiebrown/scProject/CD_Multiome/Multiome_07_2022/seurat_obj_FINAL_LABELS_112022.rds")

# add T cell subset labels:
ibd_seurat <- SetIdent(ibd_seurat, value = "Status_Donor")

tcell_subsets <- c(
  "Group_C", "Group_C", "Group_C", "Group_C", "PossibleIBD_1", "Group_B",
  "Group_B", "Group_C", "Group_D", "Group_D", "Group_C", "Group_C", 
  "Group_C", "Group_D", "Group_B", "Group_B", "Group_A", "Unknown_1",
  "Group_D", "Group_D", "Group_B", "Group_B", "Group_C", "Group_B",
  "Group_D", "Group_A", "Group_D", "Group_D", "Group_B", "Group_B"
)

names(tcell_subsets) <- levels(ibd_seurat)
ibd_seurat <- RenameIdents(ibd_seurat, tcell_subsets)

ibd_seurat$T_Cell_Subsets <- Idents(ibd_seurat)

# Proportion of T cells per Subset
ibd_seurat <- SetIdent(ibd_seurat, value = "T_Cell_Subsets")
ibd_seurat <- subset(ibd_seurat, idents = c("Group_A", "Group_B", "Group_C", "Group_D"))

ibd_seurat <- SetIdent(ibd_seurat, value = "Celltypes")
tcells_seurat <- subset(ibd_seurat, idents = c("Tcm", "Naive CD4+ T Cells", 
                                                   "GdT Cells", "Th1/Th17 Cells",
                                                   "IFN Responding T Cells", "CD8+ Cytotoxic T Cells",
                                                   "MAIT Cells"))

barplot_data <- data.frame(tcells_seurat$T_Cell_Subsets, tcells_seurat$Celltypes, 1)
colnames(barplot_data) <- c("Subset", "Celltypes", "X1")

barplot_data$Sample <- factor(barplot_data$Subset, 
                              levels = c("Group_A", "Group_B", "Group_C", "Group_D"))

# Stacked + percent
p1 <- ggplot(barplot_data, aes(fill=Celltypes, y=X1, x=Sample)) +
  ggtitle("Proportion of T Cell Types per Subset") +
  xlab("Subset Label") +
  ylab("Proportion of Cells") +
  geom_bar(position="fill", stat="identity") +
  theme(text = element_text(size=14)) +
  scale_fill_manual(values = c(
    "darkgreen", "green3", "darkslategrey", "yellow3", 
    "aquamarine2", "lightseagreen", "springgreen1"
  )) +
  scale_x_discrete(limits = rev(c(
    "Group_A", "Group_B", "Group_C", "Group_D"
  ))) +
  RotatedAxis()
p1 <- p1 + guides(fill=guide_legend(title = "Cell Types")) + coord_flip()

pdf(file = "/Users/maggiebrown/scProject/CD_Multiome/Multiome_07_2022/Pseudobulk/proportion_tcells_per_subset.pdf", 
    width = 5, height = 3)
p1
dev.off()

### Use Speckle to assess celltype proportion differences ###

## B cells ##

# Have to subset metadata:
bcell_df <- data.frame(bcells_seurat$Celltypes, bcells_seurat$Status_Donor, bcells_seurat$B_Cell_Subsets)
colnames(bcell_df) <- c("Celltypes", "Donor", "Subset")

bcell_df$Celltypes <- factor(bcell_df$Celltypes, 
                            levels = c("Transitional B Cells", "Resting Naive B Cells", 
                                       "Activated B Cells", "Plasma B Cells"))

keep <- setdiff(unique(bcell_df$Donor), c("Healthy_6", "Crohns_5", "PossibleIBD_1", "Unknown_1"))
bcell_df$Donor <- factor(bcell_df$Donor, levels = keep)

bcell_df$Subset <- factor(bcell_df$Subset, 
                             levels = c("Group_1", "Group_2", "Group_3", "Group_4"))

# Run propeller testing for cell type proportion differences between groups
prop <- propeller(clusters = bcell_df$Celltypes, sample = bcell_df$Donor, 
          group = bcell_df$Subset)

# BaselineProp PropMean.Group_1 PropMean.Group_2 PropMean.Group_3 PropMean.Group_4 Fstatistic
# Activated B Cells       0.08934579      0.962299040       0.05770271        0.0477157       0.03870267  23.908269
# Transitional B Cells    0.54591900      0.004065041       0.44447740        0.5796811       0.53425938  16.405106
# Resting Naive B Cells   0.32373832      0.017987041       0.35943035        0.3559148       0.39031489  10.596129
# Plasma B Cells          0.04099688      0.015648878       0.13838953        0.0166884       0.03672305   1.791559
#                           P.Value          FDR
# Activated B Cells     3.340378e-08 1.336151e-07
# Transitional B Cells  1.452543e-06 2.905086e-06
# Resting Naive B Cells 5.924140e-05 7.898854e-05
# Plasma B Cells        1.693262e-01 1.693262e-01


# exclude group 1

subset_bcell_df <- bcell_df[bcell_df$Subset %in% c("Group_2", "Group_3", "Group_4"),]

subset_bcell_df$Subset <- factor(subset_bcell_df$Subset, 
                                 levels = c("Group_2", "Group_3", "Group_4"))

# Run propeller testing for cell type proportion differences between groups
prop <- propeller(clusters = subset_bcell_df$Celltypes, sample = subset_bcell_df$Donor, 
                  group = subset_bcell_df$Subset)



## T Cells ##

# Have to subset metadata:
tcell_df <- data.frame(tcells_seurat$Celltypes, tcells_seurat$Status_Donor, tcells_seurat$T_Cell_Subsets)
colnames(tcell_df) <- c("Celltypes", "Donor", "Subset")

tcell_df$Celltypes <- factor(tcell_df$Celltypes, 
                             levels = c("Tcm", "Naive CD4+ T Cells", 
                                        "GdT Cells", "Th1/Th17 Cells",
                                        "IFN Responding T Cells", "CD8+ Cytotoxic T Cells",
                                        "MAIT Cells"))

keep <- setdiff(unique(tcell_df$Donor), c("PossibleIBD_1", "Unknown_1"))
tcell_df$Donor <- factor(tcell_df$Donor, levels = keep)

tcell_df$Subset <- factor(tcell_df$Subset, 
                          levels = c("Group_A", "Group_B", "Group_C", "Group_D"))

# Run propeller testing for cell type proportion differences between groups
prop <- propeller(clusters = tcell_df$Celltypes, sample = tcell_df$Donor, 
                  group = tcell_df$Subset)

# BaselineProp PropMean.Group_A PropMean.Group_B PropMean.Group_C
# Th1/Th17 Cells           0.04053587                1       0.00000000    -3.489094e-17
# Tcm                      0.38147508                0       0.37145116     4.878183e-01
# Naive CD4+ T Cells       0.28150389                0       0.31834213     2.905371e-01
# CD8+ Cytotoxic T Cells   0.20740421                0       0.21341354     1.443391e-01
# MAIT Cells               0.03935465                0       0.02978838     5.243624e-02
# GdT Cells                0.03025065                0       0.03213426     2.486933e-02
# IFN Responding T Cells   0.01947566                0       0.03487053    -2.891206e-19
# PropMean.Group_D  Fstatistic      P.Value          FDR
# Th1/Th17 Cells               0.00000000 222.1685740 4.675520e-19 3.272864e-18
# Tcm                          0.28963139  52.8950907 2.135408e-11 7.473928e-11
# Naive CD4+ T Cells           0.30744147  35.8526842 1.566945e-09 3.656206e-09
# CD8+ Cytotoxic T Cells       0.34394882  22.2436380 1.916300e-07 3.353524e-07
# MAIT Cells                   0.03223722   1.6395040 2.037178e-01 2.852049e-01
# GdT Cells                    0.02674110   1.0111878 4.031310e-01 4.703194e-01
# IFN Responding T Cells       0.00000000   0.8445737 4.816072e-01 4.816072e-01

# exclude group A

subset_tcell_df <- tcell_df[tcell_df$Subset %in% c("Group_B", "Group_C", "Group_D"),]

subset_tcell_df$Subset <- factor(subset_tcell_df$Subset, 
                          levels = c("Group_B", "Group_C", "Group_D"))

# Run propeller testing for cell type proportion differences between groups
prop <- propeller(clusters = subset_tcell_df$Celltypes, sample = subset_tcell_df$Donor, 
                  group = subset_tcell_df$Subset)
