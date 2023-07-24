#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(Dune)
library(ggplot2)
library(Rmagic)
library(Seurat)
library(stringr)

# 07-15-2022
# Compute RandIndex for Harmonyv2 ONLY
# 1.) ATAC only, with Harmony (n = 108 total)  
# 2.) ATAC only, without Harmony (n = 36 total)
# 3.) ATAC+RNA, with Harmony (n = 108 total) -- Only this for this script
# 4.) ATAC+RNA, without Harmony (n = 36 total)
# Plots are generate separately

# Previously, was computing RandIndex for all 108 instances -- this is taking waaay too long.
# So, will be computing for each version of Harmony separately.

##############################################
### 1.0 ATAC+GEX, with Harmonyv2 (n = 36) ###
##############################################
setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/Filtered-Final-projIBD-ATACandGEX-07-2022/CLUSTS_UMAPS/Harmony_Combo")

### 1.1 load in data, build a dataframe ###

files <- list.files(path = "/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/Filtered-Final-projIBD-ATACandGEX-07-2022/CLUSTS_UMAPS/Harmony_Combo",
                    pattern = ".csv")

# Only keep Harmonyv1 Files:
files <- str_subset(files, "Harmonyv2", negate = FALSE)

data <- do.call(cbind, sapply(files,data.table::fread, simplify = FALSE))
# clean dataframe
data = data[-1,]
data <- as.data.frame(data)
rownames(data) <- data[,1]
df_new <- data %>% select(-contains(".V1"))
colnames(df_new) <- gsub('.{7}$', '', colnames(df_new))

### 1.2 do RandIndex ###

atac_df <- lapply(df_new, gsub, pattern="C", replacement="")
atac_df_final <- as.data.frame(atac_df)

merger <- Dune(clusMat = atac_df_final, verbose = TRUE)

setwd("/storage/home/mfisher42/scProjects/CD_Subra/ArchR_Project/Multiome_07_2022/ArchR_07_2022/Filtered-Final-projIBD-ATACandGEX-07-2022/CLUSTS_UMAPS/Harmony_Combo/RandIndex/Harmonyv2")

write.csv(merger$initialMat, file = "merger_initialMat.csv")
write.csv(merger$currentMat, file = "merger_currentMat.csv")
write.csv(merger$merges, file = "merger_merges.csv")
write.csv(merger$ImpMetric, file = "merger_ImpMetric.csv")
write.csv(merger$metric, file = "merger_metric.csv")

