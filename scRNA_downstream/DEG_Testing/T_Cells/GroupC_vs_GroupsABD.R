#!/usr/bin/env Rscript

.libPaths("/storage/home/mfisher42/bin/R-4.2.2/library")

library(MAST)
library(Seurat)

# 03-31-2023
# Run MAST on Group3 vs Groups 1,2,4

##================================================================================##
##                     Detecting DEGs using MAST method                           ##
##  It can't test interaction terms, eg. Disease condition * Tissue type,         ##
##  can only test if gene expression is significantly different between two       ##
##  groups, eg. disease vs control, tissue A vs B                                 ##
##================================================================================##


MAST_runfun <- function(object,assay = "RNA",slot = "data", features = NULL, ident.1 = NULL, ident.2 = NULL, fc.name = NULL,
                        cov = NULL,min.cells.group = 3, reduction = NULL,pseudocount.use = 1,method = "bayesglm",base = 2){
  ## cov: selected covariance variables for fitting in the hurdle model, which need to have the same name as
  ##      in Seurat object meta data slot.
  ## base: used to set base of log transformation, default value is 2, alternative is exp(1)
  ## method:

  ## Getting the cell ids in ident.1 and ident.2,respectively
  cells.1 = row.names(object@meta.data)[which(Idents(object) %in% ident.1)]
  cells.2 = row.names(object@meta.data)[which(Idents(object) %in% ident.2)]
  if(is.null(features)){
          features = row.names(object)
  }
  ## error checking
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", cells.1)
  } else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", cells.2)
    return(NULL)
  } else if (length(x = cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, " cells")
  } else if (length(x = cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, " cells")
  } else if (any(!cells.1 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.1) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.1 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  } else if (any(!cells.2 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.2) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.2 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  }


  ## constructing group information datasframe
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])

  ## Adding customed covariables
  if(!is.null(cov)){
    group.info <- cbind(group.info, object@meta.data[row.names(group.info),cov])
    colnames(group.info) <- c("group",cov)

  }else{
    colnames(group.info) <- "group"
  }

  ##  creating fomula for the hurdle model, but bayesglm fomula is slightly different from glmer fomula to set up random effect.
  if(method=="bayesglm"){
    fmla <- as.formula(object = paste0(" ~ cngeneson + ", paste(colnames(group.info), collapse = "+")))
  }else if(method=="glmer"){
    if(ncol(group.info)>1){
      if(ncol(group.info)==2){
        fmla <- as.formula(
          object = paste0(" ~ group + cngeneson + (1|",colnames(group.info)[2],")"))

      }else{
        part2 = c()
        for (col.index in c(2:ncol(group.info))){
          part2 = paste0(part2," + (1|",colnames(group.info)[col.index],")")
        }
        fmla <- as.formula(object = paste0(" ~ group + cngeneson ",part2))
      }
    }else{
      fmla <- as.formula(object = paste0(" ~ group + cngeneson"))
    }
  }


  ## Constructing single cell dataset
  fdata = data.frame(primerid = features)
  fit_data <- MAST::FromMatrix(
    exprsArray = as.matrix(object[[assay]]@data[features,row.names(group.info)]),
    cData = group.info,
    fData = fdata
  )

  # recalculating the cellular detection rate
  cdr2 <-colSums(SummarizedExperiment::assay(fit_data)>0)
  colData(fit_data)$cngeneson <- scale(cdr2)

  data = as.matrix(object[[assay]]@data)

  thresh.min = 0
  pct.1 <- round(
   rowSums(data[features, cells.1, drop = FALSE] > thresh.min) /
      length(cells.1),
    digits = 3
  )

  pct.2 <- round(
    rowSums(data[features, cells.2, drop = FALSE] > thresh.min) /
      length(cells.2),
    digits = 3
  )

  # feature selection (based on average difference)
  # if selecting data slot is scale.data, then calculating the rowMeans as avg gene expression
  # if selecting data slot is data, then using formula log(x = rowMeans(x = expm1(x = x)) + pseudocount.use)
  # if selecting data slot is count, then using formula log(x = rowMeans(x = x) + pseudocount.use)
  mean.fxn <- if (is.null(x = reduction) && slot != "scale.data") {
    switch(
      EXPR = slot,
      'data' = function(x) {
        return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use, base = base))
      },
      function(x) {
        return(log(x = rowMeans(x = x) + pseudocount.use, base = base))
      }
    )
  } else {
    rowMeans
  }

  base.text <- ifelse(
    test = base == exp(1),
    yes = "",
    no = base
  )
  fc.name <- fc.name %||% ifelse(
    test = slot == "scale.data",
    yes = "avg_diff",
    no = paste0("avg_log", base.text, "FC")
  )

  ## Calculating the avg_logFC between two groups
  data.1 <- mean.fxn(data[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(data[features, cells.2, drop = FALSE])
  total.diff <- (data.1 - data.2)

  ## Constructing zlm model and using Group1 as reference group
  cond<-factor(colData(fit_data)$group)
  cond<-relevel(cond,"Group1")
  colData(fit_data)$group<-cond
  # fmla <- as.formula(
  #   object = paste0(" ~ ", paste(colnames(group.info), collapse = "+"))
  # )
  if(method == "glmer"){
    zlmCond <- MAST::zlm(formula = fmla, fit_data,method = method,ebayes = F,strictConvergence = FALSE,fitArgsD = list(nAGQ = 0))
  }else{
    zlmCond <- MAST::zlm(formula = fmla, fit_data,method = method)
  }

  ## Running a likelihood ratio test here, testing for differences when we drop the group factor
  ## only test the group coefficient.
  summaryCond <- summary(object = zlmCond, doLRT = 'groupGroup2')
  summaryDt <- summaryCond$datatable
  # fcHurdle <- merge(
  #   summaryDt[contrast=='groupGroup2' & component=='H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
  #   summaryDt[contrast=='groupGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid'
  # ) #logFC coefficients
  # fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle <- merge(summaryDt[contrast=='groupGroup2' & component=='H',.(primerid, `Pr(>Chisq)`)],
                    summaryDt[contrast=='groupGroup2' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #hurdle P values
  fcHurdle = fcHurdle[features,]
  p_val <- fcHurdle[,`Pr(>Chisq)`]
  p_val_adj = p.adjust(
    p = p_val,
    method = "bonferroni",
    n = nrow(object)
  )
  Coef <- fcHurdle[,`coef`]
  genes.return <- fcHurdle[,`primerid`]

  to.return <- data.frame(pct.1 = pct.1, pct.2 = pct.2, Coef = Coef, avg_logFC = total.diff, p_val = p_val, p_val_adj = p_val_adj, row.names = genes.return)
  names(to.return)[4] = fc.name
  return(to.return)
}

`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

### end of function ###

### Load Seurat object and run MAST ###
ibd_seurat <- readRDS("/storage/home/mfisher42/scProjects/CD_Subra/Seurat_Object/seurat_obj_FINAL_LABELS_with_PEAKS_122022.rds")
ibd_seurat <- SetIdent(ibd_seurat, value = "Coarse_Celltypes")

DefaultAssay(ibd_seurat) <- "OriginalRNA"
ibd_seurat <- NormalizeData(ibd_seurat)

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

## Get T cells ##
tcells <- subset(sub_ibd_seurat, idents = c("T_Cells"))

# remove unknowns
tcells <- SetIdent(tcells, value = "Status_Donor")
keep <- setdiff(unique(tcells$Status_Donor), c("PossibleIBD_1", "Unknown_1"))
tcells <- subset(tcells, idents = keep)

ident1 <- c("Crohns_1", "Crohns_18", "Crohns_11", # GroupC
	    "Crohns_19", "Crohns_14", "Healthy_7",
	    "Crohns_9", "Crohns_16", "Healthy_9")

ident2 <- setdiff(unique(tcells$Status_Donor), ident1)

degs <- MAST_runfun(tcells, assay = "OriginalRNA",
                    ident.1 = ident1,
                    ident.2 = ident2,
                    cov = "nCount_OriginalRNA"
                    )

# only keep DEGs passing adj p-value threshold
sig_degs <- degs[degs$p_val_adj <= 0.05,]
sig_degs <- sig_degs[abs(sig_degs$avg_log2FC) >= 0.25,]

# filter for pct cells expressing > 10%
sig_degs <- sig_degs[sig_degs$pct.1 >= 0.1 | sig_degs$pct.2 >= 0.1,]

setwd("/storage/home/mfisher42/scProjects/CD_Subra/Subset_Specific_DEGs_DARs/All_DEGs_03312023/T_Cells")

write.csv(sig_degs, file = "GroupC_DEGs.csv")






