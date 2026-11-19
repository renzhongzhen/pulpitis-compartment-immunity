# -*- coding: utf-8 -*-
# 20_prep_GSE185222.R : 读 4 个 CSV (gene x cell, normalized values) -> Seurat -> Harmony -> cluster -> UMAP
suppressMessages({
  library(Seurat); library(data.table); library(Matrix); library(harmony); library(dplyr)
})

sc_dir  <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/GSE185222"
out_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- c(DTP_01 = "GSM5608427_DTP_01.csv.gz",
           DTP_02 = "GSM5608428_DTP_02.csv.gz",
           DTP_04 = "GSM5608429_DTP_04.csv.gz",
           DTP_05 = "GSM5608430_DTP_05.csv.gz")

mats <- list()
for (s in names(files)) {
  dt <- fread(file.path(sc_dir, files[[s]]), showProgress = FALSE)
  genes <- make.unique(as.character(dt[[1]]))
  m <- as.matrix(dt[, -1, with = FALSE])
  rownames(m) <- genes
  mats[[s]] <- as(m, "dgCMatrix")
  cat(sprintf("%s: %d genes x %d cells\n", s, nrow(mats[[s]]), ncol(mats[[s]])))
  rm(dt, m); gc(verbose = FALSE)
}

common <- Reduce(intersect, lapply(mats, rownames))
cat("Common genes:", length(common), "\n")
mats <- lapply(names(mats), function(s) {
  m <- mats[[s]][common, , drop = FALSE]
  colnames(m) <- paste0(s, "_", colnames(m))
  m
})
big <- do.call(cbind, unname(mats))
cat("Merged matrix:", nrow(big), "genes x", ncol(big), "cells\n")

obj <- CreateSeuratObject(counts = big, project = "GSE185222")
obj[["RNA"]]$data <- obj[["RNA"]]$counts   # values already normalized by authors
obj$sample <- sub("^(DTP_\\d+)_.*$", "\\1", colnames(obj))
cat("Sample labels:", paste(unique(obj$sample), collapse = ", "), "\n")
cat("Seurat object created\n")

obj <- FindVariableFeatures(obj, selection.method = "mvp", nfeatures = 2000, verbose = FALSE)
obj <- ScaleData(obj, verbose = FALSE)
obj <- RunPCA(obj, npcs = 30, verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "sample", verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20, verbose = FALSE)

saveRDS(obj, file.path(out_dir, "GSE185222_seurat.rds"))
cat("\nCluster x sample table:\n")
print(table(obj$seurat_clusters, obj$orig.ident))
cat("\nCells per sample:\n"); print(table(obj$orig.ident))
cat("DONE\n")
