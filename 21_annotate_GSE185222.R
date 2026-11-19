# -*- coding: utf-8 -*-
# 21_annotate_GSE185222.R : marker-based annotation, cluster composition, UMAP plots
suppressMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
})

sc_dir  <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc"
out_dir <- file.path(sc_dir, "results")
obj <- readRDS(file.path(out_dir, "GSE185222_seurat.rds"))
cat("Loaded:", ncol(obj), "cells,", nlevels(obj$seurat_clusters), "clusters\n")

marker_sets <- list(
  DPSC_Stromal = c("THY1","PDGFRA","COL1A1","COL1A2","VIM","MCAM","NGFR","ALPL"),
  Odontoblast  = c("DSPP","DMP1","IBSP","SPP1","OMD","COL10A1"),
  Endothelial  = c("PECAM1","VWF","CLDN5","SELE"),
  Pericyte_SMC = c("ACTA2","PDGFRB","RGS5"),
  Macrophage   = c("CD68","CD14","CSF1R","LYZ","MRC1","CD163","FCGR3A"),
  Neutrophil   = c("FCGR3B","S100A8","S100A9","ELANE","MPO"),
  Tcell        = c("CD3D","CD3E","CD2","IL7R"),
  NK           = c("NKG7","KLRD1","GNLY"),
  Bcell        = c("CD79A","MS4A1","CD19"),
  Plasma       = c("MZB1","JCHAIN","SDC1"),
  DC           = c("CD1C","FCER1A","LAMP3"),
  Mast         = c("TPSAB1","CPA3","KIT"),
  Nerve        = c("PLP1","S100B","MPZ","SOX10"),
  Epithelial   = c("KRT5","KRT14","KRT19"),
  Erythrocyte  = c("HBB","HBA1")
)
all_markers <- unique(unlist(marker_sets))
all_markers <- all_markers[all_markers %in% rownames(obj)]
cat("Markers found:", length(all_markers), "/", length(unique(unlist(marker_sets))), "\n")

avg <- AverageExpression(obj, features = all_markers, group.by = "seurat_clusters",
                         assays = "RNA", layer = "data")$RNA
avg <- as.matrix(avg)
write.csv(avg, file.path(out_dir, "GSE185222_cluster_marker_avg.csv"))

# cluster x sample composition
comp <- as.data.frame(table(cluster = obj$seurat_clusters, sample = obj$sample))
write.csv(comp, file.path(out_dir, "GSE185222_cluster_sample_counts.csv"), row.names = FALSE)
cat("\nCluster x sample:\n"); print(table(obj$seurat_clusters, obj$sample))

# UMAP plots
p1 <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("Clusters")
p2 <- DimPlot(obj, group.by = "sample") + ggtitle("Samples")
ggsave(file.path(out_dir, "GSE185222_UMAP_clusters.png"), p1 + p2, width = 14, height = 6, dpi = 200)

# DotPlot of marker sets
p3 <- DotPlot(obj, features = marker_sets, group.by = "seurat_clusters") + RotatedAxis()
ggsave(file.path(out_dir, "GSE185222_DotPlot_markers.png"), p3, width = 20, height = 6, dpi = 200)

cat("DONE\n")
