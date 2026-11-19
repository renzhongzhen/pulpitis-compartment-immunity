# -*- coding: utf-8 -*-
# 24_trajectory_GSE185222.R : slingshot trajectory for the odontoblast lineage
suppressMessages({
  library(Seurat); library(slingshot); library(SingleCellExperiment)
  library(ggplot2); library(dplyr); library(tidyr); library(princurve)
})

out_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"
obj <- readRDS(file.path(out_dir, "GSE185222_seurat_annotated.rds"))

keep <- c("Fibroblasts/MSC", "Fibroblasts", "Odontoblast-like (OMD+)", "Odontoblasts (DMP1+/IBSP+)")
sub <- subset(obj, celltype %in% keep)
cat("Cells in lineage subset:", ncol(sub), "\n")
sub <- RunUMAP(sub, reduction = "harmony", dims = 1:20, reduction.name = "umap2", verbose = FALSE)

counts <- GetAssayData(sub, assay = "RNA", layer = "counts")
sce <- SingleCellExperiment(assays = list(counts = counts), colData = sub@meta.data)
reducedDims(sce) <- list(UMAP = Embeddings(sub, "umap2"))
sce <- slingshot(sce, clusterLabels = "celltype", reducedDim = "UMAP",
                 start.clus = "Fibroblasts/MSC", approx_points = 150)

pt <- slingPseudotime(sce)
cat("Lineages:", ncol(pt), "\n")
pt_df <- as.data.frame(pt)
pt_df$celltype <- as.character(sce$celltype)
pt_df$condition <- as.character(sce$condition)
write.csv(pt_df, file.path(out_dir, "GSE185222_trajectory_pseudotime.csv"), row.names = FALSE)

# ---- plot curves ----
emb <- as.data.frame(Embeddings(sub, "umap2"))
colnames(emb) <- c("UMAP1", "UMAP2")
emb$celltype <- as.character(sub$celltype)
crv <- slingCurves(sce)
crv_df <- do.call(rbind, lapply(seq_along(crv), function(i) {
  d <- as.data.frame(crv[[i]]$s[crv[[i]]$ord, , drop = FALSE]); colnames(d) <- c("UMAP1","UMAP2"); d$lineage <- paste0("Lineage", i); d
}))
p1 <- ggplot(emb, aes(UMAP1, UMAP2)) + geom_point(aes(color = celltype), size = 0.4, alpha = 0.5) +
  geom_path(data = crv_df, aes(group = lineage), color = "black", linewidth = 1.2) +
  labs(title = "Odontoblast lineage trajectory (slingshot)") + theme_bw() +
  theme(legend.text = element_text(size = 8))

# pseudotime plot (lineage 1)
emb$pt1 <- pt[, 1]
p2 <- ggplot(emb %>% filter(!is.na(pt1)), aes(UMAP1, UMAP2, color = pt1)) +
  geom_point(size = 0.5) + scale_color_viridis_c() +
  labs(title = "Pseudotime (lineage 1)", color = "pseudotime") + theme_bw()
ggsave(file.path(out_dir, "GSE185222_trajectory.png"), p1 + p2, width = 12, height = 5.5, dpi = 200)

# ---- gene trends along pseudotime ----
genes <- c("THY1","PDGFRA","ALPL","OMD","DMP1","IBSP","SPP1","COL1A1","COL1A2","NGFR","MMP9","VDR")
genes <- genes[genes %in% rownames(sce)]
expr <- as.matrix(GetAssayData(sub, assay = "RNA", layer = "data")[genes, ])
tr <- lapply(genes, function(g) {
  d <- data.frame(celltype = as.character(sce$celltype), pt = pt[, 1], expr = expr[g, ]) %>%
        filter(!is.na(pt)) %>% mutate(bin = cut(pt, breaks = 10, labels = FALSE)) %>%
        group_by(bin) %>% summarise(pt = mean(pt), expr = mean(expr), .groups = "drop") %>% mutate(gene = g)
  d
})
tr <- do.call(rbind, tr)
p3 <- ggplot(tr, aes(pt, expr, color = gene)) + geom_line(linewidth = 0.9) + geom_point(size = 1) +
  facet_wrap(~ gene, scales = "free_y", ncol = 4) +
  labs(x = "Pseudotime (lineage 1)", y = "Mean expression (log-normalized)", title = "Gene dynamics along odontoblast lineage") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(out_dir, "GSE185222_trajectory_genes.png"), p3, width = 12, height = 8, dpi = 200)

cat("\nLineage composition:\n"); print(table(sce$celltype, sce$condition))
cat("DONE\n")
