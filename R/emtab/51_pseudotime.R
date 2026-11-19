suppressPackageStartupMessages({
  library(Seurat); library(slingshot); library(Matrix)
  library(org.Hs.eg.db); library(AnnotationDbi)
  library(ggplot2); library(pheatmap); library(RColorBrewer)
})
sc_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc"
outdir <- file.path(sc_dir, "results")
p <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/ea960c355a014573b0e904afc349b633processedmatrix.rds"
x <- readRDS(p)
md <- x@meta.data
hs <- factor(as.character(md$Health_state), levels = c("Healthy", "Early Stage", "Advanced Stage"))
rd <- Embeddings(x, reduction = "umap")
colnames(rd) <- c("UMAP_1", "UMAP_2")

set.seed(1)
sds <- slingshot(rd, clusterLabels = hs, start.clus = "Healthy", reducedDim = "UMAP", approx_points = 400)
n_line <- ncol(slingPseudotime(sds))
cat("lineages:", n_line, "\n")
pt <- slingPseudotime(sds, na = FALSE)
pt_avg <- rowMeans(pt, na.rm = TRUE)
pt_avg[!is.finite(pt_avg)] <- NA
md$pseudotime <- pt_avg
write.csv(data.frame(cell = rownames(md), stage = hs, pseudotime = pt_avg),
          file.path(outdir, "EMTAB_pseudotime.csv"), row.names = FALSE)

ok <- !is.na(pt_avg)
kw <- kruskal.test(pt_avg[ok] ~ hs[ok])
cat("pseudotime by Health_state: Kruskal p =", kw$p.value, "\n")
m <- tapply(pt_avg[ok], hs[ok], median, na.rm = TRUE)
cat("median pseudotime:", paste(names(m), round(m, 2), collapse = ", "), "\n")
rho <- cor(as.numeric(hs[ok]), pt_avg[ok], method = "spearman")
cat("Spearman rho (stage order vs pseudotime):", round(rho, 3), "\n")

# UMAP colored by pseudotime
df <- data.frame(UMAP_1 = rd[, 1], UMAP_2 = rd[, 2], pt = pt_avg)
df$stage <- hs
set.seed(2)
idx <- sample(which(ok), 30000)
p1 <- ggplot(df[idx, ], aes(x = UMAP_1, y = UMAP_2, color = pt)) +
  geom_point(size = 0.12) + scale_color_viridis_c(option = "C", na.value = "grey85") +
  labs(title = "UMAP colored by slingshot pseudotime (start: Healthy)", color = "pseudotime") +
  theme_bw(base_size = 10) + theme(legend.position = "right")
ggsave(file.path(outdir, "EMTAB_pseudotime_umap.png"), p1, width = 7, height = 5.5, dpi = 160)

p2 <- ggplot(df[idx, ], aes(x = UMAP_1, y = UMAP_2, color = stage)) +
  geom_point(size = 0.12) + scale_color_manual(values = c(Healthy = "#4C72B0", `Early Stage` = "#DD8452", `Advanced Stage` = "#C44E52")) +
  labs(title = "UMAP by disease stage", color = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "right")
ggsave(file.path(outdir, "EMTAB_pseudotime_umap_stage.png"), p2, width = 7, height = 5.5, dpi = 160)

# gene dynamics along pseudotime
genes <- c("SELE", "ACKR1", "CCL2", "CCL5", "CXCL8", "IL1B", "STAT1", "NFKB1", "SPP1", "COL1A1", "DCN", "DSPP", "NES", "MYC")
cnt <- LayerData(x, assay = "RNA", layer = "counts")
ens <- sub("\\..*$", "", rownames(cnt))
sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
keep <- !is.na(sym) & !duplicated(sym)
cnt2 <- cnt[keep, , drop = FALSE]; rownames(cnt2) <- sym[keep]
cs <- colSums(cnt2); cs[cs == 0] <- 1
genes <- genes[genes %in% rownames(cnt2)]
Esub <- as.matrix(cnt2[genes, , drop = FALSE])
Elog <- log2(sweep(Esub, 2, cs, "/") * 1e4 + 1)
rm(cnt, cnt2); gc()

bins <- cut(pt_avg, breaks = quantile(pt_avg[ok], probs = seq(0, 1, 0.1)), include.lowest = TRUE)
bins_ok <- !is.na(bins)
bin_means <- do.call(cbind, lapply(levels(bins), function(b) {
  idx <- which(bins == b)
  if (length(idx) == 0) return(rep(NA, nrow(Elog)))
  rowMeans(Elog[, idx, drop = FALSE])
}))
colnames(bin_means) <- paste0("B", seq_len(ncol(bin_means)))
z <- t(apply(bin_means, 1, function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)))
z[!is.finite(z)] <- 0
write.csv(data.frame(gene = rownames(z), z), file.path(outdir, "EMTAB_pseudotime_gene_trend.csv"), row.names = FALSE)

annot_col <- data.frame(stage = factor(c(rep("Healthy", 2), rep("Early", 4), rep("Advanced", 4))))
rownames(annot_col) <- colnames(z)
ann_colors <- list(stage = c(Healthy = "#4C72B0", Early = "#DD8452", Advanced = "#C44E52"))
png(file.path(outdir, "EMTAB_pseudotime_gene_trend.png"), width = 1800, height = 1200, res = 160)
pheatmap(z, cluster_cols = FALSE, cluster_rows = FALSE, color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
         annotation_col = annot_col, annotation_colors = ann_colors, main = "Gene dynamics along pseudotime (10 bins, Healthy to Advanced)")
dev.off()
cat("\nDONE\n")
