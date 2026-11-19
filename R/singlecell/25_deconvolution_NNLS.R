# -*- coding: utf-8 -*-
# 25_deconvolution_NNLS.R : signature-matrix NNLS deconvolution of GSE77459 bulk
# signature from GSE185222 cell-type marker expression (linear space, per-celltype scaled)
suppressMessages({ library(Seurat); library(dplyr) })

sc_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc"
res_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"
obj <- readRDS(file.path(sc_dir, "results", "GSE185222_seurat_annotated.rds"))

# ---- cell-type mean expression in log-normalized space ----
obj <- subset(obj, celltype != "Mixed/low-quality")
data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
ct <- as.character(obj$celltype)
lev <- sort(unique(ct))
avg_log <- sapply(lev, function(g) Matrix::rowMeans(data_mat[, ct == g, drop = FALSE]))
colnames(avg_log) <- lev
avg_lin <- expm1(avg_log)                 # inverse of log1p normalization

# ---- marker-based signature: top 200 specificity genes per cell type ----
spec <- sweep(avg_lin, 1, rowSums(avg_lin) + 1e-8, "/")
markers <- unique(unlist(lapply(lev, function(g) {
  o <- order(spec[, g], decreasing = TRUE)
  rownames(spec)[o[1:200]]
})))
sig <- avg_lin[markers, , drop = FALSE]
sig <- apply(sig, 2, function(x) x / sum(x))   # each celltype column sums to 1
colnames(sig) <- lev
cat("Signature:", nrow(sig), "marker genes x", ncol(sig), "cell types\n")

# ---- bulk ----
expr <- as.matrix(read.csv(file.path(res_dir, "GSE77459_expr_gene.csv"), row.names = 1, check.names = FALSE))
groups <- read.csv(file.path(res_dir, "GSE77459_groups.csv"))
bulk_lin <- 2^expr
common <- intersect(rownames(sig), rownames(bulk_lin))
cat("Common genes:", length(common), "\n")
sig <- sig[common, , drop = FALSE]; bulk_lin <- bulk_lin[common, , drop = FALSE]

# ---- NNLS per sample ----
snames <- colnames(bulk_lin)
props <- sapply(snames, function(s) {
  y <- bulk_lin[, s]; y <- y / sum(y)
  fit <- nnls::nnls(sig, y)
  w <- fit$x; if (sum(w) == 0) w <- rep(1, length(w))
  w / sum(w)
})
props <- t(props)
rownames(props) <- snames
colnames(props) <- lev
write.csv(props, file.path(sc_dir, "results", "GSE77459_deconv_proportions.csv"))

# ---- compare groups ----
df <- as.data.frame(props)
df$sample <- snames
df$group <- groups$group[match(snames, groups$sample)]
write.csv(df, file.path(sc_dir, "results", "GSE77459_deconv_with_group.csv"), row.names = FALSE)
cat("\nMean proportions by group:\n")
print(df %>% group_by(group) %>% summarise(across(all_of(lev), mean), .groups = "drop") %>% as.data.frame(), digits = 3)

cat("\nWilcoxon tests (Normal vs Pulpitis):\n")
pvals <- c()
for (ctp in lev) {
  x <- df[[ctp]][df$group == "Pulpitis"]; yv <- df[[ctp]][df$group == "Normal"]
  p <- tryCatch(wilcox.test(x, yv, exact = FALSE)$p.value, error = function(e) NA)
  pvals[ctp] <- p
  cat(sprintf("%-30s pulpitis=%.4f normal=%.4f p=%.4g\n", ctp, mean(x), mean(yv), p))
}
write.csv(data.frame(celltype = lev, p_wilcox = as.numeric(pvals)),
          file.path(sc_dir, "results", "GSE77459_deconv_stats.csv"), row.names = FALSE)

# barplot
png(file.path(sc_dir, "results", "GSE77459_deconv_barplot.png"), width = 1800, height = 1100, res = 200)
op <- par(mar = c(5, 14, 3, 6))
mat <- t(props)                      # celltype x sample
bp <- barplot(mat, beside = FALSE, col = rainbow(nrow(mat)), border = NA,
              las = 2, ylab = "Estimated cell-type proportion", xlab = "",
              main = "NNLS deconvolution of GSE77459 (single-cell marker signature)",
              cex.names = 0.9, ylim = c(0, 1))
legend("topright", legend = rownames(mat), fill = rainbow(nrow(mat)), bty = "n", cex = 0.7, inset = c(-0.02, 0))
par(op); dev.off()
cat("DONE\n")
