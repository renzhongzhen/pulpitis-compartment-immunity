suppressPackageStartupMessages({library(Seurat); library(Matrix); library(nnls)})
sc_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc"
res_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"
outdir <- file.path(sc_dir, "results")

obj <- readRDS(file.path(sc_dir, "results", "GSE185222_seurat_annotated.rds"))
obj <- subset(obj, celltype != "Mixed/low-quality")
data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
ct <- as.character(obj$celltype)
lev <- sort(unique(ct))
avg_log <- sapply(lev, function(g) Matrix::rowMeans(data_mat[, ct == g, drop = FALSE]))
colnames(avg_log) <- lev
avg_lin <- expm1(avg_log)

spec <- sweep(avg_lin, 1, rowSums(avg_lin) + 1e-8, "/")

expr <- as.matrix(read.csv(file.path(res_dir, "GSE77459_expr_gene.csv"), row.names = 1, check.names = FALSE))
groups <- read.csv(file.path(res_dir, "GSE77459_groups.csv"))
bulk_lin <- 2^expr
snames <- colnames(bulk_lin)
g <- groups$group[match(snames, groups$sample)]

for (topN in c(100, 200)) {
  markers <- unique(unlist(lapply(lev, function(cc) {
    o <- order(spec[, cc], decreasing = TRUE)
    rownames(spec)[o[seq_len(topN)]]
  })))
  sig <- avg_lin[markers, , drop = FALSE]
  sig <- apply(sig, 2, function(x) x / sum(x))
  common <- intersect(rownames(sig), rownames(bulk_lin))
  cat("topN", topN, "| signature:", length(markers), "genes | common:", length(common), "\n")
  sig <- sig[common, , drop = FALSE]
  b2 <- bulk_lin[common, , drop = FALSE]
  props <- sapply(snames, function(s) {
    y <- b2[, s]; y <- y / sum(y)
    fit <- nnls::nnls(sig, y)
    w <- fit$x; if (sum(w) == 0) w <- rep(1, length(w))
    w / sum(w)
  })
  props <- t(props)
  colnames(props) <- lev
  df <- as.data.frame(props)
  df$sample <- snames
  df$group <- g
  write.csv(df, file.path(outdir, paste0("GSE77459_deconv_top", topN, ".csv")), row.names = FALSE)
  st <- do.call(rbind, lapply(lev, function(cc) {
    pv <- tryCatch(wilcox.test(df[[cc]] ~ df$group, exact = FALSE)$p.value, error = function(e) NA)
    data.frame(celltype = cc,
               mean_normal = mean(df[[cc]][df$group != "Pulpitis"]),
               mean_pulpitis = mean(df[[cc]][df$group == "Pulpitis"]),
               p = signif(pv, 3))
  }))
  write.csv(st, file.path(outdir, paste0("GSE77459_deconv_top", topN, "_stats.csv")), row.names = FALSE)
  cat("== topN =", topN, "==\n"); print(st)
}
cat("\nDONE\n")
