suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(DESeq2); library(fgsea)
  library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2); library(pheatmap)
})
args <- commandArgs(trailingOnly = TRUE)
supplement_mode <- length(args) > 0 && nzchar(args[1])
p <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/ea960c355a014573b0e904afc349b633processedmatrix.rds"
outdir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"
gmt_file <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/deep/h.all.v2024.1.Hs.symbols.gmt"

x <- readRDS(p)
md <- x@meta.data
cnt <- LayerData(x, assay = "RNA", layer = "counts")
cts <- levels(Idents(x))
if (supplement_mode) cts <- strsplit(args[1], "|", fixed = TRUE)[[1]]
lv <- levels(factor(as.character(md$samples)))
stage_of <- function(s) ifelse(grepl("^Healthy", s), "Healthy",
                        ifelse(grepl("^Early", s), "Early", "Advanced"))
stage <- factor(stage_of(lv), levels = c("Healthy", "Early", "Advanced"))
cat("mode:", if (supplement_mode) paste("supplement:", paste(cts, collapse = ";")) else "full",
    " celltypes:", length(cts), "\n")

gmt_lines <- readLines(gmt_file)
hall <- setNames(lapply(strsplit(gmt_lines, "\t"), function(v) v[-(1:2)]),
                 sapply(strsplit(gmt_lines, "\t"), `[`, 1))
ens <- sub("\\..*$", "", rownames(cnt))
sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

cellcounts <- matrix(0, nrow = length(cts), ncol = 12, dimnames = list(cts, lv))
de_merged <- list()
gsea_all <- list()

for (ct in cts) {
  cells <- which(Idents(x) == ct)
  pat <- factor(as.character(md$samples[cells]), levels = lv)
  cellcounts[ct, ] <- as.integer(table(pat))
  if (length(cells) < 100) { cat(ct, ": too few cells, skip\n"); next }
  ncell_pat <- table(pat)
  use <- names(ncell_pat)[ncell_pat >= 10]
  if (length(use) < 6) { cat(ct, ": <6 usable patients (n>=30), skip\n"); next }
  st_use <- stage[match(use, lv)]
  if (any(table(st_use) < 2)) { cat(ct, ": a stage has <2 patients, skip\n"); next }
  ind <- sparseMatrix(i = seq_along(pat), j = as.integer(pat), x = 1, dims = c(length(pat), 12))
  pb <- as.matrix(cnt[, cells, drop = FALSE] %*% ind)
  colnames(pb) <- lv
  pb <- pb[, use, drop = FALSE]
  keep <- rowSums(pb >= 5) >= 3
  dds <- DESeqDataSetFromMatrix(countData = pb[keep, , drop = FALSE],
                                colData = data.frame(stage = st_use), design = ~ stage)
  ok_des <- tryCatch({
    dds <- estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq(dds, quiet = TRUE)
    TRUE
  }, error = function(e) { cat(ct, "DESeq error:", conditionMessage(e), "\n"); FALSE })
  if (!ok_des) next
  r1 <- results(dds, contrast = c("stage", "Advanced", "Healthy"), independentFiltering = FALSE)
  r1$symbol <- ifelse(is.na(sym[rownames(r1)]), rownames(r1), sym[rownames(r1)])
  rk <- r1[!is.na(r1$symbol) & !is.na(r1$log2FoldChange) & !is.na(r1$stat), ]
  rk <- rk[order(abs(rk$stat), decreasing = TRUE), ]
  rk <- rk[!duplicated(rk$symbol), ]
  ranks <- setNames(rk$log2FoldChange, rk$symbol)
  g <- fgseaMultilevel(hall, ranks, minSize = 10, maxSize = 500)
  g2 <- as.data.frame(g)[, c("pathway", "pval", "padj", "NES")]
  g2$celltype <- ct
  gsea_all[[ct]] <- g2
  sig <- r1[which(r1$padj < 0.05 & abs(r1$log2FoldChange) >= 1),
            c("symbol", "log2FoldChange", "padj", "stat")]
  if (nrow(sig) > 0) { sig$celltype <- ct; de_merged[[ct]] <- sig }
  cat(sprintf("%-32s cells=%6d usable=%2d sigDE=%5d gsea_sig=%3d\n", ct, length(cells),
              length(use), nrow(sig), sum(g2$padj < 0.05, na.rm = TRUE)))
}

de_path <- file.path(outdir, "EMTAB_celltype_DE_Advanced_vs_Healthy.csv")
gs_path <- file.path(outdir, "EMTAB_celltype_GSEA_Hallmark.csv")
cc_path <- file.path(outdir, "EMTAB_celltype_pb_cellcounts.csv")

de_df <- if (length(de_merged) > 0) do.call(rbind, de_merged) else NULL
gs <- if (length(gsea_all) > 0) do.call(rbind, gsea_all) else NULL

if (supplement_mode) {
  write.csv(as.data.frame(cellcounts), file.path(outdir, "EMTAB_celltype_pb_cellcounts_supp.csv"))
  if (!is.null(de_df)) write.csv(de_df, file.path(outdir, "EMTAB_celltype_DE_supp.csv"), row.names = FALSE)
  if (!is.null(gs)) write.csv(gs, file.path(outdir, "EMTAB_celltype_GSEA_supp.csv"), row.names = FALSE)
} else {
  write.csv(as.data.frame(cellcounts), cc_path)
  if (!is.null(de_df)) write.csv(de_df, de_path, row.names = FALSE)
  if (!is.null(gs)) write.csv(gs, gs_path, row.names = FALSE)
}

if (!supplement_mode && !is.null(gs)) {
  mat <- reshape(gs[, c("pathway", "celltype", "NES")], idvar = "pathway",
                 timevar = "celltype", direction = "wide")
  rownames(mat) <- mat$pathway; mat$pathway <- NULL
  colnames(mat) <- sub("^NES\\.", "", colnames(mat))
  mat[is.na(mat)] <- 0
  mat <- mat[order(-rowSums(abs(mat))), , drop = FALSE]
  mat <- mat[1:min(25, nrow(mat)), , drop = FALSE]
  png(file.path(outdir, "EMTAB_celltype_GSEA_heatmap.png"), width = 2000, height = 1400, res = 150)
  pheatmap(mat, cluster_cols = TRUE, cluster_rows = TRUE, angle_col = 45,
           color = colorRampPalette(c("#4575B4", "#FFFFBF", "#D73027"))(100),
           fontsize = 9, fontsize_col = 8,
           main = "Hallmark GSEA per cell type (Advanced vs Healthy, E-MTAB-14896)")
  dev.off()
  for (ct in c("Endothelial cells 3", "Fibroblast & Stromal cells 1", "Blood & Immune cells 1")) {
    gct <- gs[gs$celltype == ct, ]
    gct <- gct[order(-abs(gct$NES)), ][1:min(12, nrow(gct)), ]
    gct$pathway <- sub("^HALLMARK_", "", gct$pathway)
    gct$dir <- factor(ifelse(gct$NES > 0, "Up in Advanced", "Down in Advanced"),
                      levels = c("Up in Advanced", "Down in Advanced"))
    gct$pathway <- factor(gct$pathway, levels = rev(gct$pathway[order(gct$NES)]))
    pp <- ggplot(gct, aes(x = NES, y = pathway, fill = dir)) +
      geom_col(width = 0.65) +
      scale_fill_manual(values = c("Up in Advanced" = "#D73027", "Down in Advanced" = "#4575B4")) +
      labs(title = paste("Hallmark GSEA:", ct), x = "NES", y = NULL) +
      theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.title = element_blank())
    fn <- file.path(outdir, paste0("EMTAB_GSEA_", gsub("[ &]", "_", ct), ".png"))
    ggsave(fn, pp, width = 6.8, height = 4.2, dpi = 200)
  }
}
cat("\nDONE\n")
