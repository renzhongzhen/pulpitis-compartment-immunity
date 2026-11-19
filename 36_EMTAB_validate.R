suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2)
})
p <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/ea960c355a014573b0e904afc349b633processedmatrix.rds"
outdir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
x <- readRDS(p)
grp <- factor(x$samples)
lv <- levels(grp)
cat("samples:", length(lv), "\n"); print(lv)

# --- pseudobulk counts ---
cnt <- LayerData(x, assay = "RNA", layer = "counts")
ind <- sparseMatrix(i = seq_along(grp), j = as.integer(grp), x = 1, dims = c(length(grp), length(lv)))
pb <- as.matrix(cnt %*% ind)
colnames(pb) <- lv

ens <- sub("\\..*$", "", rownames(pb))
sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
gn <- ifelse(is.na(sym), ens, sym)
rownames(pb) <- make.unique(gn)

cpm <- sweep(pb, 2, colSums(pb), "/") * 1e6
lcpm <- log2(cpm + 1)

stage_of <- function(s) ifelse(grepl("^Healthy", s), "Healthy",
                        ifelse(grepl("^Early", s), "Early", "Advanced"))
stage <- factor(stage_of(lv), levels = c("Healthy", "Early", "Advanced"))
cat("\n== pseudobulk size ==\n"); print(dim(pb))

# --- 7-gene signature ---
g7 <- c("ADGRL4", "CCL21", "FOSB", "ICAM2", "KIT", "PI3", "SELE")
g7 <- g7[g7 %in% rownames(lcpm)]
cat("signature genes found:", paste(g7, collapse = ","), "\n")
sig <- as.data.frame(t(lcpm[g7, , drop = FALSE]))
sig$sample <- rownames(sig); sig$stage <- stage[match(rownames(sig), lv)]
write.csv(sig, file.path(outdir, "EMTAB_signature7_by_sample.csv"), row.names = FALSE)

# --- M1/M2 module scores (pseudobulk mean of markers) ---
m1 <- c("IL1B", "CXCL10", "PTGS2", "IL6", "CXCL9", "CD86", "NOS2", "TNF", "CCL3", "CCL4")
m2 <- c("CD163", "MRC1", "MSR1", "IL10", "ARG1", "CCL18", "TGFB1", "CD209")
m1 <- m1[m1 %in% rownames(lcpm)]; m2 <- m2[m2 %in% rownames(lcpm)]
cat("M1 markers:", paste(m1, collapse = ","), "\nM2 markers:", paste(m2, collapse = ","), "\n")
mod <- data.frame(sample = lv, stage = stage,
                  M1 = colMeans(lcpm[m1, , drop = FALSE]),
                  M2 = colMeans(lcpm[m2, , drop = FALSE]))
mod$M1_minus_M2 <- mod$M1 - mod$M2
write.csv(mod, file.path(outdir, "EMTAB_module_scores.csv"), row.names = FALSE)
cat("\n== module scores by stage ==\n"); print(aggregate(cbind(M1, M2, M1_minus_M2) ~ stage, mod, function(v) round(mean(v), 3)))
cat("Kruskal M1 p=", kruskal.test(M1 ~ stage, mod)$p.value,
    " M2 p=", kruskal.test(M2 ~ stage, mod)$p.value,
    " M1-M2 p=", kruskal.test(M1_minus_M2 ~ stage, mod)$p.value, "\n", sep = "")

# --- immune / endothelial / stromal composition per sample ---
comp <- as.data.frame(prop.table(table(x$samples, Idents(x)), 1) * 100)
colnames(comp) <- c("sample", "celltype", "pct")
comp$stage <- stage[match(comp$sample, lv)]
write.csv(comp, file.path(outdir, "EMTAB_composition_by_sample.csv"), row.names = FALSE)
agg <- aggregate(pct ~ celltype + stage, comp, function(v) round(mean(v), 2))
cat("\n== mean pct by celltype x stage ==\n"); print(reshape(agg, idvar = "celltype", timevar = "stage", direction = "wide"))

# trend test per celltype (Kruskal over 3 stages)
tr <- do.call(rbind, lapply(split(comp, comp$celltype), function(d) {
  k <- tryCatch(kruskal.test(pct ~ stage, d)$p.value, error = function(e) NA)
  data.frame(celltype = d$celltype[1], p_kruskal = signif(k, 3),
             Healthy = round(mean(d$pct[d$stage == "Healthy"]), 2),
             Early = round(mean(d$pct[d$stage == "Early"]), 2),
             Advanced = round(mean(d$pct[d$stage == "Advanced"]), 2))
}))
write.csv(tr, file.path(outdir, "EMTAB_composition_trend.csv"), row.names = FALSE)
cat("\n== composition trend ==\n"); print(tr[order(tr$p_kruskal), ])

# --- key genes by celltype (average lognorm) ---
dat <- LayerData(x, assay = "RNA", layer = "data")
kg <- unique(c(g7, m1, m2, "CXCL2", "CXCL8", "CCL2", "ACKR1", "SELE", "ICAM1", "ALPL", "RUNX2", "VEGFB", "FGFR2", "VDR", "MMP8"))
key <- kg[kg %in% rownames(dat)]
ct <- Idents(x)
avg <- sapply(levels(ct), function(cc) Matrix::rowMeans(dat[key, ct == cc, drop = FALSE]))
write.csv(round(as.data.frame(avg), 3), file.path(outdir, "EMTAB_keygene_celltype_avg.csv"))

# --- plots ---
stage_cols <- c(Healthy = "#4C72B0", Early = "#DD8452", Advanced = "#C44E52")

p1 <- ggplot(agg, aes(x = stage, y = pct, fill = celltype)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_manual(values = rep(c("#4C72B0", "#55A868", "#C44E52", "#8172B2", "#937860",
                                   "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD", "#B07AA1",
                                   "#76B7B2", "#FF9DA7"), 2)) +
  labs(x = NULL, y = "Mean proportion (%)", title = "E-MTAB-14896: cell composition across disease stages") +
  theme_bw(base_size = 10) + theme(legend.text = element_text(size = 7))
ggsave(file.path(outdir, "EMTAB_composition_stage.png"), p1, width = 7.2, height = 4.6, dpi = 200)

hs <- scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027", midpoint = 3)
ph <- as.data.frame(as.table(lcpm[g7, , drop = FALSE]))
colnames(ph) <- c("gene", "sample", "expr")
ph$stage <- stage[match(ph$sample, lv)]
p2 <- ggplot(ph, aes(x = sample, y = gene, fill = expr)) + geom_tile(colour = "white", linewidth = 0.4) +
  hs + facet_grid(~stage, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL, title = "Seven-gene signature (pseudobulk log2 CPM)", fill = "log2 CPM") +
  theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))
ggsave(file.path(outdir, "EMTAB_signature7_heatmap.png"), p2, width = 7.6, height = 3.4, dpi = 200)

md <- rbind(data.frame(stage = mod$stage, module = "M1", score = mod$M1),
            data.frame(stage = mod$stage, module = "M2", score = mod$M2))
p3 <- ggplot(md, aes(x = stage, y = score, fill = stage)) + geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 1.3, alpha = 0.8) + facet_wrap(~module, scales = "free_y") +
  scale_fill_manual(values = stage_cols) +
  labs(x = NULL, y = "Mean marker log2 CPM", title = "Macrophage polarisation modules across stages") +
  theme_bw(base_size = 10) + theme(legend.position = "none")
ggsave(file.path(outdir, "EMTAB_M1M2_stage.png"), p3, width = 6.2, height = 3.6, dpi = 200)

cat("\nDONE\n")
