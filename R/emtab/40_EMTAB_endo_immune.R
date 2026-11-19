suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(org.Hs.eg.db); library(AnnotationDbi)
  library(ggplot2)
})
p <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/ea960c355a014573b0e904afc349b633processedmatrix.rds"
outdir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"

x <- readRDS(p)
md <- x@meta.data
cnt <- LayerData(x, assay = "RNA", layer = "counts")
dat <- LayerData(x, assay = "RNA", layer = "data")
ens <- sub("\\..*$", "", rownames(dat))
sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
keep <- !is.na(sym) & !duplicated(sym)
dat2 <- dat[keep, , drop = FALSE]; rownames(dat2) <- sym[keep]
cnt2 <- cnt[keep, , drop = FALSE]; rownames(cnt2) <- sym[keep]

lv <- levels(factor(as.character(md$samples)))
stage_of <- function(s) ifelse(grepl("^Healthy", s), "Healthy",
                        ifelse(grepl("^Early", s), "Early", "Advanced"))
stage_pat <- factor(stage_of(lv), levels = c("Healthy", "Early", "Advanced"))
stage_cell <- stage_pat[factor(as.character(md$samples), levels = lv)]

# ================= 2. endothelial subtype annotation =================
endo_cts <- paste0("Endothelial cells", c(" 1", " 2", " 3", " 4"))
if (!all(endo_cts %in% levels(Idents(x)))) endo_cts <- levels(Idents(x))[grepl("ndothelial", levels(Idents(x)))]
e_markers <- c("ACKR1", "SELE", "VCAM1", "ICAM1", "CCL2", "CXCL8",
               "CA4", "RGCC", "VWA1", "PLVAP",
               "GJA5", "HEY1", "SEMA3G", "EFNB2", "FBLN5",
               "ESM1", "KDR", "DLL4", "FLT1", "PGF",
               "PROX1", "LYVE1", "PDPN", "CCL21")
e_markers <- e_markers[e_markers %in% rownames(dat2)]
em <- sapply(endo_cts, function(cc) Matrix::rowMeans(dat2[e_markers, Idents(x) == cc, drop = FALSE]))
write.csv(round(em, 3), file.path(outdir, "EMTAB_endothelial_subtype_markers.csv"))

eml <- as.data.frame(as.table(em)); colnames(eml) <- c("marker", "cluster", "expr")
eml$expr <- as.numeric(eml$expr)
eml$marker <- factor(eml$marker, levels = rev(e_markers))
eml$cluster <- factor(eml$cluster, levels = endo_cts)
pp <- ggplot(eml, aes(x = cluster, y = marker, fill = expr)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#4575B4", mid = "#FFFFBF", high = "#D73027", midpoint = median(eml$expr)) +
  labs(x = NULL, y = NULL, title = "Endothelial subtype markers (E-MTAB-14896)", fill = "mean lognorm") +
  theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(outdir, "EMTAB_endothelial_marker_heatmap.png"), pp, width = 5.6, height = 6.4, dpi = 200)

endo_prop <- lapply(endo_cts, function(cc) {
  d <- data.frame(sample = lv, stage = stage_pat, n = as.integer(table(factor(md$samples[Idents(x) == cc], levels = lv))))
  d$total <- as.integer(table(md$samples))[match(d$sample, names(table(md$samples)))]
  d$pct <- 100 * d$n / d$total
  d$cluster <- cc
  d
})
ep <- do.call(rbind, endo_prop)
agg_ep <- aggregate(pct ~ cluster + stage, ep, mean)
cat("== endothelial subtype pct by stage ==\n")
print(reshape(agg_ep, idvar = "cluster", timevar = "stage", direction = "wide"))
write.csv(ep, file.path(outdir, "EMTAB_endothelial_subtype_pct.csv"), row.names = FALSE)

# ================= 3. immune M1/M2 scoring =================
imm_cts <- levels(Idents(x))[grepl("Blood & Immune", levels(Idents(x)))]
m1 <- c("IL1B", "CXCL10", "PTGS2", "IL6", "CXCL9", "CD86", "NOS2", "TNF", "CCL3", "CCL4")
m2 <- c("CD163", "MRC1", "MSR1", "IL10", "ARG1", "CCL18", "TGFB1", "CD209")
m1 <- m1[m1 %in% rownames(dat2)]; m2 <- m2[m2 %in% rownames(dat2)]
imm_cells <- Idents(x) %in% imm_cts
sub <- dat2[c(m1, m2), imm_cells, drop = FALSE]
subm <- as.matrix(sub)
z <- t(scale(t(subm)))
score <- data.frame(cell = colnames(subm), celltype = as.character(Idents(x)[imm_cells]),
                    stage = stage_cell[imm_cells], sample = as.character(md$samples[imm_cells]))
score$M1 <- colMeans(z[m1, , drop = FALSE])
score$M2 <- colMeans(z[m2, , drop = FALSE])
score$M1_minus_M2 <- score$M1 - score$M2
write.csv(score, file.path(outdir, "EMTAB_immune_M1M2_cellscore.csv"), row.names = FALSE)
cat("== cell-level M1/M2 by stage (mean) ==\n")
print(aggregate(cbind(M1, M2, M1_minus_M2) ~ stage, score, mean))
cat("cell-level Wilcoxon M1-M2 Advanced vs Healthy p=",
    wilcox.test(score$M1_minus_M2[score$stage == "Advanced"],
                score$M1_minus_M2[score$stage == "Healthy"])$p.value, "\n")
cat("cell-level Wilcoxon M1 Advanced vs Healthy p=",
    wilcox.test(score$M1[score$stage == "Advanced"],
                score$M1[score$stage == "Healthy"])$p.value, "\n")

# patient-level pseudobulk within immune compartment
cells <- which(imm_cells)
pat <- factor(as.character(md$samples[cells]), levels = lv)
ind <- sparseMatrix(i = seq_along(pat), j = as.integer(pat), x = 1, dims = c(length(pat), 12))
pbi <- as.matrix(cnt2[, cells, drop = FALSE] %*% ind)
cpm_i <- sweep(pbi, 2, colSums(pbi), "/") * 1e6
lcpm_i <- log2(cpm_i + 1)
imm_mod <- data.frame(sample = lv, stage = stage_pat,
                      M1 = colMeans(lcpm_i[m1, , drop = FALSE]),
                      M2 = colMeans(lcpm_i[m2, , drop = FALSE]))
imm_mod$M1_minus_M2 <- imm_mod$M1 - imm_mod$M2
write.csv(imm_mod, file.path(outdir, "EMTAB_immune_M1M2_by_patient.csv"), row.names = FALSE)
cat("== patient-level immune M1/M2 by stage ==\n")
print(aggregate(cbind(M1, M2, M1_minus_M2) ~ stage, imm_mod, mean))
cat("patient-level Kruskal M1-M2 p=", kruskal.test(M1_minus_M2 ~ stage, imm_mod)$p.value, "\n")

mdl <- rbind(data.frame(stage = score$stage, module = "M1", value = score$M1),
             data.frame(stage = score$stage, module = "M2", value = score$M2))
pp2 <- ggplot(mdl, aes(x = stage, y = value, fill = stage)) +
  geom_boxplot(outlier.shape = NA, width = 0.55) + facet_wrap(~module, scales = "free_y") +
  scale_fill_manual(values = c(Healthy = "#4C72B0", Early = "#DD8452", Advanced = "#C44E52")) +
  labs(x = NULL, y = "z-scored module score", title = "M1/M2 module scores within immune cells (E-MTAB-14896)") +
  theme_bw(base_size = 10) + theme(legend.position = "none")
ggsave(file.path(outdir, "EMTAB_immune_M1M2_stage.png"), pp2, width = 6.2, height = 3.6, dpi = 200)

# ================= 4b. downsampling robustness of composition =================
set.seed(20261119)
nsim <- 20; n_sub <- 1000
cts_all <- levels(Idents(x))
sub_res <- list()
for (i in seq_len(nsim)) {
  keep_cells <- unlist(lapply(lv, function(s) {
    w <- which(md$samples == s)
    if (length(w) <= n_sub) w else sample(w, n_sub)
  }))
  tab <- table(factor(Idents(x)[keep_cells], levels = cts_all),
               factor(as.character(md$samples[keep_cells]), levels = lv))
  pct <- 100 * prop.table(tab, 2)
  grp <- stage_of(colnames(pct))
  d <- rowMeans(pct[, grp == "Advanced"]) - rowMeans(pct[, grp == "Healthy"])
  sub_res[[i]] <- d
}
submat <- do.call(rbind, sub_res)
colnames(submat) <- cts_all
subdf <- data.frame(celltype = cts_all,
                    mean_diff = colMeans(submat),
                    lo = apply(submat, 2, quantile, 0.025),
                    hi = apply(submat, 2, quantile, 0.975),
                    sign_consistent = 100 * colMeans(sign(submat) == sign(colMeans(submat))))
subdf <- subdf[order(-subdf$mean_diff), ]
write.csv(subdf, file.path(outdir, "EMTAB_downsampling_robustness.csv"), row.names = FALSE)
cat("== downsampling robustness (mean Advanced-Healthy pct diff, 95% CI) ==\n")
print(subdf, row.names = FALSE)

cat("\nDONE\n")
