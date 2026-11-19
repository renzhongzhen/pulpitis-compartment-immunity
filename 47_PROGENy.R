suppressPackageStartupMessages({
  library(progeny); library(dplyr)
  library(Seurat); library(Matrix); library(org.Hs.eg.db); library(AnnotationDbi)
  library(ggplot2)
})
res_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"
sc_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc"
outdir <- file.path(sc_dir, "results")

# ---- 1. bulk GSE77459 ----
expr <- read.csv(file.path(res_dir, "GSE77459_expr_gene.csv"), row.names = 1, check.names = FALSE)
grp <- read.csv(file.path(res_dir, "GSE77459_groups.csv"))
g <- grp$group[match(colnames(expr), grp$sample)]
pa <- progeny(as.matrix(expr), scale = TRUE, organism = "Human", top = 100, perm = 1)
pa <- as.data.frame(pa)
pa$sample <- rownames(pa)
pa$group <- g[match(pa$sample, grp$sample)]
paths <- setdiff(colnames(pa), c("sample", "group"))
cmp <- do.call(rbind, lapply(paths, function(pth) {
  x <- pa[[pth]][pa$group == "Pulpitis"]; y <- pa[[pth]][pa$group == "Normal"]
  data.frame(pathway = pth, normal = mean(y), pulpitis = mean(x),
             diff = mean(x) - mean(y),
             p = tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA))
}))
cmp$padj <- p.adjust(cmp$p, "BH")
cmp <- cmp[order(-cmp$diff), ]
write.csv(cmp, file.path(res_dir, "deep", "PROGENy_bulk_compare.csv"), row.names = FALSE)
cat("== bulk PROGENy pathway activity (by diff) ==\n")
print(cmp, row.names = FALSE)

cmp$pathway <- factor(cmp$pathway, levels = rev(cmp$pathway[order(cmp$diff)]))
cmp$dir <- ifelse(cmp$diff > 0, "Up in pulpitis", "Down in pulpitis")
pp <- ggplot(cmp, aes(x = diff, y = pathway, fill = dir)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Up in pulpitis" = "#D73027", "Down in pulpitis" = "#4575B4")) +
  labs(title = "PROGENy pathway activity (pulpitis minus normal, GSE77459)", x = "activity difference", y = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(res_dir, "deep", "PROGENy_bulk.png"), pp, width = 6.5, height = 4.5, dpi = 160)

# ---- 2. E-MTAB three stages ----
p <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/ea960c355a014573b0e904afc349b633processedmatrix.rds"
x <- readRDS(p)
md <- x@meta.data
cnt <- LayerData(x, assay = "RNA", layer = "counts")
ens <- sub("\\..*$", "", rownames(cnt))
sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
keep <- !is.na(sym) & !duplicated(sym)
cnt2 <- cnt[keep, , drop = FALSE]; rownames(cnt2) <- sym[keep]
lv <- levels(factor(as.character(md$samples)))
stage_of <- function(s) ifelse(grepl("^Healthy", s), "Healthy",
                        ifelse(grepl("^Early", s), "Early", "Advanced"))
stage <- factor(stage_of(lv), levels = c("Healthy", "Early", "Advanced"))
pat <- factor(as.character(md$samples), levels = lv)
ind <- sparseMatrix(i = seq_along(pat), j = as.integer(pat), x = 1, dims = c(length(pat), 12))
pb <- as.matrix(cnt2 %*% ind)
colnames(pb) <- lv

pa2 <- progeny(pb, scale = TRUE, organism = "Human", top = 100, perm = 1)
pa2 <- as.data.frame(pa2)
pa2$stage <- stage[match(rownames(pa2), lv)]
write.csv(pa2, file.path(outdir, "EMTAB_PROGENy_activity.csv"), row.names = FALSE)

tr <- do.call(rbind, lapply(paths, function(pth) {
  v <- pa2[[pth]]
  m <- tapply(v, pa2$stage, mean)
  data.frame(pathway = pth, Healthy = m[1], Early = m[2], Advanced = m[3],
             p = kruskal.test(v ~ pa2$stage)$p.value)
}))
tr$padj <- p.adjust(tr$p, "BH")
tr$mono_up <- tr$Advanced > tr$Early & tr$Early > tr$Healthy
tr <- tr[order(-tr$padj), ]
write.csv(tr, file.path(outdir, "EMTAB_PROGENy_trend.csv"), row.names = FALSE)
cat("== E-MTAB PROGENy by stage ==\n")
print(tr, row.names = FALSE)

focus <- tr$pathway[tr$padj < 0.15 | tr$mono_up][1:min(8, sum(tr$padj < 0.15 | tr$mono_up))]
pdl <- pa2[, c("stage", focus)]
pdl <- reshape(pdl, direction = "long", varying = list(focus), v.names = "act",
               times = focus, timevar = "pathway")
colnames(pdl)[2] <- "pathway"
pp2 <- ggplot(pdl, aes(x = stage, y = act, group = pathway, color = pathway)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  labs(title = "PROGENy pathway activity across stages (E-MTAB-14896)", x = NULL, y = "activity (scaled)") +
  theme_bw(base_size = 10) + theme(legend.position = "right", legend.title = element_blank())
ggsave(file.path(outdir, "EMTAB_PROGENy_stage.png"), pp2, width = 7.5, height = 5, dpi = 160)
cat("\nDONE\n")
