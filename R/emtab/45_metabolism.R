suppressPackageStartupMessages({library(Seurat); library(Matrix); library(org.Hs.eg.db); library(AnnotationDbi); library(ggplot2)})
outdir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"
gmt_file <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/deep/h.all.v2024.1.Hs.symbols.gmt"

gmt_lines <- readLines(gmt_file)
hall <- setNames(lapply(strsplit(gmt_lines, "\t"), function(v) v[-(1:2)]),
                 sapply(strsplit(gmt_lines, "\t"), `[`, 1))
met_sets <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_GLYCOLYSIS",
              "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_ADIPOGENESIS",
              "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_HEME_METABOLISM",
              "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_XENOBIOTIC_METABOLISM",
              "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_FATTY_ACID_METABOLISM")
met_sets <- unique(met_sets)
hall <- hall[met_sets]

# ---- E-MTAB-14896: patient-level pseudobulk per cell type ----
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
cts <- levels(Idents(x))

score_set <- function(mat_logcpm, genes) {
  g <- genes[genes %in% rownames(mat_logcpm)]
  if (length(g) < 5) return(rep(NA, ncol(mat_logcpm)))
  colMeans(mat_logcpm[g, , drop = FALSE])
}

out_rows <- list()
for (cc in cts) {
  cells <- which(Idents(x) == cc)
  if (length(cells) < 200) next
  pat <- factor(as.character(md$samples[cells]), levels = lv)
  ind <- sparseMatrix(i = seq_along(pat), j = as.integer(pat), x = 1, dims = c(length(pat), 12))
  pb <- as.matrix(cnt2[, cells, drop = FALSE] %*% ind)
  lcpm <- log2(sweep(pb, 2, colSums(pb), "/") * 1e6 + 1)
  for (ms in names(hall)) {
    sc <- score_set(lcpm, hall[[ms]])
    out_rows[[paste(cc, ms)]] <- data.frame(celltype = cc, set = sub("^HALLMARK_", "", ms),
                                            sample = lv, stage = stage, score = sc)
  }
}
met_pb <- do.call(rbind, out_rows)
write.csv(met_pb, file.path(outdir, "EMTAB_metabolism_scores.csv"), row.names = FALSE)

agg <- aggregate(score ~ celltype + set + stage, met_pb, mean)
cat("== metabolism scores by cell type and stage (E-MTAB) ==\n")
print(reshape(agg, idvar = c("celltype", "set"), timevar = "stage", direction = "wide")[order(agg$celltype), ],
      row.names = FALSE)

krus <- do.call(rbind, lapply(split(met_pb, paste(met_pb$celltype, met_pb$set)), function(d) {
  data.frame(celltype = d$celltype[1], set = d$set[1],
             p = tryCatch(kruskal.test(score ~ stage, d)$p.value, error = function(e) NA))
}))
write.csv(krus, file.path(outdir, "EMTAB_metabolism_kruskal.csv"), row.names = FALSE)
cat("== Kruskal by cell type/set (top 12 significant) ==\n")
print(head(krus[order(krus$p), ], 12), row.names = FALSE)

agg2 <- aggregate(score ~ celltype + set, met_pb, mean, na.rm = TRUE)
agg2 <- agg2[complete.cases(agg2), ]
pp <- ggplot(agg2, aes(x = reorder(set, score, FUN = mean), y = score)) +
  geom_boxplot(fill = "#9ECAE1", outlier.shape = NA) +
  coord_flip() + labs(x = NULL, y = "mean score across cell types") +
  theme_bw(base_size = 9)
ggsave(file.path(outdir, "EMTAB_metabolism_summary.png"), pp, width = 6, height = 4, dpi = 160)

cat("\nDONE\n")
