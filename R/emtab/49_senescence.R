suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(org.Hs.eg.db); library(AnnotationDbi)
  library(ggplot2)
})
sc_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc"
outdir <- file.path(sc_dir, "results")
# SenMayo senescence gene set (Saul et al. 2022, Nat Commun) - canonical list
senmayo <- c("ACVR1B", "ANG", "ANGLPTL4", "AREG", "AXL", "BEX3", "BMP2", "BMP6", "BRIX1",
  "BTG2", "CAPN2", "CAV1", "CAVIN3", "CBLB", "CCL2", "CCL24", "CCL26", "CCL3", "CCL4", "CCL5",
  "CCL7", "CCL8", "CCL13", "CCL20", "CCL25", "CCL27", "CCL28", "CD55", "CDKN1A", "CDKN1B",
  "CDKN2B", "COL1A2", "COL3A1", "COL4A1", "COL4A2", "COL5A1", "COL5A2", "COL6A1", "COL6A2",
  "COL6A3", "CSF1", "CSF2", "CSF2RB", "CTSB", "CTSL1", "CXCL1", "CXCL2", "CXCL3", "CXCL5",
  "CXCL8", "CXCL12", "CXCL14", "CXCL16", "CXCR2", "DDIT4", "DNAJB4", "EDN1", "EGF", "EGFR",
  "EREG", "ESM1", "ETV5", "FAS", "FGF1", "FGF2", "FGF7", "FGF13", "FLT4", "FOXM1", "GDF15",
  "GEM", "GPC4", "HGF", "HMGB1", "HNRNPD", "HSPA1A", "HSPA1B", "HSPA2", "HSPA8", "ICAM1",
  "ICAM2", "ICAM3", "ICAM4", "ICAM5", "IFNG", "IGF1", "IGF2", "IGFBP1", "IGFBP2", "IGFBP3",
  "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IGF1R", "IL1A", "IL1B", "IL6", "IL6ST", "IL7",
  "IL10", "IL13", "IL15", "IL18", "IL32", "INHBA", "IRF1", "IRF3", "IRF7", "ITGAV", "ITGB1",
  "ITGB3", "ITGB4", "JUN", "KITLG", "LCP1", "MAP2K1", "MAP2K3", "MAP2K6", "MAP3K5", "MAPK14",
  "MIF", "MMP1", "MMP2", "MMP3", "MMP9", "MMP10", "MMP12", "MMP13", "MMP14", "NAMPT", "NFKB1",
  "NFKB2", "NOX4", "NRG1", "PAPPA", "PDGFA", "PDGFB", "PDGFD", "PDGFRA", "PDGFRB", "PECAM1",
  "PLAT", "PLAUR", "PLAU", "PTGER2", "PTGES", "PTGS2", "SERPINB2", "SERPINE1", "SERPINE2",
  "SIRT1", "SIRT3", "SOD1", "SOD2", "SPARC", "TIMP1", "TIMP2", "TIMP3", "TIMP4", "TNF", "TNFRSF1A",
  "TNFRSF1B", "TNFSF10", "TNFRSF10A", "TNFRSF10B", "TNFSF13", "TNFRSF13B", "TNFSF15", "TGFB1",
  "TGFB2", "TGFB3", "TGFBR1", "TGFBR2", "TGFBR3", "THBS1", "THBS2", "TP53", "TP63", "TRAF3IP2",
  "TXN", "VEGFA", "VEGFB", "VEGFC", "VEGFD", "VIM", "WNT16", "WNT2", "YAP1")
senmayo <- unique(senmayo)

p <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/ea960c355a014573b0e904afc349b633processedmatrix.rds"
x <- readRDS(p)
md <- x@meta.data
ct <- factor(as.character(Idents(x)))
cnt <- LayerData(x, assay = "RNA", layer = "counts")
ens <- sub("\\..*$", "", rownames(cnt))
sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
keep <- !is.na(sym) & !duplicated(sym)
cnt2 <- cnt[keep, , drop = FALSE]; rownames(cnt2) <- sym[keep]
cs <- colSums(cnt2); cs[cs == 0] <- 1
genes <- senmayo[senmayo %in% rownames(cnt2)]
cat("SenMayo genes present:", length(genes), "of", length(senmayo), "\n")
Esub <- as.matrix(cnt2[genes, , drop = FALSE])
score <- colMeans(log2(sweep(Esub, 2, cs, "/") * 1e4 + 1))
md$senmayo <- as.numeric(score)

stage_of_pat <- function(s) ifelse(grepl("^Healthy", s), "Healthy",
                            ifelse(grepl("^Early", s), "Early", "Advanced"))
md$stage <- factor(stage_of_pat(as.character(md$samples)), levels = c("Healthy", "Early", "Advanced"))
md$ct <- as.character(ct)
write.csv(md[, c("ct", "stage", "samples", "senmayo")], file.path(outdir, "EMTAB_senmayo_scores.csv"), row.names = FALSE)

agg <- aggregate(senmayo ~ ct + stage, md, mean)
agg$stage <- factor(agg$stage, levels = c("Healthy", "Early", "Advanced"))
cat("== SenMayo score by cell type and stage ==\n")
print(agg[order(agg$ct, agg$stage), ], row.names = FALSE)

krus <- do.call(rbind, lapply(unique(md$ct), function(cc) {
  v <- md[md$ct == cc, ]
  data.frame(ct = cc, p = kruskal.test(senmayo ~ stage, v)$p.value,
             H = mean(v$senmayo[v$stage == "Healthy"]),
             E = mean(v$senmayo[v$stage == "Early"]),
             A = mean(v$senmayo[v$stage == "Advanced"]))
}))
krus$padj <- p.adjust(krus$p, "BH")
write.csv(krus, file.path(outdir, "EMTAB_senmayo_stats.csv"), row.names = FALSE)
cat("== SenMayo Kruskal by cell type ==\n")
print(krus[order(krus$p), ], row.names = FALSE)

sub_ct <- krus$ct[krus$padj < 0.2 | krus$ct %in% c("Endothelial cells 3", "Fibroblast & Stromal cells 1")]
if (length(sub_ct) < 3) sub_ct <- head(krus$ct[order(krus$p)], 6)
agg2 <- agg[agg$ct %in% sub_ct, ]
pp <- ggplot(agg2, aes(x = stage, y = senmayo, group = ct, color = ct)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
  labs(title = "SenMayo senescence score across stages (E-MTAB-14896)", x = NULL, y = "SenMayo score") +
  theme_bw(base_size = 10) + theme(legend.position = "right", legend.title = element_blank())
ggsave(file.path(outdir, "EMTAB_senmayo_stage.png"), pp, width = 8, height = 5, dpi = 160)
cat("\nDONE\n")
