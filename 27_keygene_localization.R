# 27_keygene_localization.R : localization of PBM-panel / biomarker / cross-check genes across cell types
suppressMessages({ library(Seurat); library(ggplot2) })

sc_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"
obj <- readRDS(file.path(sc_dir, "GSE185222_seurat_annotated.rds"))
obj <- subset(obj, celltype != "Mixed/low-quality")

genes <- c(
  # PBM panel core
  "VDR", "MMP8", "MMP10", "BMP4", "VEGFB", "ALPL", "TWIST1", "RUNX2", "TNFSF11", "FGFR2",
  # cross-check
  "SPP1", "MMP9", "ITGA2", "BMP5", "ICAM1", "COL10A1", "CD36",
  # M1/M2 markers
  "IL1B", "CXCL10", "PTGS2", "IL6", "CXCL9", "CD86", "CD163", "MRC1", "MSR1", "IL10", "ARG1",
  # ROC candidates
  "PI3", "CCL21", "SELE", "ACP5", "SOCS3", "CH25H", "ATP6V0D2", "SLC37A2", "FOSB",
  "NR1D1", "GAP43", "SPON1", "DBP", "NR1D2", "RSPO2", "ECM2"
)
present <- genes[genes %in% rownames(obj)]
missing <- setdiff(genes, present)
cat("Present:", length(present), "genes\n")
cat("Missing:", paste(missing, collapse = ", "), "\n")

data_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
ct <- as.character(obj$celltype)
lev <- sort(unique(ct))
avg <- sapply(lev, function(g) Matrix::rowMeans(data_mat[present, ct == g, drop = FALSE]))
colnames(avg) <- lev
avg <- round(avg, 3)
write.csv(avg, file.path(sc_dir, "GSE185222_allkeygene_celltype_avg.csv"))

# top cell type per gene
top_ct <- apply(avg, 1, function(x) names(which.max(x)))
peak <- apply(avg, 1, max)
tab <- data.frame(gene = present, top_celltype = top_ct, max_avg = round(peak, 3))
write.csv(tab, file.path(sc_dir, "GSE185222_keygene_topcelltype.csv"), row.names = FALSE)
print(tab, row.names = FALSE)

# DotPlot
p <- DotPlot(obj, features = present, group.by = "celltype", cols = "RdYlBu") +
  RotatedAxis() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
ggsave(file.path(sc_dir, "GSE185222_allkeygene_dotplot.png"), p, width = 7, height = 9, dpi = 300)
cat("DONE\n")
