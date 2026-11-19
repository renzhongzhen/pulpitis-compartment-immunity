# -*- coding: utf-8 -*-
# 22_composition_genes_GSE185222.R : annotation, composition, key gene localization, M1/M2 score
suppressMessages({
  library(Seurat); library(dplyr); library(ggplot2); library(tidyr)
})

out_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"
obj <- readRDS(file.path(out_dir, "GSE185222_seurat.rds"))

anno <- c("0" = "Odontoblast-like (OMD+)", "1" = "Fibroblasts", "2" = "Fibroblasts/MSC",
          "3" = "Endothelial", "4" = "T/NK cells", "5" = "Macrophages/Monocytes",
          "6" = "Pericytes/SMC", "7" = "B cells", "8" = "Odontoblasts (DMP1+/IBSP+)",
          "9" = "Neural/Schwann (SPP1+)", "10" = "Plasma cells", "11" = "Schwann cells (MPZ+)",
          "12" = "Mixed/low-quality")
obj$celltype <- unname(anno[as.character(obj$seurat_clusters)])
obj$condition <- ifelse(obj$sample %in% c("DTP_01", "DTP_04"), "Healthy/Enamel caries", "Deep caries")

# ---- composition ----
comp <- as.data.frame(table(celltype = obj$celltype, condition = obj$condition))
comp$pct <- ave(comp$Freq, comp$condition, FUN = function(x) 100 * x / sum(x))
write.csv(comp, file.path(out_dir, "GSE185222_celltype_composition.csv"), row.names = FALSE)
print(comp %>% group_by(condition) %>% summarise(total = sum(Freq)))

df <- as.data.frame(table(sample = obj$sample, celltype = obj$celltype, condition = obj$condition))
df$pct <- ave(df$Freq, df$sample, FUN = function(x) 100 * x / sum(x))
p_bar <- ggplot(df, aes(x = sample, y = pct, fill = celltype)) +
  geom_col(width = 0.7) + facet_wrap(~condition, scales = "free_x") +
  labs(y = "Percent of cells (%)", x = "", title = "Cell type composition per sample") +
  theme_bw() + theme(legend.text = element_text(size = 7))
ggsave(file.path(out_dir, "GSE185222_composition_bar.png"), p_bar, width = 9, height = 6, dpi = 200)

# ---- key gene dotplot ----
genes <- list(
  M1_markers = c("IL1B","CXCL10","PTGS2","IL6","CXCL9","CD86"),
  M2_markers = c("CD163","MRC1","MSR1","IL10","ARG1","CCL18"),
  ROC_up     = c("PI3","CCL21","SELE","ACP5","SOCS3","CH25H","ATP6V0D2","SLC37A2","FOSB"),
  ROC_down   = c("NR1D1","GAP43","CRISPLD1","SPON1","MOXD1","DBP","NR1D2","RSPO2","ECM2","RN7SK"),
  PBM_genes  = c("VDR","MMP8","MMP10","BMP4","VEGFB","ALPL","TWIST1","RUNX2","TNFSF11","CD36"),
  Cross_genes= c("SPP1","MMP9","ITGA2","BMP5","COL10A1","ICAM1","FGFR2")
)
for (g in names(genes)) genes[[g]] <- genes[[g]][genes[[g]] %in% rownames(obj)]

p_dot <- DotPlot(obj, features = genes, group.by = "celltype", assay = "RNA") +
  RotatedAxis() + labs(title = "Key genes across cell types (GSE185222)", x = "", y = "")
ggsave(file.path(out_dir, "GSE185222_keygene_dotplot.png"), p_dot, width = 22, height = 7, dpi = 200)

# average expression per celltype -> csv
avg <- AverageExpression(obj, features = unique(unlist(genes)), group.by = "celltype",
                         assays = "RNA", layer = "data")$RNA
write.csv(as.matrix(avg), file.path(out_dir, "GSE185222_keygene_celltype_avg.csv"))

# ---- macrophage M1/M2 score ----
obj_mac <- subset(obj, celltype == "Macrophages/Monocytes")
obj_mac <- AddModuleScore(obj_mac, features = list(
  M1 = c("IL1B","CXCL10","PTGS2","IL6","CXCL9","CD86"),
  M2 = c("CD163","MRC1","MSR1","IL10","ARG1","CCL18")), name = "M1M2", seed = 1)
md <- obj_mac@meta.data[, c("sample","condition","M1M21","M1M22")]
colnames(md) <- c("sample","condition","M1_score","M2_score")
write.csv(md, file.path(out_dir, "GSE185222_mac_M1M2_scores.csv"), row.names = FALSE)
cat("\nM1/M2 score by condition:\n")
print(md %>% group_by(condition) %>% summarise(M1 = mean(M1_score), M2 = mean(M2_score), n = n()))

pm <- ggplot(md, aes(x = condition, y = M1_score, fill = condition)) + geom_boxplot() +
  labs(title = "Macrophage M1 module score", y = "M1 score", x = "") + theme_bw() +
  theme(legend.position = "none")
pm2 <- ggplot(md, aes(x = condition, y = M2_score, fill = condition)) + geom_boxplot() +
  labs(title = "Macrophage M2 module score", y = "M2 score", x = "") + theme_bw() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "GSE185222_mac_M1M2_boxplot.png"), pm + pm2, width = 8, height = 5, dpi = 200)

cat("\nCelltype x condition table:\n"); print(table(obj$celltype, obj$condition))
saveRDS(obj, file.path(out_dir, "GSE185222_seurat_annotated.rds"))
cat("DONE\n")
