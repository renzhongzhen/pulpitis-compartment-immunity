suppressMessages({
  library(GSVA)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

data_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/data"
res_dir  <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"

# 优先用联合校正矩阵；没有就用 GSE77459 单队列
if (file.exists(file.path(data_dir, "combined.RData"))) {
  load(file.path(data_dir, "combined.RData"))   # E_c, group, batch
  E_mat <- as.matrix(E_c)
  group <- factor(group, levels = c("Normal", "Pulpitis"))
  cat("Using combined ComBat matrix\n")
} else {
  load(file.path(data_dir, "GSE77459_genelevel.RData"))  # expr_gene, grp
  E_mat <- as.matrix(expr_gene)
  group <- factor(grp, levels = c("Normal", "Pulpitis"))
  cat("Using GSE77459 gene-level matrix (no combined yet)\n")
}

# ---------- 免疫细胞标志基因集（文献经典面板，M1/M2 为核心） ----------
marks <- list(
  M1_macrophage = c("IL1B","TNF","IL6","IL12B","CXCL9","CXCL10","NOS2","CD80","CD86","PTGS2"),
  M2_macrophage = c("CD163","MRC1","MSR1","IL10","TGFB1","ARG1","CCL18","CCL22","TGM2"),
  Monocyte      = c("CD14","LYZ","CSF1R","FCGR3A","S100A12","FPR1"),
  Neutrophil    = c("FCGR3B","CSF3R","S100A8","S100A9","CEACAM8","CXCR2"),
  Mast_cell     = c("TPSAB1","TPSB2","CPA3","MS4A2","HDC"),
  NK_cell       = c("NKG7","KLRD1","KLRB1","GNLY","NCAM1","NCR1"),
  CD8_Tcell     = c("CD8A","CD8B","GZMA","GZMB","PRF1"),
  CD4_Tcell     = c("CD4","IL7R","CD40LG","CD3E"),
  Treg          = c("FOXP3","CTLA4","IL2RA","TNFRSF18","CCR8"),
  B_cell        = c("CD19","CD79A","MS4A1","CD22","BLK"),
  Plasma_cell   = c("SDC1","TNFRSF17","SLAMF7","XBP1","MZB1"),
  Dendritic_cell= c("CD1C","ITGAX","FCER1A","CLEC10A","LAMP3")
)

gsets <- lapply(marks, function(g) intersect(g, rownames(E_mat)))
gsets <- gsets[lengths(gsets) >= 3]
cat("Gene sets used:", paste(names(gsets), collapse=" | "), "\n")
# ---------- ssGSEA ----------
param <- ssgseaParam(E_mat, gsets)
ss <- gsva(param)
scores <- as.data.frame(t(ss))
scores$sample <- rownames(scores)
scores$group <- group[rownames(scores)]
write.csv(scores, file.path(res_dir, "immune_ssgsea_scores.csv"), row.names = FALSE)

# ---------- 组间比较（Wilcoxon + limma P 值） ----------
long <- pivot_longer(scores, cols = -c(sample, group), names_to = "cell", values_to = "score")
test_df <- long %>%
  group_by(cell) %>%
  summarise(p_wilcox = wilcox.test(score[group == "Pulpitis"], score[group == "Normal"])$p.value,
            mean_pulpitis = mean(score[group == "Pulpitis"]),
            mean_normal = mean(score[group == "Normal"])) %>%
  mutate(direction = ifelse(mean_pulpitis > mean_normal, "Up", "Down"),
         p_fdr = p.adjust(p_wilcox, "BH")) %>%
  arrange(p_wilcox)
write.csv(as.data.frame(test_df), file.path(res_dir, "immune_cell_comparison.csv"), row.names = FALSE)
print(as.data.frame(test_df))

# ---------- 箱线图 ----------
p_box <- ggplot(long, aes(x = cell, y = score, fill = group)) +
  geom_boxplot(outlier.size = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Immune cell infiltration (ssGSEA)", y = "ssGSEA score", x = "") +
  scale_fill_manual(values = c("Normal" = "#4575b4", "Pulpitis" = "#d73027"))
ggsave(file.path(res_dir, "immune_boxplot.png"), p_box, width = 9, height = 5.5, dpi = 300)

# ---------- 热图 ----------
library(pheatmap)
ann_col <- data.frame(Group = group, row.names = rownames(scores))
mat <- t(as.matrix(scores[, setdiff(colnames(scores), c("sample","group"))]))
p_hm <- pheatmap(mat, annotation_col = ann_col, scale = "row",
                 annotation_colors = list(Group = c(Normal = "#4575b4", Pulpitis = "#d73027")),
                 main = "Immune cell ssGSEA scores", show_colnames = FALSE,
                 fontsize = 9, silent = TRUE)
ggsave(file.path(res_dir, "immune_heatmap.png"), p_hm, width = 9, height = 5.5, dpi = 300)

cat("\nDONE\n")
