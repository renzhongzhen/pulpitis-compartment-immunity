suppressMessages({
  library(dplyr)
  library(pheatmap)
  library(ggplot2)
})

data_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/data"
res_dir  <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"

if (!requireNamespace("pROC", quietly = TRUE)) {
  install.packages("pROC", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
}
library(pROC)

load(file.path(data_dir, "GSE77459_genelevel.RData"))  # expr_gene, grp
deg <- read.csv(file.path(res_dir, "GSE77459_DEG_full.csv"))
grp01 <- ifelse(grp == "Pulpitis", 1, 0)

# ---------- 候选基因 ROC（最显著的 10 上 + 10 下） ----------
cand_up <- deg %>% filter(change == "Up") %>% arrange(adj.P.Val) %>% head(10)
cand_dn <- deg %>% filter(change == "Down") %>% arrange(adj.P.Val) %>% head(10)
cand <- rbind(cand_up, cand_dn)

auc_res <- data.frame()
roc_list <- list()
for (g in cand$gene) {
  if (!(g %in% rownames(expr_gene))) next
  r <- roc(grp01, as.numeric(expr_gene[g, ]), quiet = TRUE)
  auc_res <- rbind(auc_res, data.frame(gene = g,
                                       logFC = cand$logFC[cand$gene == g],
                                       AUC = as.numeric(auc(r)),
                                       stringsAsFactors = FALSE))
  roc_list[[g]] <- r
}
auc_res <- auc_res[order(-auc_res$AUC), ]
write.csv(auc_res, file.path(res_dir, "ROC_AUC_candidates.csv"), row.names = FALSE)
print(auc_res)

# ---------- ROC 曲线图（选 AUC 最高的 6 个 + ICAM1 等机制基因） ----------
show_genes <- unique(c(auc_res$gene[1:6], intersect(c("ICAM1","IL6","CXCL10","PTGS2","CD86","TNF"), rownames(expr_gene))))
roc_df <- do.call(rbind, lapply(show_genes, function(g) {
  r <- roc(grp01, as.numeric(expr_gene[g, ]), quiet = TRUE)
  data.frame(gene = g, sens = r$sensitivities, spec = r$specificities)
}))
p_roc <- ggplot(roc_df, aes(x = 1 - spec, y = sens, color = gene)) +
  geom_line(linewidth = 0.9) +
  geom_abline(linetype = "dashed", color = "grey50") +
  labs(title = "ROC curves of candidate pulpitis biomarkers (GSE77459)",
       x = "1 - Specificity", y = "Sensitivity") +
  theme_bw() + theme(legend.position = "right")
ggsave(file.path(res_dir, "ROC_curves.png"), p_roc, width = 7.5, height = 5.5, dpi = 300)

# ---------- Top 50 DEG 热图 ----------
top50 <- deg %>% filter(change != "NS") %>% arrange(adj.P.Val) %>% head(50)
mat <- expr_gene[top50$gene, ]
mat <- t(scale(t(mat)))
ann_col <- data.frame(Group = grp, row.names = colnames(mat))
p_hm <- pheatmap(mat, annotation_col = ann_col, show_rownames = TRUE,
                 fontsize_row = 7, fontsize_col = 8,
                 annotation_colors = list(Group = c(Normal = "#4575b4", Pulpitis = "#d73027")),
                 main = "Top 50 DEGs (GSE77459)", fontsize = 8, silent = TRUE)
ggsave(file.path(res_dir, "top50DEG_heatmap.png"), p_hm, width = 6.5, height = 9, dpi = 300)

cat("\nDONE\n")
