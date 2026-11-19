suppressMessages({ library(ggplot2); library(dplyr) })

data_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/data"
res_dir  <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"

p <- read.csv(file.path(res_dir, "PBM_panel_dCt.csv"))
d <- p[!is.na(p$LPS_dCt) & !is.na(p$PBM_dCt), ]
cat("Genes with both effects:", nrow(d), "\n")

# 逆转判定：PBM 与 LPS 效应方向相反
d$rev <- sign(d$PBM_dCt) * sign(d$LPS_dCt) < 0
n_rev <- sum(d$rev)
cat("Reversed:", n_rev, "/", nrow(d), "\n")

# 二项检验（H0: 逆转概率 0.5）
bt <- binom.test(n_rev, nrow(d), p = 0.5)
cat("Binomial test p =", bt$p.value, "\n")

# 象限分类
d$class <- ifelse(d$LPS_dCt < 0 & d$PBM_dCt > 0, "PBM suppresses LPS-upregulated genes",
           ifelse(d$LPS_dCt > 0 & d$PBM_dCt < 0, "PBM restores LPS-downregulated genes",
           ifelse(d$LPS_dCt < 0 & d$PBM_dCt < 0, "PBM further increases", "PBM further decreases")))
write.csv(d[order(-d$rev), c("gene","LPS_dCt","PBM_dCt","class","rev")],
          file.path(res_dir, "PBM_LPS_reversal.csv"), row.names = FALSE)
print(table(d$class))

# ---------- 图1：逆转散点 ----------
d$cls2 <- ifelse(d$rev, "Reversed", "Not reversed")
p_scat <- ggplot(d, aes(x = LPS_dCt, y = PBM_dCt, color = cls2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("Reversed" = "#2ca02c", "Not reversed" = "grey70")) +
  annotate("text", x = -3.5, y = 3.2, size = 3.6, hjust = 0,
           label = "PBM suppresses LPS-upregulated genes") +
  annotate("text", x = -3.5, y = -3.4, size = 3.6, hjust = 0,
           label = "PBM further increases") +
  annotate("text", x = 0.8, y = 3.4, size = 3.6, hjust = 0,
           label = "PBM further decreases") +
  annotate("text", x = 0.8, y = -3.2, size = 3.6, hjust = 0,
           label = "PBM restores LPS-downregulated genes") +
  labs(title = paste0("PBM reverses LPS-induced changes in ", n_rev, "/", nrow(d),
                      " genes (", round(100 * n_rev / nrow(d)), "%)\nbinomial test p = ", signif(bt$p.value, 3)),
       x = expression(Delta*Ct[LPS]~"(LPS vs Control)"),
       y = expression(Delta*Ct[PBM]~"(LPS+Laser vs LPS)")) +
  theme_bw() + theme(legend.position = "top")
ggsave(file.path(res_dir, "PBM_LPS_reversal_scatter.png"), p_scat, width = 7.5, height = 6.5, dpi = 300)

# ---------- 图2：关键基因柱状（LPS 效应 vs PBM 效应） ----------
key <- c("VDR","MMP8","MMP10","BMP4","VEGFB","ALPL","TWIST1","RUNX2","CTSK","TNFSF11","FGFR2","NFKB1")
kd <- d[d$gene %in% key, ]
kd$gene <- factor(kd$gene, levels = rev(key))
kl <- tidyr::pivot_longer(kd[, c("gene","LPS_dCt","PBM_dCt")], cols = c(LPS_dCt, PBM_dCt),
                          names_to = "effect", values_to = "dCt")
kl$effect <- factor(kl$effect, levels = c("LPS_dCt","PBM_dCt"),
                    labels = c("LPS effect","PBM effect"))
p_bar <- ggplot(kl, aes(x = gene, y = dCt, fill = effect)) +
  geom_col(position = "dodge", width = 0.65) +
  coord_flip() +
  geom_hline(yintercept = 0, color = "grey40") +
  scale_fill_manual(values = c("LPS effect" = "#bdbdbd", "PBM effect" = "#e377c2")) +
  labs(title = "Effects of LPS and PBM on key pulpitis-related genes (DPSCs, 24h)",
       y = expression(Delta*Ct), x = "") +
  theme_bw()
ggsave(file.path(res_dir, "PBM_key_genes_bar.png"), p_bar, width = 7, height = 5.5, dpi = 300)

cat("\nKey gene table:\n")
print(kd[, c("gene","LPS_dCt","PBM_dCt","class")])
cat("\nDONE\n")
