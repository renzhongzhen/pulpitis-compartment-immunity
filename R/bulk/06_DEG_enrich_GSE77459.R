suppressMessages({
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
})

data_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/data"
res_dir  <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"
load(file.path(data_dir, "GSE77459_genelevel.RData"))  # expr_gene, grp

grp <- factor(grp, levels = c("Normal", "Pulpitis"))

# ---------- limma DEG ----------
design <- model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
fit <- lmFit(expr_gene, design)
contr <- makeContrasts(Pulpitis - Normal, levels = design)
fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)
deg <- topTable(fit2, number = Inf, sort.by = "none")
deg$gene <- rownames(deg)

fc_cut <- 1; padj_cut <- 0.05
deg$change <- "NS"
deg$change[deg$logFC >  fc_cut & deg$adj.P.Val < padj_cut] <- "Up"
deg$change[deg$logFC < -fc_cut & deg$adj.P.Val < padj_cut] <- "Down"
deg$change <- factor(deg$change, levels = c("Up", "Down", "NS"))
cat("DEG counts:\n"); print(table(deg$change))
write.csv(deg, file.path(res_dir, "GSE77459_DEG_full.csv"), row.names = FALSE)

# ---------- Volcano ----------
top_label <- deg %>% filter(change != "NS") %>% arrange(adj.P.Val) %>% head(10)
p_vol <- ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
  geom_point(size = 0.6, alpha = 0.7) +
  scale_color_manual(values = c("Up" = "#d73027", "Down" = "#4575b4", "NS" = "grey80")) +
  geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "grey40") +
  geom_text(data = top_label, aes(label = gene), size = 2.5, color = "black",
            nudge_y = 0.5, check_overlap = TRUE) +
  labs(title = "GSE77459  Pulpitis vs Normal",
       x = "log2 fold change", y = "-log10(adjusted P)") +
  theme_bw() + theme(legend.position = "top")
ggsave(file.path(res_dir, "volcano_GSE77459.png"), p_vol, width = 7, height = 6, dpi = 300)

up_sym <- deg$gene[deg$change == "Up"]; dn_sym <- deg$gene[deg$change == "Down"]

# ---------- GO BP ----------
go_up <- enrichGO(up_sym, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
go_dn <- enrichGO(dn_sym, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)

if (!is.null(go_up) && nrow(as.data.frame(go_up)) > 0) {
  write.csv(as.data.frame(go_up), file.path(res_dir, "GO_BP_up.csv"), row.names = FALSE)
  ggsave(file.path(res_dir, "GO_BP_up.png"),
         dotplot(go_up, showCategory = 15) + ggtitle("GO BP - Upregulated in Pulpitis"),
         width = 8, height = 7, dpi = 300)
}
if (!is.null(go_dn) && nrow(as.data.frame(go_dn)) > 0) {
  write.csv(as.data.frame(go_dn), file.path(res_dir, "GO_BP_down.csv"), row.names = FALSE)
  ggsave(file.path(res_dir, "GO_BP_down.png"),
         dotplot(go_dn, showCategory = 15) + ggtitle("GO BP - Downregulated in Pulpitis"),
         width = 8, height = 7, dpi = 300)
}

# ---------- KEGG：在线优先，失败则离线 org.Hs.egPATH2EG ----------
eg_up <- bitr(up_sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
eg_dn <- bitr(dn_sym, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

do_kegg <- function(egdf, label) {
  res <- NULL
  online <- tryCatch(
    enrichKEGG(egdf$ENTREZID, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) { cat("  KEGG online failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(online)) { res <- online; cat("  KEGG online OK\n") }
  # 离线兜底
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    p2e <- as.list(org.Hs.egPATH2EG)
    t2g <- data.frame(term = rep(names(p2e), lengths(p2e)),
                      gene = unlist(p2e), stringsAsFactors = FALSE)
    off <- enricher(egdf$ENTREZID, TERM2GENE = t2g, pvalueCutoff = 0.05,
                    pAdjustMethod = "BH", qvalueCutoff = 0.2)
    res <- off; cat("  KEGG offline fallback used\n")
  }
  res
}

kegg_up <- do_kegg(eg_up, "up")
kegg_dn <- do_kegg(eg_dn, "down")

if (!is.null(kegg_up) && nrow(as.data.frame(kegg_up)) > 0) {
  write.csv(as.data.frame(kegg_up), file.path(res_dir, "KEGG_up.csv"), row.names = FALSE)
  ggsave(file.path(res_dir, "KEGG_up.png"),
         dotplot(kegg_up, showCategory = 15) + ggtitle("KEGG - Upregulated in Pulpitis"),
         width = 8, height = 7, dpi = 300)
}
if (!is.null(kegg_dn) && nrow(as.data.frame(kegg_dn)) > 0) {
  write.csv(as.data.frame(kegg_dn), file.path(res_dir, "KEGG_down.csv"), row.names = FALSE)
  ggsave(file.path(res_dir, "KEGG_down.png"),
         dotplot(kegg_dn, showCategory = 15) + ggtitle("KEGG - Downregulated in Pulpitis"),
         width = 8, height = 7, dpi = 300)
}

saveRDS(list(go_up = go_up, go_dn = go_dn, kegg_up = kegg_up, kegg_dn = kegg_dn),
        file.path(res_dir, "enrichment_GSE77459.rds"))

cat("Top 20 up genes:\n"); print(head(deg$gene[deg$change == "Up"], 20))
cat("Top 20 down genes:\n"); print(head(deg$gene[deg$change == "Down"], 20))
cat("\nDONE\n")
