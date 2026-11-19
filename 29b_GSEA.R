# 29b_GSEA.R : genome-wide GSEA with local Hallmark GMT + GO BP gene sets
suppressMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(GO.db); library(AnnotationDbi)
  library(dplyr); library(ggplot2); library(enrichplot)
})

res_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"
deep_dir <- file.path(res_dir, "deep")

# ---- ranking ----
deg <- read.csv(file.path(res_dir, "GSE77459_DEG_full.csv"))
deg <- deg[!is.na(deg$gene) & !duplicated(deg$gene), ]
geneList <- deg$logFC
names(geneList) <- deg$gene
geneList <- sort(geneList, decreasing = TRUE)
cat("Ranked genes:", length(geneList), "\n")

# ---- Hallmark GMT ----
gmt <- readLines(file.path(deep_dir, "h.all.v2024.1.Hs.symbols.gmt"))
split_gmt <- strsplit(gmt, "\t")
t2g_h <- do.call(rbind, lapply(split_gmt, function(x) {
  data.frame(term = x[1], gene = x[-(1:2)], stringsAsFactors = FALSE)
}))
cat("Hallmark terms:", length(unique(t2g_h$term)), " pairs:", nrow(t2g_h), "\n")

gsea_h <- GSEA(geneList, TERM2GENE = t2g_h, minGSSize = 10, maxGSSize = 500,
               pvalueCutoff = 0.05, eps = 0, verbose = FALSE)
dh <- as.data.frame(gsea_h)
write.csv(dh, file.path(deep_dir, "GSEA_Hallmark.csv"), row.names = FALSE)
cat("Hallmark significant:", nrow(dh), "\n")

# ---- GO BP ----
go_terms <- AnnotationDbi::select(GO.db, keys = keys(GO.db, keytype = "GOID"),
                                  columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
go_terms <- go_terms[go_terms$ONTOLOGY == "BP", ]
t2g_go <- AnnotationDbi::select(org.Hs.eg.db,
                                keys = keys(org.Hs.eg.db, keytype = "ENTREZID"),
                                columns = c("GO", "SYMBOL"), keytype = "ENTREZID")
t2g_go <- t2g_go[!is.na(t2g_go$GO) & !is.na(t2g_go$SYMBOL), ]
t2g_go <- t2g_go[t2g_go$GO %in% go_terms$GOID, c("GO", "SYMBOL")]
colnames(t2g_go) <- c("term", "gene")
t2g_go <- unique(t2g_go)
cat("GO BP pairs:", nrow(t2g_go), "\n")

gsea_go <- GSEA(geneList, TERM2GENE = t2g_go, minGSSize = 15, maxGSSize = 500,
                pvalueCutoff = 0.05, eps = 0, verbose = FALSE)
dg <- as.data.frame(gsea_go)
write.csv(dg, file.path(deep_dir, "GSEA_GOBP.csv"), row.names = FALSE)
cat("GO BP significant:", nrow(dg), "\n")

# ---- figures ----
dh <- dh[order(-dh$NES), ]
top_h <- rbind(head(dh, 8), tail(dh, 8))
top_h$label <- gsub("^HALLMARK_", "", top_h$ID)
top_h$label <- factor(top_h$label, levels = top_h$label[order(top_h$NES)])
p1 <- ggplot(top_h, aes(label, NES, fill = NES > 0)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#2980b9"),
                    breaks = c("TRUE", "FALSE"),
                    labels = c("Activated in pulpitis", "Suppressed in pulpitis")) +
  labs(x = "", y = "Normalized enrichment score (NES)", fill = "",
       title = "Hallmark gene sets in irreversible pulpitis (GSEA)") +
  theme_bw() + theme(axis.text.y = element_text(size = 7.5), legend.position = "bottom")
ggsave(file.path(deep_dir, "GSEA_Hallmark_bar.png"), p1, width = 8, height = 6.5, dpi = 300)

dg <- dg[order(-dg$NES), ]
top_g <- rbind(head(dg, 10), tail(dg, 10))
top_g$label <- substr(top_g$Description, 1, 55)
top_g$label <- factor(top_g$label, levels = top_g$label[order(top_g$NES)])
p2 <- ggplot(top_g, aes(label, NES, fill = NES > 0)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#2980b9"),
                    breaks = c("TRUE", "FALSE"),
                    labels = c("Activated in pulpitis", "Suppressed in pulpitis")) +
  labs(x = "", y = "Normalized enrichment score (NES)", fill = "",
       title = "GO biological processes in irreversible pulpitis (GSEA)") +
  theme_bw() + theme(axis.text.y = element_text(size = 6.5), legend.position = "bottom")
ggsave(file.path(deep_dir, "GSEA_GOBP_bar.png"), p2, width = 8.5, height = 7, dpi = 300)

if (nrow(dh) >= 3) {
  ids <- dh$ID[order(-abs(dh$NES))][1:3]
  p3 <- gseaplot2(gsea_h, geneSetID = ids, title = "Top Hallmark gene sets", pvalue_table = TRUE)
  ggsave(file.path(deep_dir, "GSEA_top3_curves.png"), p3, width = 10, height = 8, dpi = 300)
}
cat("DONE\n")
