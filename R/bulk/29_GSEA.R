# 29_GSEA.R : genome-wide GSEA of GSE77459 (Hallmark + KEGG), ranked by limma logFC
suppressMessages({
  library(clusterProfiler); library(msigdbr); library(dplyr)
  library(ggplot2); library(enrichplot)
})

res_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"
deep_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/deep"
dir.create(deep_dir, showWarnings = FALSE, recursive = TRUE)

deg <- read.csv(file.path(res_dir, "GSE77459_DEG_full.csv"))
deg <- deg[!is.na(deg$gene) & !duplicated(deg$gene), ]
geneList <- deg$logFC
names(geneList) <- deg$gene
geneList <- sort(geneList, decreasing = TRUE)
cat("Ranked genes:", length(geneList), "\n")

hs_h <- msigdbr(species = "Homo sapiens", collection = "H")
cat("Hallmark terms:", length(unique(hs_h$gs_name)), "\n")

hs_k <- tryCatch(
  msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY"),
  error = function(e) msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS")
)
cat("KEGG terms:", length(unique(hs_k$gs_name)), "\n")

gsea_h <- GSEA(geneList, TERM2GENE = hs_h[, c("gs_name", "gene_symbol")],
               minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0, verbose = FALSE)
gsea_k <- GSEA(geneList, TERM2GENE = hs_k[, c("gs_name", "gene_symbol")],
               minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, eps = 0, verbose = FALSE)

write.csv(as.data.frame(gsea_h), file.path(deep_dir, "GSEA_Hallmark.csv"), row.names = FALSE)
write.csv(as.data.frame(gsea_k), file.path(deep_dir, "GSEA_KEGG.csv"), row.names = FALSE)
cat("Hallmark significant:", nrow(as.data.frame(gsea_h)), " KEGG significant:", nrow(as.data.frame(gsea_k)), "\n")

dh <- as.data.frame(gsea_h)
dh <- dh[order(-dh$NES), ]
top_h <- rbind(head(dh, 8), tail(dh, 8))
top_h$gs_name <- factor(top_h$gs_name, levels = top_h$gs_name[order(top_h$NES)])
p1 <- ggplot(top_h, aes(x = gs_name, y = NES, fill = NES > 0)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#2980b9"),
                    breaks = c("TRUE", "FALSE"),
                    labels = c("Activated in pulpitis", "Suppressed in pulpitis")) +
  labs(x = "", y = "Normalized enrichment score (NES)", fill = "",
       title = "Hallmark gene sets enriched in irreversible pulpitis (GSEA)") +
  theme_bw() + theme(axis.text.y = element_text(size = 7), legend.position = "bottom")
ggsave(file.path(deep_dir, "GSEA_Hallmark_bar.png"), p1, width = 8, height = 6, dpi = 300)

dk <- as.data.frame(gsea_k)
dk <- dk[order(-dk$NES), ]
top_k <- rbind(head(dk, 8), tail(dk, 8))
top_k$gs_name <- gsub("^KEGG_LEGACY_|^KEGG_MEDICUS_", "", top_k$gs_name)
top_k$gs_name <- factor(top_k$gs_name, levels = top_k$gs_name[order(top_k$NES)])
p2 <- ggplot(top_k, aes(x = gs_name, y = NES, fill = NES > 0)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#2980b9"),
                    breaks = c("TRUE", "FALSE"),
                    labels = c("Activated in pulpitis", "Suppressed in pulpitis")) +
  labs(x = "", y = "Normalized enrichment score (NES)", fill = "",
       title = "KEGG pathways enriched in irreversible pulpitis (GSEA)") +
  theme_bw() + theme(axis.text.y = element_text(size = 7), legend.position = "bottom")
ggsave(file.path(deep_dir, "GSEA_KEGG_bar.png"), p2, width = 8, height = 6, dpi = 300)

if (nrow(dh) >= 3) {
  ids <- dh$ID[order(-abs(dh$NES))][1:3]
  p3 <- gseaplot2(gsea_h, geneSetID = ids, title = "Top Hallmark gene sets", pvalue_table = TRUE)
  ggsave(file.path(deep_dir, "GSEA_top3_curves.png"), p3, width = 10, height = 8, dpi = 300)
}
cat("DONE\n")
