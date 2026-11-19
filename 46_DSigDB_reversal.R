suppressPackageStartupMessages({library(fgsea); library(ggplot2)})
outdir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/deep"
res_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"

gmt_lines <- readLines(file.path(outdir, "DSigDB.gmt"))
dsig <- setNames(lapply(strsplit(gmt_lines, "\t"), function(v) v[-1]),
                 sapply(strsplit(gmt_lines, "\t"), `[`, 1))
cat("DSigDB sets:", length(dsig), "\n")

deg <- read.csv(file.path(res_dir, "GSE77459_DEG_full.csv"), stringsAsFactors = FALSE)
deg <- deg[!is.na(deg$logFC) & !is.na(deg$t), ]
deg <- deg[order(-abs(deg$t)), ]
deg <- deg[!duplicated(deg$gene), ]
rk <- setNames(deg$logFC, deg$gene)

fg <- fgseaMultilevel(dsig, rk, minSize = 15, maxSize = 500)
res <- as.data.frame(fg)[, c("pathway", "pval", "padj", "NES", "size")]
res$pval[is.na(res$pval)] <- 1; res$padj[is.na(res$padj)] <- 1
write.csv(res, file.path(outdir, "DSigDB_drug_fgsea.csv"), row.names = FALSE)

sig <- res[res$padj < 0.05, ]
cat("significant drugs:", nrow(sig), "\n")
top_neg <- head(sig[order(sig$NES), ], 20)
top_pos <- head(sig[order(-sig$NES), ], 20)
cat("== top negative NES (drug signature opposite to pulpitis program) ==\n")
print(top_neg[, c("pathway", "NES", "padj", "size")], row.names = FALSE)
cat("== top positive NES (drug signature aligned with pulpitis program) ==\n")
print(top_pos[, c("pathway", "NES", "padj", "size")], row.names = FALSE)

sel <- rbind(data.frame(name = top_neg$pathway, NES = top_neg$NES, dir = "Opposite to pulpitis"),
             data.frame(name = top_pos$pathway, NES = top_pos$NES, dir = "Aligned with pulpitis"))
sel <- sel[order(sel$NES), ]
sel$name <- factor(sel$name, levels = sel$name)
pp <- ggplot(sel, aes(x = NES, y = name, fill = dir)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Opposite to pulpitis" = "#4575B4", "Aligned with pulpitis" = "#D73027")) +
  labs(title = "DSigDB drug perturbation signatures versus pulpitis transcriptome (GSE77459)",
       x = "NES", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(outdir, "DSigDB_drug_top40.png"), pp, width = 10, height = 10, dpi = 160)
cat("\nDONE\n")
