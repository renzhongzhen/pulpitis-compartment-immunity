suppressPackageStartupMessages({library(ggplot2)})
outdir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/deep"
res <- read.csv(file.path(outdir, "DSigDB_drug_fgsea.csv"), stringsAsFactors = FALSE)

res$is_up <- grepl(" UP$", res$pathway)
res$is_down <- grepl(" DOWN$", res$pathway)
# opposing = drug effect opposite to the pulpitis program:
#   drug-DOWN sets enriched at the disease UP end (NES>0), drug-UP sets enriched at the disease DOWN end (NES<0)
# aligned  = drug effect parallel to the pulpitis program
opposing <- res[(res$is_down & res$NES > 1.5) | (res$is_up & res$NES < -1.5) & res$padj < 0.05, ]
opposing <- res[((res$is_down & res$NES > 1.5) | (res$is_up & res$NES < -1.5)) & res$padj < 0.05, ]
aligned  <- res[((res$is_up & res$NES > 1.5) | (res$is_down & res$NES < -1.5)) & res$padj < 0.05, ]
opposing$drug <- sub(" (MCF7|PC3|HL60|A549|HEPG2|HCC515|HA1E|SKB).*$", "", opposing$pathway)
opposing$drug <- sub(" (CTD|TTD|LINCS|Kinome Scan).*$", "", opposing$drug)
aligned$drug <- sub(" (MCF7|PC3|HL60|A549|HEPG2|HCC515|HA1E|SKB).*$", "", aligned$pathway)
aligned$drug <- sub(" (CTD|TTD|LINCS|Kinome Scan).*$", "", aligned$drug)

opp_agg <- aggregate(cbind(NES, padj) ~ drug, opposing, function(v) v[which.max(abs(v))])
opp_agg$padj <- aggregate(padj ~ drug, opposing, min)$padj
opp_top <- head(opp_agg[order(-opp_agg$NES), ], 20)
ali_top <- head(aligned[order(-aligned$NES), ], 20)
cat("== drugs OPPOSING the pulpitis program (candidate therapeutic direction) ==\n")
print(opp_top, row.names = FALSE)
cat("== drugs ALIGNED with the pulpitis program (pro-inflammatory / caution) ==\n")
print(ali_top[, c("drug", "NES", "padj", "size", "pathway")], row.names = FALSE)

write.csv(opp_agg, file.path(outdir, "DSigDB_opposing_drugs.csv"), row.names = FALSE)
ali_out <- do.call(rbind, lapply(split(aligned, aligned$drug), function(d) d[which.max(abs(d$NES)), c("drug", "NES", "padj", "size", "pathway")]))
write.csv(ali_out, file.path(outdir, "DSigDB_aligned_drugs.csv"), row.names = FALSE)

sel <- rbind(data.frame(name = opp_top$drug, NES = opp_top$NES, dir = "Opposing pulpitis"),
             data.frame(name = ali_out$drug[1:20], NES = ali_out$NES[1:20], dir = "Aligned with pulpitis"))
sel <- sel[!duplicated(sel$name), ]
sel <- sel[order(sel$NES), ]
sel$name <- factor(sel$name, levels = sel$name)
pp <- ggplot(sel, aes(x = NES, y = name, fill = dir)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Opposing pulpitis" = "#4575B4", "Aligned with pulpitis" = "#D73027")) +
  labs(title = "DSigDB drug signatures versus the pulpitis transcriptome (direction-aware)",
       x = "NES", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(outdir, "DSigDB_drug_direction_top40.png"), pp, width = 10, height = 10, dpi = 160)
cat("\nDONE\n")
