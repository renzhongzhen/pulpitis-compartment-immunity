suppressPackageStartupMessages({library(ggplot2)})
outdir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"
delta <- read.csv(file.path(outdir, "EMTAB_LR_pair_delta.csv"), stringsAsFactors = FALSE)

# focus axes: immune <-> Endothelial 3 (ACKR1+ venule), Endo3 -> matrix, immune -> matrix
axes <- list(
  "Immune -> Endo3"  = delta[grepl("Blood & Immune", delta$from) & delta$to == "Endothelial cells 3", ],
  "Endo3 -> Immune"  = delta[delta$from == "Endothelial cells 3" & grepl("Blood & Immune", delta$to), ],
  "Immune -> Endo1"  = delta[grepl("Blood & Immune", delta$from) & delta$to == "Endothelial cells 1", ],
  "Endo1 -> Immune"  = delta[delta$from == "Endothelial cells 1" & grepl("Blood & Immune", delta$to), ],
  "Immune -> Fib1"   = delta[grepl("Blood & Immune", delta$from) & delta$to == "Fibroblast & Stromal cells 1", ],
  "Fib1 -> Immune"   = delta[delta$from == "Fibroblast & Stromal cells 1" & grepl("Blood & Immune", delta$to), ],
  "Fib1 -> Endo3"    = delta[delta$from == "Fibroblast & Stromal cells 1" & delta$to == "Endothelial cells 3", ],
  "Immune -> Glial2" = delta[grepl("Blood & Immune", delta$from) & delta$to == "Glial cells 2", ]
)
ax <- do.call(rbind, lapply(names(axes), function(nm) {
  d <- axes[[nm]]
  if (nrow(d) == 0) return(NULL)
  data.frame(axis = nm, healthy = sum(d$healthy), early = sum(d$early), advanced = sum(d$advanced),
             delta = sum(d$delta_adv))
}))
ax <- ax[order(-ax$delta), ]
write.csv(ax, file.path(outdir, "EMTAB_LR_focus_axes.csv"), row.names = FALSE)
cat("== communication share by axis (normalized), Advanced minus Healthy ==\n")
print(ax, row.names = FALSE)

ld <- reshape(ax, direction = "long", varying = list(c("healthy", "early", "advanced")),
              v.names = "share", times = c("Healthy", "Early", "Advanced"), timevar = "stage")
ld$stage <- factor(ld$stage, levels = c("Healthy", "Early", "Advanced"))
ld$axis <- factor(ld$axis, levels = ax$axis)
pp <- ggplot(ld, aes(x = stage, y = share, group = axis, color = axis)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  labs(title = "Ligand-receptor communication share on key axes (E-MTAB-14896, normalized)",
       x = NULL, y = "communication share") +
  theme_bw(base_size = 10) + theme(legend.position = "right", legend.title = element_blank())
ggsave(file.path(outdir, "EMTAB_LR_focus_axes.png"), pp, width = 8.5, height = 5, dpi = 160)
cat("\nDONE\n")
