suppressPackageStartupMessages(library(dplyr))
sc_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc"
outdir <- file.path(sc_dir, "results")
d <- read.csv(file.path(outdir, "EMTAB_composition_by_sample.csv"), stringsAsFactors = FALSE,
              check.names = FALSE)
cat("samples:", length(unique(d$sample)), "celltypes:", length(unique(d$celltype)), "\n")
endo_ct <- grep("^Endothelial cells", unique(d$celltype), value = TRUE)
strom_ct <- grep("^Fibroblast & Stromal", unique(d$celltype), value = TRUE)
cat("endothelial subclusters:", paste(endo_ct, collapse = " | "), "\n")
cat("stromal subclusters:", paste(strom_ct, collapse = " | "), "\n")

pp <- d %>% group_by(sample, stage) %>% summarise(
  endo = sum(pct[celltype %in% endo_ct]),
  stroma = sum(pct[celltype %in% strom_ct]),
  endo3 = pct[celltype == "Endothelial cells 3"],
  fib1 = pct[celltype == "Fibroblast & Stromal cells 1"], .groups = "drop") %>%
  mutate(ratio = endo / stroma,
         stage = factor(stage, levels = c("Healthy", "Early", "Advanced")))
print(as.data.frame(pp), digits = 4)

kw <- function(v, g) { g <- droplevels(g); if (length(unique(g)) < 2) return(NA_real_)
  unname(kruskal.test(v ~ g)$p.value) }
stages <- levels(pp$stage)
means <- function(df, v) sapply(stages, function(s) mean(df[[v]][df$stage == s]))

res <- data.frame(metric = c("Endothelial total (%)", "Endothelial/Stromal ratio",
                             "Endothelial 3 (%)", "Fibroblast/Stromal 1 (%)"),
                  p_kruskal = c(kw(pp$endo, pp$stage), kw(pp$ratio, pp$stage),
                                kw(pp$endo3, pp$stage), kw(pp$fib1, pp$stage)))
res$mean_healthy <- c(means(pp, "endo")[1], means(pp, "ratio")[1], means(pp, "endo3")[1], means(pp, "fib1")[1])
res$mean_early <- c(means(pp, "endo")[2], means(pp, "ratio")[2], means(pp, "endo3")[2], means(pp, "fib1")[2])
res$mean_advanced <- c(means(pp, "endo")[3], means(pp, "ratio")[3], means(pp, "endo3")[3], means(pp, "fib1")[3])
# Wilcoxon advanced vs healthy (two-group complement to KW)
ah <- pp[pp$stage %in% c("Healthy", "Advanced"), ]
res$p_wilcox_adv_vs_healthy <- c(wilcox.test(endo ~ stage, data = ah)$p.value,
                                 wilcox.test(ratio ~ stage, data = ah)$p.value,
                                 wilcox.test(endo3 ~ stage, data = ah)$p.value,
                                 wilcox.test(fib1 ~ stage, data = ah)$p.value)
print(res, digits = 4)
write.csv(res, file.path(outdir, "EMTAB_endo_ratio_stats.csv"), row.names = FALSE)

# ---- leave-one-patient-out ----
vars <- c("endo", "ratio", "endo3")
lopo <- do.call(rbind, lapply(pp$sample, function(p) {
  s <- pp[pp$sample != p, ]
  s$stage <- droplevels(s$stage)
  data.frame(removed = p,
             p_endo = kw(s$endo, s$stage), p_ratio = kw(s$ratio, s$stage), p_endo3 = kw(s$endo3, s$stage),
             endo_H = mean(s$endo[s$stage == "Healthy"]), endo_A = mean(s$endo[s$stage == "Advanced"]),
             endo3_H = mean(s$endo3[s$stage == "Healthy"]), endo3_A = mean(s$endo3[s$stage == "Advanced"]),
             ratio_H = mean(s$ratio[s$stage == "Healthy"]), ratio_A = mean(s$ratio[s$stage == "Advanced"]),
             stringsAsFactors = FALSE)
}))
print(lopo, digits = 4)
write.csv(lopo, file.path(outdir, "EMTAB_LOPO_endo.csv"), row.names = FALSE)
cat("\nLOPO summary:\n")
cat("endo P range:", range(lopo$p_endo), " | ratio P range:", range(lopo$p_ratio),
    " | endo3 P range:", range(lopo$p_endo3), "\n")
cat("endo advanced-healthy difference range:", range(lopo$endo_A - lopo$endo_H), "\n")
cat("endo3 advanced-healthy difference range:", range(lopo$endo3_A - lopo$endo3_H), "\n")
cat("ratio advanced/healthy fold range:", range(lopo$ratio_A / lopo$ratio_H), "\n")
cat("\nDONE\n")
