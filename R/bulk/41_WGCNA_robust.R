suppressPackageStartupMessages(library(WGCNA))
options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 2)
outdir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/deep"

expr <- read.csv("C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/GSE77459_expr_gene.csv",
                 row.names = 1, check.names = FALSE)
datExpr <- t(expr)                      # 12 samples x genes
group <- factor(ifelse(grepl("pulpitis|inflamed|P", rownames(datExpr), ignore.case = TRUE), 1, 0))
grp <- read.csv("C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/GSE77459_groups.csv", stringsAsFactors = FALSE)
trait <- as.numeric(grp$group[match(rownames(datExpr), grp$sample)] == "Pulpitis")
names(trait) <- rownames(datExpr)

ass <- read.csv("C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/deep/WGCNA_module_assignment.csv", stringsAsFactors = FALSE)
m1_genes <- ass$gene[ass$module == "M1"]
m2_genes <- ass$gene[ass$module == "M2"]
cat("reference M1:", length(m1_genes), "M2:", length(m2_genes), "\n")

run_wgcna <- function(datExpr, trait, n_top, power, minSize, deepSplit, mergeH) {
  mads <- apply(datExpr, 2, mad)
  keep <- order(-mads)[1:n_top]
  net <- blockwiseModules(datExpr[, keep, drop = FALSE], power = power,
                          TOMType = "signed", minModuleSize = minSize,
                          mergeCutHeight = mergeH, deepSplit = deepSplit,
                          numericLabels = TRUE, pamRespectsDendro = FALSE,
                          saveTOMs = FALSE, saveTOMFileBase = tempfile(),
                          verbose = 0)
  mods <- net$colors
  res <- list()
  for (m in setdiff(sort(unique(mods)), 0)) {
    g <- colnames(datExpr)[keep][mods == m]
    j1 <- length(intersect(g, m1_genes)) / length(union(g, m1_genes))
    j2 <- length(intersect(g, m2_genes)) / length(union(g, m2_genes))
    me <- net$MEs[, paste0("ME", m)]
    r <- cor(me, trait); p <- corPvalueStudent(r, nrow(datExpr))
    res[[as.character(m)]] <- data.frame(module = m, ngenes = length(g), jac_M1 = j1, jac_M2 = j2, r = r, p = p)
  }
  do.call(rbind, res)
}

variants <- list(
  original = list(6000, 5, 30, 3, 0.30),
  power7   = list(6000, 7, 30, 3, 0.30),
  min50    = list(6000, 5, 50, 3, 0.30),
  deep2    = list(6000, 5, 30, 2, 0.30),
  top5500  = list(5500, 5, 30, 3, 0.30)
)

out <- list()
for (v in names(variants)) {
  cat("running variant:", v, "\n")
  prm <- variants[[v]]
  r <- tryCatch(run_wgcna(datExpr, trait, prm[[1]], prm[[2]], prm[[3]], prm[[4]], prm[[5]]),
                error = function(e) { cat("error:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) next
  m1_like <- r[which.max(r$jac_M1), ]
  r_neg <- r[r$r < 0, , drop = FALSE]
  if (nrow(r_neg) > 0) {
    m2_like <- r_neg[which.max(r_neg$jac_M2), ]
  } else {
    m2_like <- data.frame(module = NA, ngenes = NA, jac_M2 = NA, r = NA, p = NA)
  }
  out[[v]] <- data.frame(variant = v,
                         n_modules = nrow(r),
                         M1like_module = m1_like$module, M1like_ngenes = m1_like$ngenes,
                         M1like_jac = round(m1_like$jac_M1, 3), M1like_r = round(m1_like$r, 3), M1like_p = signif(m1_like$p, 3),
                         M2like_module = m2_like$module, M2like_ngenes = m2_like$ngenes,
                         M2like_jac = round(m2_like$jac_M2, 3), M2like_r = round(m2_like$r, 3), M2like_p = signif(m2_like$p, 3))
}
res <- do.call(rbind, out)
write.csv(res, file.path(outdir, "WGCNA_robustness_variants.csv"), row.names = FALSE)
print(res)
cat("\nDONE\n")
