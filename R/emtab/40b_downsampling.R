suppressPackageStartupMessages({library(Seurat); library(Matrix)})
p <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/ea960c355a014573b0e904afc349b633processedmatrix.rds"
outdir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"

x <- readRDS(p)
md <- x@meta.data
lv <- levels(factor(as.character(md$samples)))
stage_of <- function(s) ifelse(grepl("^Healthy", s), "Healthy",
                        ifelse(grepl("^Early", s), "Early", "Advanced"))

set.seed(20261119)
nsim <- 20; n_sub <- 1000
cts_all <- levels(Idents(x))
sub_res <- list()
for (i in seq_len(nsim)) {
  keep_cells <- unlist(lapply(lv, function(s) {
    w <- which(md$samples == s)
    if (length(w) <= n_sub) w else sample(w, n_sub)
  }))
  tab <- table(factor(Idents(x)[keep_cells], levels = cts_all),
               factor(as.character(md$samples[keep_cells]), levels = lv))
  pct <- 100 * prop.table(tab, 2)
  grp <- stage_of(colnames(pct))
  sub_res[[i]] <- rowMeans(pct[, grp == "Advanced"]) - rowMeans(pct[, grp == "Healthy"])
}
submat <- do.call(rbind, sub_res)
colnames(submat) <- cts_all
sm <- colMeans(submat)
subdf <- data.frame(celltype = cts_all,
                    mean_diff = round(sm, 3),
                    lo = round(apply(submat, 2, quantile, 0.025), 3),
                    hi = round(apply(submat, 2, quantile, 0.975), 3),
                    sign_consistent = 100 * colMeans(sign(submat) == matrix(rep(sign(sm), each = nsim), nrow = nsim)))
subdf <- subdf[order(-subdf$mean_diff), ]
write.csv(subdf, file.path(outdir, "EMTAB_downsampling_robustness.csv"), row.names = FALSE)
cat("== downsampling robustness (mean Advanced-Healthy pct diff, 95% CI, 20 x 1000 cells/patient) ==\n")
print(subdf, row.names = FALSE)
cat("\nDONE\n")
