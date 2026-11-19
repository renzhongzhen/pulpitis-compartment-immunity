suppressPackageStartupMessages({
  library(dorothea); library(decoupleR); library(dplyr)
  library(Seurat); library(Matrix); library(org.Hs.eg.db); library(AnnotationDbi)
  library(ggplot2)
})
res_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"
sc_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc"

net <- as.data.frame(dorothea_hs)
net <- net[net$confidence %in% c("A", "B"), c("tf", "target", "mor")]
colnames(net) <- c("source", "target", "mor")
cat("regulon: TFs", length(unique(net$source)), " edges", nrow(net), "\n")

# ---- 1. GSE77459 bulk ----
expr <- read.csv(file.path(res_dir, "GSE77459_expr_gene.csv"), row.names = 1, check.names = FALSE)
grp <- read.csv(file.path(res_dir, "GSE77459_groups.csv"))
g <- grp$group[match(colnames(expr), grp$sample)]
cat("bulk samples:", ncol(expr), "\n")

tf_act <- run_ulm(mat = as.matrix(expr), network = net,
                  .source = "source", .target = "target", .mor = "mor", minsize = 5)
tf_bulk <- as.data.frame(tf_act)
write.csv(tf_bulk, file.path(res_dir, "deep", "TF_activity_GSE77459.csv"), row.names = FALSE)

comp <- do.call(rbind, lapply(unique(tf_bulk$source), function(tf) {
  v <- tf_bulk$score[tf_bulk$source == tf]
  if (length(v) != 12) return(NULL)
  x <- v[g == "Pulpitis"]; y <- v[g == "Normal"]
  data.frame(TF = tf, mean_normal = mean(y), mean_pulpitis = mean(x),
             diff = mean(x) - mean(y),
             p = tryCatch(wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA))
}))
comp$padj <- p.adjust(comp$p, "BH")
comp <- comp[order(-abs(comp$diff)), ]
write.csv(comp, file.path(res_dir, "deep", "TF_activity_compare.csv"), row.names = FALSE)
cat("== bulk TF activity (top 15 by |diff|) ==\n")
print(head(comp[, c("TF", "mean_normal", "mean_pulpitis", "diff", "padj")], 15), row.names = FALSE)

top_tf <- head(comp[order(-comp$diff), ], 10)$TF
bot_tf <- head(comp[order(comp$diff), ], 10)$TF
sel <- comp[comp$TF %in% c(top_tf, bot_tf), ]
sel$TF <- factor(sel$TF, levels = rev(sel$TF[order(sel$diff)]))
sel$dir <- ifelse(sel$diff > 0, "Up in pulpitis", "Down in pulpitis")
pp <- ggplot(sel, aes(x = diff, y = TF, fill = dir)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Up in pulpitis" = "#D73027", "Down in pulpitis" = "#4575B4")) +
  labs(title = "TF activity difference (pulpitis minus normal, GSE77459; DoRothEA A+B, ULM)",
       x = "activity difference", y = NULL) +
  theme_bw(base_size = 10) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(res_dir, "deep", "TF_activity_top20.png"), pp, width = 6.5, height = 5, dpi = 160)

# ---- 2. E-MTAB-14896 three stages (patient pseudobulk) ----
p <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/ea960c355a014573b0e904afc349b633processedmatrix.rds"
x <- readRDS(p)
md <- x@meta.data
cnt <- LayerData(x, assay = "RNA", layer = "counts")
ens <- sub("\\..*$", "", rownames(cnt))
sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
keep <- !is.na(sym) & !duplicated(sym)
cnt2 <- cnt[keep, , drop = FALSE]; rownames(cnt2) <- sym[keep]
lv <- levels(factor(as.character(md$samples)))
stage_of <- function(s) ifelse(grepl("^Healthy", s), "Healthy",
                        ifelse(grepl("^Early", s), "Early", "Advanced"))
stage <- factor(stage_of(lv), levels = c("Healthy", "Early", "Advanced"))
pat <- factor(as.character(md$samples), levels = lv)
ind <- sparseMatrix(i = seq_along(pat), j = as.integer(pat), x = 1, dims = c(length(pat), 12))
pb <- as.matrix(cnt2 %*% ind)
colnames(pb) <- lv
lcpm <- log2(sweep(pb, 2, colSums(pb), "/") * 1e6 + 1)

tf2 <- run_ulm(mat = lcpm, network = net, .source = "source", .target = "target", .mor = "mor", minsize = 5)
tf_emt <- as.data.frame(tf2)
tf_emt$stage <- stage[match(tf_emt$condition, lv)]
write.csv(tf_emt, file.path(sc_dir, "results", "EMTAB_TF_activity.csv"), row.names = FALSE)

focus <- c("NFKB1", "RELA", "STAT1", "STAT3", "IRF1", "IRF7", "SPI1", "CEBPB", "JUN", "FOS", "JUNB", "ATF3", "MYC", "FOXP3", "TCF7")
focus <- focus[focus %in% tf_emt$source]
cat("== E-MTAB TF activity by stage (focus TFs) ==\n")
for (tf in focus) {
  v <- tf_emt[tf_emt$source == tf, ]
  m <- tapply(v$score, v$stage, mean)
  pv <- kruskal.test(score ~ stage, v)$p.value
  cat(sprintf("%-8s H=%5.2f E=%5.2f A=%5.2f  Kruskal p=%.3g\n", tf, m[1], m[2], m[3], pv))
}

fdf <- tf_emt[tf_emt$source %in% focus, ]
fdf$source <- factor(fdf$source, levels = focus)
pp2 <- ggplot(fdf, aes(x = stage, y = score, fill = stage)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  facet_wrap(~source, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(Healthy = "#4C72B0", Early = "#DD8452", Advanced = "#C44E52")) +
  labs(title = "TF activity across disease stages (E-MTAB-14896, DoRothEA A+B, ULM)", x = NULL, y = "ULM activity") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", legend.title = element_blank())
ggsave(file.path(sc_dir, "results", "EMTAB_TF_activity_stages.png"), pp2, width = 9, height = 6.5, dpi = 160)

cat("\nDONE\n")
