# 30_WGCNA.R : signed co-expression network of GSE77459 and module-trait association
suppressMessages({
  library(WGCNA); library(dplyr); library(ggplot2); library(pheatmap)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 2)

res_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"
deep_dir <- file.path(res_dir, "deep")

expr <- read.csv(file.path(res_dir, "GSE77459_expr_gene.csv"), row.names = 1, check.names = FALSE)
grp <- read.csv(file.path(res_dir, "GSE77459_groups.csv"))
datExpr <- as.data.frame(t(as.matrix(expr)))
datExpr <- datExpr[grp$sample, , drop = FALSE]
cat("Samples:", nrow(datExpr), " Genes:", ncol(datExpr), "\n")

# ---- filter genes: expressed and most variable ----
keep <- colMeans(datExpr) > 2.5
cat("Expressed genes:", sum(keep), "\n")
mad <- apply(datExpr[, keep, drop = FALSE], 2, mad)
datExpr <- datExpr[, keep, drop = FALSE][, order(-mad)][, 1:min(6000, sum(keep))]
cat("Genes for network:", ncol(datExpr), "\n")

gsg <- goodSamplesGenes(datExpr, verbose = 0)
cat("goodSamplesGenes allOK:", gsg$allOK, "\n")
if (!gsg$allOK) datExpr <- datExpr[, gsg$goodGenes]

# ---- sample clustering ----
trait <- data.frame(Pulpitis = as.integer(grp$group == "Pulpitis"), row.names = grp$sample)
sampleTree <- hclust(dist(datExpr), method = "average")
png(file.path(deep_dir, "WGCNA_sample_tree.png"), width = 1400, height = 900, res = 200)
plot(sampleTree, main = "Sample clustering", sub = "", xlab = "", cex.lab = 1, cex.axis = 1)
dev.off()

# ---- soft threshold ----
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", RsquaredCut = 0.90,
                        verbose = 0)
write.csv(sft$fitIndices, file.path(deep_dir, "WGCNA_softthreshold_fit.csv"), row.names = FALSE)
png(file.path(deep_dir, "WGCNA_softthreshold.png"), width = 1600, height = 800, res = 200)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft threshold (power)", ylab = "Scale-free topology R^2 (signed)",
     type = "n", main = "Scale independence"); text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, cex = 0.9)
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], xlab = "Soft threshold (power)",
     ylab = "Mean connectivity", type = "n", main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9)
dev.off()

power <- sft$powerEstimate
if (is.na(power)) power <- 12
cat("Selected power:", power, "\n")

# ---- network ----
net <- blockwiseModules(datExpr, power = power, networkType = "signed", TOMType = "signed",
                        minModuleSize = 30, maxBlockSize = ncol(datExpr),
                        mergeCutHeight = 0.30, numericLabels = TRUE, pamRespectsDendro = FALSE,
                        reassignThreshold = 0, saveTOMs = FALSE, deepSplit = 3, verbose = 0)
MEs <- net$MEs
mods <- net$colors
cat("Modules:", length(unique(mods)), "\n")

mod_labels <- paste0("M", as.integer(names(table(mods))))
names(mod_labels) <- names(table(mods))
mods_named <- mod_labels[as.character(mods)]

mod_df <- data.frame(gene = colnames(datExpr), module = mods_named, kME = NA)
for (m in unique(mods_named)) {
  idx <- mods_named == m
  MEcol <- paste0("ME", sub("M", "", m))
  mod_df$kME[idx] <- cor(datExpr[, idx, drop = FALSE], MEs[, MEcol])[, 1]
}
mod_df <- mod_df[order(mod_df$module, -mod_df$kME), ]
write.csv(mod_df, file.path(deep_dir, "WGCNA_module_assignment.csv"), row.names = FALSE)

# ---- module-trait ----
ME_cols <- MEs[, paste0("ME", sub("M", "", unique(mods_named))), drop = FALSE]
colnames(ME_cols) <- unique(mods_named)
mt <- data.frame(
  module = colnames(ME_cols),
  r = as.numeric(cor(ME_cols, trait$Pulpitis, use = "p")),
  p = as.numeric(corPvalueStudent(cor(ME_cols, trait$Pulpitis, use = "p"), nrow(datExpr)))
)
mt <- mt[order(-abs(mt$r)), ]
write.csv(mt, file.path(deep_dir, "WGCNA_module_trait.csv"), row.names = FALSE)
print(mt, row.names = FALSE, digits = 3)

png(file.path(deep_dir, "WGCNA_module_trait_heatmap.png"), width = 1500, height = 1300, res = 200)
pheatmap(t(matrix(mt$r, dimnames = list(mt$module, "Pulpitis"))), cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#2980b9", "white", "#c0392b"))(50),
         display_numbers = TRUE, number_format = "%.2f", fontsize = 9,
         main = "Module eigengene - pulpitis correlation")
dev.off()

# ---- hub genes ----
hub <- mod_df %>% group_by(module) %>% slice_max(kME, n = 25) %>% ungroup()
write.csv(hub, file.path(deep_dir, "WGCNA_hub_genes.csv"), row.names = FALSE)

# ---- dendrogram ----
png(file.path(deep_dir, "WGCNA_dendrogram.png"), width = 1800, height = 1000, res = 200)
plotDendroAndColors(net$dendrograms[[1]], labels2colors(mods[net$blockGenes[[1]]]),
                    "Module", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
                    guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

# ---- eigengene boxplot by group ----
top_mods <- head(mt$module, 5)
long <- do.call(rbind, lapply(top_mods, function(m) {
  data.frame(module = m, ME = ME_cols[, m], group = grp$group)
}))
p <- ggplot(long, aes(group, ME, fill = group)) + geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2) + facet_wrap(~ module, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(Normal = "#2980b9", Pulpitis = "#c0392b")) +
  labs(x = "", y = "Module eigengene", title = "Top pulpitis-associated module eigengenes") +
  theme_bw() + theme(legend.position = "none")
ggsave(file.path(deep_dir, "WGCNA_eigengene_box.png"), p, width = 10, height = 3.5, dpi = 300)
cat("DONE\n")
