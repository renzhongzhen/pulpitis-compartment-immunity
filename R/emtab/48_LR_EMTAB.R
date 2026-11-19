suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(org.Hs.eg.db); library(AnnotationDbi)
  library(dplyr); library(ggplot2)
})
sc_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc"
outdir <- file.path(sc_dir, "results")
# CellChatDB human: secreted signaling, single-subunit ligand-receptor pairs (as in GSE185222 analysis)
load(file.path(sc_dir, "pkgs", "CellChat-master", "CellChat-master", "data", "CellChatDB.human.rda"))
db <- CellChatDB.human$interaction
sec <- db[db$annotation == "Secreted Signaling", ]
lr <- sec[!grepl("\\+", sec$ligand) & !grepl("\\+", sec$receptor), c("ligand", "receptor", "pathway_name")]
cat("secreted single-subunit LR pairs:", nrow(lr), "\n")

p <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/ea960c355a014573b0e904afc349b633processedmatrix.rds"
x <- readRDS(p)
md <- x@meta.data
ct <- factor(as.character(Idents(x)))
cnt <- LayerData(x, assay = "RNA", layer = "counts")
ens <- sub("\\..*$", "", rownames(cnt))
sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
keep <- !is.na(sym) & !duplicated(sym)
cnt2 <- cnt[keep, , drop = FALSE]; rownames(cnt2) <- sym[keep]
cs <- colSums(cnt2); cs[cs == 0] <- 1
rm(cnt, x); gc()

cts <- levels(ct)
stage_of_pat <- function(s) ifelse(grepl("^Healthy", s), "Healthy",
                            ifelse(grepl("^Early", s), "Early", "Advanced"))
stage <- stage_of_pat(as.character(md$samples))
stages <- c("Healthy", "Early", "Advanced")

genes_lr <- unique(c(lr$ligand, lr$receptor))
genes_lr <- genes_lr[genes_lr %in% rownames(cnt2)]
Esub <- cnt2[genes_lr, , drop = FALSE]
rm(cnt2); gc()
cat("LR genes retained:", nrow(Esub), "\n")

ligs <- unique(lr$ligand); recs <- unique(lr$receptor)
ligs <- ligs[ligs %in% rownames(Esub)]; recs <- recs[recs %in% rownames(Esub)]
cat("ligands present:", length(ligs), " receptors present:", length(recs), "\n")

mean_ct <- function(mat_sparse, cs, ct, st, stage_val, cts) {
  mm <- do.call(cbind, lapply(cts, function(cc) {
    idx <- which(ct == cc & st == stage_val)
    if (length(idx) == 0) return(rep(0, nrow(mat_sparse)))
    m <- as.matrix(mat_sparse[, idx, drop = FALSE])
    rowMeans(log2(sweep(m, 2, cs[idx], "/") * 1e4 + 1))
  }))
  colnames(mm) <- cts; mm
}

res_all <- list()
Lmean <- list(); Rmean <- list()
for (sv in stages) {
  mm <- mean_ct(Esub, cs, ct, stage, sv, cts)
  cat("stage", sv, "mm dim:", paste(dim(mm), collapse = "x"), "\n")
  Lmean[[sv]] <- mm[ligs, , drop = FALSE]
  Rmean[[sv]] <- mm[recs, , drop = FALSE]
  s <- matrix(0, length(cts), length(cts), dimnames = list(cts, cts))
  n_bad <- 0
  for (i in seq_len(nrow(lr))) {
    i1 <- match(lr$ligand[i], ligs); i2 <- match(lr$receptor[i], recs)
    if (is.na(i1) || is.na(i2)) next
    a <- outer(as.numeric(Lmean[[sv]][i1, ]), as.numeric(Rmean[[sv]][i2, ]))
    if (nrow(a) != length(cts) || ncol(a) != length(cts)) {
      n_bad <- n_bad + 1
      next
    }
    s <- s + a
  }
  if (n_bad > 0) cat("  skipped", n_bad, "non-conformable pairs\n")
  res_all[[sv]] <- s
}
tot <- sapply(res_all, sum)
cat("total communication strength by stage:", paste(names(tot), round(tot, 1), collapse = ", "), "\n")
res_norm <- lapply(res_all, function(m) m / sum(m))

# directional stage change per celltype pair (normalized share: Advanced minus Healthy)
pairs <- expand.grid(from = cts, to = cts, stringsAsFactors = FALSE)
delta <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  a <- pairs$from[i]; b <- pairs$to[i]
  data.frame(from = a, to = b,
             healthy = res_norm[["Healthy"]][a, b],
             early = res_norm[["Early"]][a, b],
             advanced = res_norm[["Advanced"]][a, b],
             delta_adv = res_norm[["Advanced"]][a, b] - res_norm[["Healthy"]][a, b])
}))
delta <- delta[order(-delta$delta_adv), ]
write.csv(delta, file.path(outdir, "EMTAB_LR_pair_delta.csv"), row.names = FALSE)
cat("== top 15 sender->receiver pairs by Advanced minus Healthy strength ==\n")
print(head(delta, 15), row.names = FALSE)

# focus: chemokine axis involving Endothelial cells 3 (ACKR1+ venule) and immune cells
imm <- grep("Blood & Immune", cts, value = TRUE)
endo3 <- "Endothelial cells 3"
chemo <- lr[grepl("CCL|CXCL|CCR|CXCR|ACKR|XCL", lr$ligand) | grepl("CCL|CXCL|CCR|CXCR|ACKR|XCL", lr$receptor), ]
cat("chemokine-related pairs:", nrow(chemo), "\n")
ch_delta <- do.call(rbind, lapply(seq_len(nrow(chemo)), function(i) {
  li <- chemo$ligand[i]; re <- chemo$receptor[i]
  i1 <- match(li, ligs); i2 <- match(re, recs)
  if (is.na(i1) || is.na(i2)) return(NULL)
  v <- sapply(stages, function(sv) {
    s <- outer(as.numeric(Lmean[[sv]][i1, ]), as.numeric(Rmean[[sv]][i2, ]))
    dimnames(s) <- list(cts, cts)
    s <- s / sum(s)
    data.frame(endo3_rec = s[endo3, ], imm_send = max(s[imm, endo3]), endo3_send = max(s[endo3, imm]))
  }, simplify = FALSE)
  vv <- do.call(rbind, v); vv$stage <- stages
  d <- data.frame(ligand = li, receptor = re, pathway = chemo$pathway_name[i],
                  endo3_total = sum(vv$endo3_rec), imm_to_endo3_adv = vv$imm_send[3], imm_to_endo3_h = vv$imm_send[1],
                  endo3_to_imm_adv = vv$endo3_send[3], endo3_to_imm_h = vv$endo3_send[1])
  d$d_imm2endo <- d$imm_to_endo3_adv - d$imm_to_endo3_h
  d$d_endo2imm <- d$endo3_to_imm_adv - d$endo3_to_imm_h
  d
}))
ch_delta <- ch_delta[order(-ch_delta$d_imm2endo), ]
write.csv(ch_delta, file.path(outdir, "EMTAB_LR_chemokine_endo3.csv"), row.names = FALSE)
cat("== chemokine interactions with Endothelial 3 (immune->Endo3 gain) ==\n")
print(head(ch_delta[, c("ligand", "receptor", "pathway", "d_imm2endo", "d_endo2imm")], 12), row.names = FALSE)

top <- head(delta[order(-delta$delta_adv), ], 15)
top$pair <- paste(top$from, "->", top$to)
top$pair <- factor(top$pair, levels = rev(top$pair))
pp <- ggplot(top, aes(x = delta_adv, y = pair)) +
  geom_col(fill = "#D73027", width = 0.7) +
  labs(title = "Top cell-cell pairs gaining ligand-receptor strength (Advanced minus Healthy, E-MTAB-14896)",
       x = "communication strength gain", y = NULL) +
  theme_bw(base_size = 9)
ggsave(file.path(outdir, "EMTAB_LR_pair_gain.png"), pp, width = 9, height = 4.8, dpi = 160)
cat("\nDONE\n")
