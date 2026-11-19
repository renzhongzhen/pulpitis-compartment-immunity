suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(org.Hs.eg.db); library(AnnotationDbi)
})
sc_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc"
outdir <- file.path(sc_dir, "results")
load(file.path(sc_dir, "pkgs", "CellChat-master", "CellChat-master", "data", "CellChatDB.human.rda"))
db <- CellChatDB.human$interaction
sec <- db[db$annotation == "Secreted Signaling", ]
lr <- sec[!grepl("\\+", sec$ligand) & !grepl("\\+", sec$receptor), c("ligand", "receptor")]

p <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/ea960c355a014573b0e904afc349b633processedmatrix.rds"
x <- readRDS(p)
md <- x@meta.data
ct <- as.character(Idents(x))
cnt <- LayerData(x, assay = "RNA", layer = "counts")
ens <- sub("\\..*$", "", rownames(cnt))
sym <- mapIds(org.Hs.eg.db, keys = ens, keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
keep <- !is.na(sym) & !duplicated(sym)
cnt2 <- cnt[keep, , drop = FALSE]; rownames(cnt2) <- sym[keep]
cs <- colSums(cnt2); cs[cs == 0] <- 1
rm(cnt, x); gc()

genes_lr <- unique(c(lr$ligand, lr$receptor)); genes_lr <- genes_lr[genes_lr %in% rownames(cnt2)]
Esub <- cnt2[genes_lr, , drop = FALSE]; rm(cnt2); gc()
ligs <- unique(lr$ligand); recs <- unique(lr$receptor)
ligs <- ligs[ligs %in% rownames(Esub)]; recs <- recs[recs %in% rownames(Esub)]
lr <- lr[lr$ligand %in% ligs & lr$receptor %in% recs, ]
i1 <- match(lr$ligand, ligs); i2 <- match(lr$receptor, recs)
cts <- sort(unique(ct))
pat <- as.character(md$samples)
patients <- unique(pat)
cat("patients:", length(patients), "cell types:", length(cts), "LR pairs:", nrow(lr), "\n")

stage_of <- function(s) ifelse(grepl("^Healthy", s), "Healthy", ifelse(grepl("^Early", s), "Early", "Advanced"))
res <- list()
for (pj in patients) {
  idx_p <- which(pat == pj)
  sub <- Esub[, idx_p, drop = FALSE]
  ctp <- ct[idx_p]; csp <- cs[idx_p]
  mm <- sapply(cts, function(cc) {
    idx <- which(ctp == cc)
    if (length(idx) == 0) return(rep(0, nrow(sub)))
    m <- as.matrix(sub[, idx, drop = FALSE])
    rowMeans(log2(sweep(m, 2, csp[idx], "/") * 1e4 + 1))
  })
  L <- mm[ligs, , drop = FALSE]; R <- mm[recs, , drop = FALSE]
  s <- matrix(0, length(cts), length(cts), dimnames = list(cts, cts))
  for (k in seq_len(nrow(lr))) s <- s + outer(as.numeric(L[i1[k], ]), as.numeric(R[i2[k], ]))
  res[[pj]] <- s
  cat("patient", pj, "cells", length(idx_p), "total strength", round(sum(s), 1), "\n")
}

long <- do.call(rbind, lapply(patients, function(pj) {
  s <- res[[pj]]
  data.frame(patient = pj, stage = stage_of(pj), from = rep(cts, each = length(cts)),
             to = rep(cts, times = length(cts)), strength = as.vector(s), stringsAsFactors = FALSE)
}))
write.csv(long, file.path(outdir, "EMTAB_LR_patient_strength.csv"), row.names = FALSE)

share_of <- function(s, a, b) s[a, b] / sum(s)
axes <- list(
  "Fib1 -> Endo3" = c("Fibroblast & Stromal cells 1", "Endothelial cells 3"),
  "Fib1 -> Immune" = c("Fibroblast & Stromal cells 1", "Blood & Immune cells 1"),
  "Immune -> Endo3" = c("Blood & Immune cells 1", "Endothelial cells 3"),
  "Immune -> Fib1" = c("Blood & Immune cells 1", "Fibroblast & Stromal cells 1"),
  "Endo1 -> Immune" = c("Endothelial cells 1", "Blood & Immune cells 1"),
  "Immune -> Endo1" = c("Blood & Immune cells 1", "Endothelial cells 1"),
  "Endo3 -> Immune" = c("Endothelial cells 3", "Blood & Immune cells 1"),
  "Immune -> Glial2" = c("Blood & Immune cells 1", "Glial cells 2"))

set.seed(20261119)
B <- 10000
out <- do.call(rbind, lapply(names(axes), function(nm) {
  ab <- axes[[nm]]
  if (!all(ab %in% cts)) return(NULL)
  sh <- sapply(patients, function(pj) share_of(res[[pj]], ab[1], ab[2]))
  st <- stage_of(patients)
  df <- data.frame(patient = patients, stage = st, share = as.numeric(sh), stringsAsFactors = FALSE)
  write.csv(df, file.path(outdir, paste0("EMTAB_LR_patient_share_", gsub("[^A-Za-z0-9]", "_", nm), ".csv")), row.names = FALSE)
  obs <- mean(df$share[df$stage == "Advanced"]) - mean(df$share[df$stage == "Healthy"])
  kw <- kruskal.test(share ~ factor(stage, levels = c("Healthy", "Early", "Advanced")), data = df)$p.value
  perm <- replicate(B, {
    lab <- sample(st)
    mean(df$share[lab == "Advanced"]) - mean(df$share[lab == "Healthy"])
  })
  p_two <- (sum(abs(perm) >= abs(obs)) + 1) / (B + 1)
  p_one <- (sum(perm >= obs) + 1) / (B + 1)
  data.frame(axis = nm, healthy = mean(df$share[df$stage == "Healthy"]),
             early = mean(df$share[df$stage == "Early"]), advanced = mean(df$share[df$stage == "Advanced"]),
             delta = obs, p_kw = kw, p_perm_two = p_two, p_perm_one_sided = p_one)
}))
out$p_kw_fdr <- p.adjust(out$p_kw, method = "BH")
out$p_perm_fdr <- p.adjust(out$p_perm_two, method = "BH")
print(out, digits = 4)
write.csv(out, file.path(outdir, "EMTAB_LR_patient_permutation.csv"), row.names = FALSE)
cat("\nDONE\n")
