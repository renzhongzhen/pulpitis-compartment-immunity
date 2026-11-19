# -*- coding: utf-8 -*-
# 23_communication_GSE185222.R : ligand-receptor communication scoring with CellChatDB
# custom implementation (no CellChat compilation needed): mean-expression co-expression
# score + label-permutation significance, computed separately for Healthy vs Deep caries
suppressMessages({
  library(Seurat); library(Matrix); library(dplyr)
})

out_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"
db_file <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/pkgs/CellChat-master/CellChat-master/data/CellChatDB.human.rda"
obj <- readRDS(file.path(out_dir, "GSE185222_seurat_annotated.rds"))

load(db_file)  # loads CellChatDB.human
inter <- CellChatDB.human$interaction
cat("Total LR interactions in CellChatDB.human:", nrow(inter), "\n")
# secreted signaling only, single-subunit ligand & receptor
inter <- inter[grepl("Secreted Signaling", inter$annotation), ]
inter <- inter[!grepl("_", inter$ligand) & !grepl("_", inter$receptor), ]
inter <- inter[, c("interaction_name", "pathway_name", "ligand", "receptor")]
inter <- inter[inter$ligand %in% rownames(obj) & inter$receptor %in% rownames(obj), ]
cat("Usable single-subunit secreted LR pairs:", nrow(inter), "\n")

lr_genes <- unique(c(inter$ligand, inter$receptor))
expr <- as.matrix(GetAssayData(obj, assay = "RNA", layer = "data")[lr_genes, ])
expr <- t(expr)  # cells x genes

group_score <- function(groups, expr_mat, inter, n_min = 10) {
  lev <- unique(groups)
  avg <- sapply(lev, function(g) {
    idx <- which(groups == g)
    if (length(idx) < n_min) return(rep(NA_real_, ncol(expr_mat)))
    colMeans(expr_mat[idx, , drop = FALSE])
  })
  avg <- t(avg)  # celltypes x genes
  rownames(avg) <- lev
  sc <- sapply(seq_len(nrow(inter)), function(i) {
    l <- avg[, inter$ligand[i]]; r <- avg[, inter$receptor[i]]
    sqrt(pmax(l, 0) * pmax(r, 0))
  })
  # sc: celltypes x pairs -> expand sender x receiver
  out <- list()
  for (i in seq_len(nrow(inter))) {
    s <- sc[, i]
    out[[i]] <- expand.grid(sender = lev, receiver = lev, stringsAsFactors = FALSE) %>%
      mutate(pair = inter$interaction_name[i], pathway = inter$pathway_name[i],
             score = as.numeric(outer(s, s, function(a, b) sqrt(pmax(a,0) * pmax(b,0)))[cbind(match(sender, lev), match(receiver, lev))]))
  }
  do.call(rbind, out)
}

# ---- observed scores per condition ----
res <- list()
for (cond in c("Healthy/Enamel caries", "Deep caries")) {
  sub <- subset(obj, condition == cond)
  grp <- as.character(sub$celltype)
  res[[cond]] <- group_score(grp, t(as.matrix(GetAssayData(sub, assay = "RNA", layer = "data")[lr_genes, ])), inter)
  cat(cond, "done\n")
}
scores <- bind_rows(lapply(names(res), function(n) { d <- res[[n]]; d$condition <- n; d }))
write.csv(scores, file.path(out_dir, "GSE185222_LR_scores_all.csv"), row.names = FALSE)

# ---- permutation test for Deep caries (top candidate pairs) ----
sub <- subset(obj, condition == "Deep caries")
grp <- as.character(sub$celltype)
expr_d <- t(as.matrix(GetAssayData(sub, assay = "RNA", layer = "data")[lr_genes, ]))
obs <- res[["Deep caries"]]
top <- obs %>% arrange(desc(score)) %>% slice(1:300)
set.seed(42)
nperm <- 200
perm_scores <- matrix(NA_real_, nrow = nrow(top), ncol = nperm)
for (p in seq_len(nperm)) {
  gl <- sample(grp)
  ps <- group_score(gl, expr_d, inter)
  ps <- ps %>% mutate(key = paste(sender, receiver, pair, sep = "|"))
  key <- paste(top$sender, top$receiver, top$pair, sep = "|")
  perm_scores[, p] <- ps$score[match(key, ps$key)]
  if (p %% 50 == 0) cat("perm", p, "\n")
}
top$p_value <- rowMeans(perm_scores >= top$score, na.rm = TRUE)
write.csv(top, file.path(out_dir, "GSE185222_LR_top_deep_perm.csv"), row.names = FALSE)

# ---- differential communication deep vs healthy ----
wide <- scores %>% mutate(key = paste(sender, receiver, pair, sep = "|")) %>%
  select(key, score, condition) %>% tidyr::pivot_wider(names_from = condition, values_from = score)
wide$delta <- wide[["Deep caries"]] - wide[["Healthy/Enamel caries"]]
wide <- wide %>% arrange(desc(delta))
write.csv(wide, file.path(out_dir, "GSE185222_LR_delta_deep_vs_healthy.csv"), row.names = FALSE)
cat("\nTop 15 up-communicated in deep caries:\n")
print(head(wide %>% filter(!is.na(delta)) %>% arrange(desc(delta)), 15))
cat("DONE\n")
