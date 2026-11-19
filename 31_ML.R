# 31_ML.R : LASSO / random forest / SVM-RFE feature selection + diagnostic model with CV
suppressMessages({
  library(glmnet); library(randomForest); library(e1071); library(pROC)
  library(ggplot2); library(dplyr); library(tidyr)
})
set.seed(20261118)

res_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"
deep_dir <- file.path(res_dir, "deep")

expr <- read.csv(file.path(res_dir, "GSE77459_expr_gene.csv"), row.names = 1, check.names = FALSE)
grp <- read.csv(file.path(res_dir, "GSE77459_groups.csv"))
deg <- read.csv(file.path(res_dir, "GSE77459_DEG_full.csv"))
deg <- deg[!is.na(deg$gene) & !duplicated(deg$gene), ]
feat <- deg$gene[deg$change != "NS"]
cat("Candidate features:", length(feat), "\n")

ord <- grp$sample
X <- t(as.matrix(expr[feat, ord]))
X <- scale(X)
y <- as.integer(grp$group == "Pulpitis")

# ---------- 1. LASSO ----------
cvfit <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = nrow(X), type.measure = "auc")
coefs <- coef(cvfit, s = "lambda.min")
lasso_genes <- rownames(coefs)[as.numeric(coefs) != 0]
lasso_genes <- setdiff(lasso_genes, "(Intercept)")
cat("LASSO genes:", length(lasso_genes), "\n")

# ---------- 2. random forest ----------
rf <- randomForest(x = X, y = factor(y), ntree = 1000, importance = TRUE)
imp <- importance(rf)[, "MeanDecreaseGini"]
rf_genes <- names(sort(imp, decreasing = TRUE))[1:30]

# ---------- 3. SVM-RFE ----------
svm_rfe <- function(x, yv, n_feat = 30, step = 0.2) {
  f <- colnames(x)
  repeat {
    if (length(f) <= n_feat) break
    m <- svm(x[, f, drop = FALSE], yv, kernel = "linear", cost = 1, scale = FALSE)
    w <- drop(t(m$coefs) %*% m$SV)
    names(w) <- colnames(m$SV)
    drop_n <- max(1, floor(length(f) * step))
    f <- names(sort(abs(w), decreasing = TRUE))[1:(length(f) - drop_n)]
  }
  f
}
rfe_genes <- svm_rfe(X, y, n_feat = 30)
cat("RF top30 / RFE top30 done\n")

# ---------- consensus ----------
votes <- data.frame(gene = unique(c(lasso_genes, rf_genes, rfe_genes)))
votes$LASSO <- votes$gene %in% lasso_genes
votes$RF <- votes$gene %in% rf_genes
votes$SVM_RFE <- votes$gene %in% rfe_genes
votes$n_methods <- rowSums(votes[, c("LASSO", "RF", "SVM_RFE")])
write.csv(votes, file.path(deep_dir, "ML_feature_votes.csv"), row.names = FALSE)

cons3 <- votes$gene[votes$n_methods == 3]
cons2 <- votes$gene[votes$n_methods >= 2]
cat("Consensus (3 methods):", length(cons3), " (>=2):", length(cons2), "\n")
final_feats <- if (length(cons3) >= 2) cons3 else if (length(cons2) >= 2) cons2 else lasso_genes[1:min(3, length(lasso_genes))]
cat("Final features:", paste(final_feats, collapse = ", "), "\n")

# ---------- model: LOOCV AUC ----------
auc_loo <- NA; roc_loo <- NULL
if (length(final_feats) >= 1) {
  pred <- numeric(nrow(X))
  for (i in seq_len(nrow(X))) {
    tr <- setdiff(seq_len(nrow(X)), i)
    df <- data.frame(X[tr, final_feats, drop = FALSE], y = y[tr])
    m <- glm(y ~ ., data = df, family = binomial)
    nd <- data.frame(X[i, final_feats, drop = FALSE])
    pred[i] <- predict(m, newdata = nd, type = "response")
  }
  roc_loo <- roc(y, pred, quiet = TRUE)
  auc_loo <- as.numeric(auc(roc_loo))
  cat("LOOCV AUC:", round(auc_loo, 3), "\n")
  p_loo <- ggplot(as.data.frame(roc_loo$specificities), aes(x = 1 - roc_loo$specificities, y = roc_loo$sensitivities)) +
    geom_line(color = "#c0392b", linewidth = 1) + geom_abline(linetype = "dashed", color = "grey50") +
    coord_equal() + labs(x = "1 - specificity", y = "Sensitivity",
                         title = sprintf("LOOCV ROC (AUC = %.3f), GSE77459", auc_loo)) + theme_bw()
  ggsave(file.path(deep_dir, "ML_LOOCV_ROC.png"), p_loo, width = 5.5, height = 5, dpi = 300)
}

# ---------- repeated CV of the whole pipeline ----------
set.seed(42)
n_rep <- 50; n_fold <- 5
freq <- setNames(integer(length(deg$gene[deg$change != "NS"])), deg$gene[deg$change != "NS"])
auc_cv <- numeric(0)
for (r in 1:n_rep) {
  folds <- sample(rep(1:n_fold, length.out = nrow(X)))
  for (k in 1:n_fold) {
    tr <- which(folds != k); te <- which(folds == k)
    if (length(unique(y[te])) < 2) next
    cvf <- tryCatch(cv.glmnet(X[tr, , drop = FALSE], y[tr], family = "binomial", alpha = 1,
                              nfolds = 4, type.measure = "auc"), error = function(e) NULL)
    if (is.null(cvf)) next
    cc <- coef(cvf, s = "lambda.min")
    g <- setdiff(rownames(cc)[as.numeric(cc) != 0], "(Intercept)")
    if (length(g) < 1) next
    freq[g] <- freq[g] + 1
    df <- data.frame(X[tr, g, drop = FALSE], y = y[tr])
    if (length(unique(y[tr])) < 2) next
    m <- tryCatch(glm(y ~ ., data = df, family = binomial), error = function(e) NULL)
    if (is.null(m)) next
    pr <- tryCatch(predict(m, newdata = data.frame(X[te, g, drop = FALSE]), type = "response"),
                   error = function(e) NULL)
    if (!is.null(pr) && length(unique(y[te])) == 2) auc_cv <- c(auc_cv, as.numeric(auc(y[te], pr)))
  }
}
cat("CV folds evaluated:", length(auc_cv), " mean AUC:", round(mean(auc_cv, na.rm = TRUE), 3),
    " sd:", round(sd(auc_cv, na.rm = TRUE), 3), "\n")
sel_freq <- sort(freq, decreasing = TRUE)[1:25]
stab <- data.frame(gene = names(sel_freq), n_selected = as.integer(sel_freq),
                   pct = round(100 * as.integer(sel_freq) / (n_rep * n_fold), 1))
write.csv(stab, file.path(deep_dir, "ML_stability_selection.csv"), row.names = FALSE)

# ---------- figures ----------
pv <- votes[votes$n_methods >= 2, ]
pv$label <- pv$gene
pv <- pv %>% arrange(desc(n_methods))
pv$label <- factor(pv$label, levels = rev(pv$label))
p1 <- ggplot(pv, aes(label, n_methods, fill = factor(n_methods))) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("2" = "#f39c12", "3" = "#c0392b")) +
  labs(x = "", y = "Number of methods selecting the gene", fill = "Methods") +
  theme_bw() + theme(axis.text.y = element_text(size = 7))
ggsave(file.path(deep_dir, "ML_feature_votes.png"), p1, width = 6, height = 8, dpi = 300)

st <- stab
st$label <- factor(st$gene, levels = st$gene[order(st$pct)])
p2 <- ggplot(st, aes(label, pct)) + geom_col(fill = "#2c3e50", width = 0.7) + coord_flip() +
  labs(x = "", y = "Selection frequency across repeated CV (%)",
       title = "LASSO stability selection (50 x 5-fold)") + theme_bw() +
  theme(axis.text.y = element_text(size = 7))
ggsave(file.path(deep_dir, "ML_stability.png"), p2, width = 6, height = 7, dpi = 300)

saveRDS(list(lasso = lasso_genes, rf = rf_genes, rfe = rfe_genes, final = final_feats,
             auc_loo = auc_loo, auc_cv = auc_cv), file.path(deep_dir, "ML_results.rds"))
cat("DONE\n")
