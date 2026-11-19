suppressMessages({
  library(GEOquery)
  library(Biobase)
  library(hugene21sttranscriptcluster.db)
  library(AnnotationDbi)
})

data_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/data"
res_dir  <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"

gse <- getGEO(filename = file.path(data_dir, "GSE77459_series_matrix.txt.gz"), getGPL = FALSE)

# 分组：source_name_ch1 "normal pulp" vs "inflammed pulp"
grp <- ifelse(grepl("inflam", pData(gse)$source_name_ch1, ignore.case = TRUE), "Pulpitis", "Normal")
names(grp) <- colnames(exprs(gse))
table(grp)

expr <- exprs(gse)                       # RMA-normalized, log2
probes <- rownames(expr)

# 探针 -> gene symbol
map <- select(hugene21sttranscriptcluster.db, keys = probes, keytype = "PROBEID", columns = "SYMBOL")
map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", ]
cat("Probes mapped to at least one symbol:", length(unique(map$PROBEID)), "/", length(probes), "\n")

# 每个探针保留第一个 symbol
map <- map[!duplicated(map$PROBEID), ]
rownames(map) <- map$PROBEID

# 只保留能映射的探针
keep <- intersect(probes, map$PROBEID)
expr <- expr[keep, , drop = FALSE]
gene <- map[keep, "SYMBOL"]

# 每个基因取平均表达最高的探针（代表探针）
rowmean <- rowMeans(expr)
df <- data.frame(gene = gene, mean = rowmean, stringsAsFactors = FALSE)
best <- tapply(seq_len(nrow(df)), df$gene, function(i) i[which.max(df$mean[i])])
best_idx <- unlist(best)
expr_gene <- expr[best_idx, , drop = FALSE]
rownames(expr_gene) <- names(best)

cat("Gene-level matrix:", nrow(expr_gene), "genes x", ncol(expr_gene), "samples\n")

# 保存基因水平矩阵 + 分组
write.csv(expr_gene, file.path(res_dir, "GSE77459_expr_gene.csv"), quote = TRUE)
write.csv(data.frame(sample = names(grp), group = grp), file.path(res_dir, "GSE77459_groups.csv"), row.names = FALSE)
save(expr_gene, grp, file = file.path(data_dir, "GSE77459_genelevel.RData"))
cat("Saved gene-level matrix and groups.\n")
