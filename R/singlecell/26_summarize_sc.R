d <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results"
m <- read.csv(file.path(d, "GSE185222_mac_M1M2_scores.csv"))
cat("=== M1/M2 by condition ===\n")
print(aggregate(cbind(M1_score, M2_score) ~ condition, data = m, mean), digits = 3)
g <- split(m$M1_score - m$M2_score, m$condition)
cat("M1-M2 means:", paste(names(g), signif(sapply(g, mean), 3), collapse = " | "), "\n")
dc <- m$M1_score[m$condition == "Deep caries"]
hc <- m$M1_score[m$condition != "Deep caries"]
cat("p(M1 deep vs healthy):", wilcox.test(dc, hc)$p.value, "\n")

tj <- read.csv(file.path(d, "GSE185222_trajectory_pseudotime.csv"))
cat("\n=== pseudotime by condition ===\n")
print(aggregate(cbind(Lineage1, Lineage2) ~ condition, data = tj, mean, na.rm = TRUE), digits = 3)
cat("=== by celltype ===\n")
print(aggregate(cbind(Lineage1, Lineage2) ~ celltype, data = tj, mean, na.rm = TRUE), digits = 3)

lr <- read.csv(file.path(d, "GSE185222_LR_delta_deep_vs_healthy.csv"))
sel <- lr[grepl("SPP1|PDGF|VEGFA|TGFB1|IL1B|TNF", lr$key), ]
cat("\n=== selected axes (", nrow(sel), " rows) ===\n")
print(head(sel[order(-sel$delta), ], 15), row.names = FALSE)

kg <- read.csv(file.path(d, "GSE185222_keygene_celltype_avg.csv"))
cat("\n=== keygene_celltype_avg cols ===\n")
cat(paste(colnames(kg), collapse = ", "), "\n")
print(head(kg, 12), row.names = FALSE)
