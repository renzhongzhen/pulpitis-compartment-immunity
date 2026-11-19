suppressMessages({ library(readxl); library(dplyr) })

data_dir <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/data"
res_dir  <- "C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results"

# ---------- GSE193448: LPS+laser(PBM) vs LPS vs Control, Ct 值 ----------
x <- read_excel(file.path(data_dir, "GSE193448_GA_RT-PCR.xls"), sheet = 1, col_names = FALSE)
genes <- x[[4]][-c(1,2)]                       # Symbol 列（第4列），跳过表头2行
ct <- as.data.frame(x[-(1:2), 5:10])
# 列映射（经原文 PMC8802107 核对修正）：S1/2=Ctrl, S3/4=LPS, S5/6=LPS+laser
colnames(ct) <- c("Ctrl1","Ctrl2","LPS1","LPS2","LPSlaser1","LPSlaser2")
ct[] <- lapply(ct, as.numeric)
ok <- !is.na(genes) & genes != ""
ct <- ct[ok, ]
genes_ok <- genes[ok]
# 去掉板卡质控孔并去重
drop_ctrl <- toupper(genes_ok) %in% c("PPC", "RTC", "NAC", "GDC", "AGCC")
ct <- ct[!drop_ctrl, ]
genes_clean <- genes_ok[!drop_ctrl]
keep_first <- !duplicated(genes_clean)
ct <- ct[keep_first, ]
rownames(ct) <- genes_clean[keep_first]
cat("PBM panel genes:", nrow(ct), "\n")

mean_ct <- data.frame(
  Ctrl = rowMeans(ct[, c("Ctrl1","Ctrl2")], na.rm = TRUE),
  LPS_laser = rowMeans(ct[, c("LPSlaser1","LPSlaser2")], na.rm = TRUE),
  LPS = rowMeans(ct[, c("LPS1","LPS2")], na.rm = TRUE))
mean_ct$gene <- rownames(mean_ct)

# dCt: 正 = 表达下降；负 = 表达上升
mean_ct$PBM_dCt   <- mean_ct$LPS_laser - mean_ct$LPS        # PBM 效应（LPS背景下）
mean_ct$LPS_dCt   <- mean_ct$LPS - mean_ct$Ctrl             # LPS 效应
mean_ct$PBM_vs_Ctrl_dCt <- mean_ct$LPS_laser - mean_ct$Ctrl
write.csv(mean_ct, file.path(res_dir, "PBM_panel_dCt.csv"), row.names = FALSE)

# ---------- 与牙髓炎 DEG 交叉 ----------
deg_path <- file.path(res_dir, "combined_DEG_full.csv")
if (!file.exists(deg_path)) deg_path <- file.path(res_dir, "GSE77459_DEG_full.csv")
deg <- read.csv(deg_path)
cat("Using DEG file:", basename(deg_path), "\n")

up_pulp <- deg$gene[deg$change == "Up"]
dn_pulp <- deg$gene[deg$change == "Down"]

cross <- merge(mean_ct, data.frame(gene = deg$gene, logFC_pulp = deg$logFC,
                                   padj_pulp = deg$adj.P.Val, change = deg$change),
               by = "gene")

# PBM 反向调控（抑制牙髓炎上调基因 / 恢复下调基因）
cross$PBM_anti <- ifelse(cross$change == "Up" & cross$PBM_dCt > 0, "PBM suppresses (anti-inflam)",
                  ifelse(cross$change == "Down" & cross$PBM_dCt < 0, "PBM restores",
                  ifelse(cross$change == "Up" & cross$PBM_dCt < 0, "PBM further increases",
                  ifelse(cross$change == "Down" & cross$PBM_dCt > 0, "PBM further decreases",
                         "no PBM effect"))))
cross$PBM_anti[is.na(cross$PBM_anti)] <- "PBM not measured"
cross <- cross[order(cross$PBM_dCt), ]
write.csv(cross, file.path(res_dir, "PBM_pulpitis_crosscheck.csv"), row.names = FALSE)

cat("\nPBM 交叉验证汇总：\n")
print(table(cross$PBM_anti))
cat("\n--- PBM 抑制的牙髓炎上调基因（候选机制基因） ---\n")
print(cross[which(cross$PBM_anti == "PBM suppresses (anti-inflam)"),
                 c("gene","logFC_pulp","PBM_dCt")])
cat("\n--- PBM 恢复的牙髓炎下调基因 ---\n")
print(cross[which(cross$PBM_anti == "PBM restores"),
                 c("gene","logFC_pulp","PBM_dCt")])
cat("\nDONE\n")
