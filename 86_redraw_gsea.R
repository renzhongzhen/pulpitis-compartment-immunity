# 86_redraw_gsea.R : 只重绘 GSEA 条形图（不重算 GSEA，规避 msigdbr 需要联网下载）
suppressMessages({ library(ggplot2) })

deep_dir <- 'C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/results/deep'

# ---- Hallmark ----
dh <- read.csv(file.path(deep_dir, 'GSEA_Hallmark.csv'), stringsAsFactors = FALSE)
dh <- dh[order(-dh$NES), ]
top_h <- rbind(head(dh, 8), tail(dh, 8))
top_h$label <- gsub('^HALLMARK_', '', top_h$ID)
top_h$label <- factor(top_h$label, levels = top_h$label[order(top_h$NES)])
p1 <- ggplot(top_h, aes(label, NES, fill = NES > 0)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c('TRUE' = '#c0392b', 'FALSE' = '#2980b9'),
                    breaks = c('TRUE', 'FALSE'),
                    labels = c('Activated in pulpitis', 'Suppressed in pulpitis')) +
  labs(x = '', y = 'Normalized enrichment score (NES)', fill = '',
       title = 'Hallmark gene sets in irreversible pulpitis (GSEA)') +
  theme_bw() + theme(axis.text.y = element_text(size = 7.5), legend.position = 'bottom')
ggsave(file.path(deep_dir, 'GSEA_Hallmark_bar.png'), p1, width = 8, height = 6.5, dpi = 300)
cat('SAVED Hallmark bar\n')

# ---- GO BP ----
dg <- read.csv(file.path(deep_dir, 'GSEA_GOBP.csv'), stringsAsFactors = FALSE)
dg <- dg[dg$NES > -Inf, ]
dg <- dg[order(-dg$NES), ]
top_g <- rbind(head(dg, 10), tail(dg, 10))
top_g$label <- substr(top_g$Description, 1, 55)
top_g$label <- factor(top_g$label, levels = top_g$label[order(top_g$NES)])
p2 <- ggplot(top_g, aes(label, NES, fill = NES > 0)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c('TRUE' = '#c0392b', 'FALSE' = '#2980b9'),
                    breaks = c('TRUE', 'FALSE'),
                    labels = c('Activated in pulpitis', 'Suppressed in pulpitis')) +
  labs(x = '', y = 'Normalized enrichment score (NES)', fill = '',
       title = 'GO biological processes in irreversible pulpitis (GSEA)') +
  theme_bw() + theme(axis.text.y = element_text(size = 6.5), legend.position = 'bottom')
ggsave(file.path(deep_dir, 'GSEA_GOBP_bar.png'), p2, width = 8.5, height = 7, dpi = 300)
cat('SAVED GOBP bar\n')

# 自检：正 NES 是否为红
cat('正 NES 条目数:', sum(top_h$NES > 0), ' 负 NES 条目数:', sum(top_h$NES < 0), '\n')
cat('DONE\n')
