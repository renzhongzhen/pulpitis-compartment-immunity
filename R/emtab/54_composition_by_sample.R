suppressMessages({library(Seurat); library(ggplot2); library(dplyr)})

obj_path <- 'C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results/GSE185222_seurat_annotated.rds'
out <- 'C:/Users/19354/Doubao/chats/2026-11-18/new-chat/pulp_bioinfo/sc/results/GSE185222_composition_bar.png'

o <- readRDS(obj_path)
md <- o@meta.data

tab <- md %>%
  count(sample, celltype) %>%
  group_by(sample) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(condition = ifelse(sample %in% c('DTP_01', 'DTP_04'),
                            'Healthy / enamel caries', 'Deep caries'),
         condition = factor(condition, levels = c('Deep caries', 'Healthy / enamel caries')))

p <- ggplot(tab, aes(x = sample, y = pct, fill = celltype)) +
  geom_col(width = 0.72, colour = NA) +
  facet_wrap(~condition, nrow = 1, scales = 'free_x') +
  labs(x = NULL, y = 'Percent of cells (%)', fill = 'Cell type',
       title = 'Cell type composition per sample') +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = 'grey92', colour = 'black'),
        strip.text = element_text(size = 9, face = 'bold'),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 8.5, colour = 'black'),
        legend.key.size = unit(0.35, 'cm'),
        legend.text = element_text(size = 7.5))

ggsave(out, p, width = 7.4, height = 5.4, dpi = 300, bg = 'white')
cat('SAVED:', out, '\n')
print(tab %>% group_by(sample, condition) %>% summarise(ncell = sum(n), .groups = 'drop'))
