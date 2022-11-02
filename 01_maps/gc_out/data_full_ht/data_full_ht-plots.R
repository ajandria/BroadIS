
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Setup -------------------------------------------------------------------
maps_table <- read_csv('01_maps/gc_out/data_full_ht_28_Oct_22/02a_f_maps_table.csv') 

ht_syn_ps <- read_csv('01_maps/gc_out/data_full_ht_31_Oct_22_v2/ht_syn_ps.csv')

# Load --------------------------------------------------------------------
# MAPS bar plot
maps_full <- ggplot(maps_table,
       aes(x = lof_csq_collapsed, y = maps, fill = lof_csq_collapsed)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggthemes::theme_hc() +
  ggtitle('MAPS per variant functional class') +
  coord_flip() +
  ylab("MAPS") +
  xlab('Collapsed lof and consequence calls')

maps_full %>% 
  plotly::ggplotly()

maps_sub <- ggplot(maps_table %>% filter(lof_csq_collapsed %in% c('LC', 'HC', 'missense_variant', 'synonymous_variant')),
                    aes(x = lof_csq_collapsed, y = maps, fill = lof_csq_collapsed)) +
  geom_bar(stat = 'identity', width = 0.3) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggthemes::theme_hc() +
  ggtitle('MAPS per variant functional class') +
  coord_flip() +
  ylab("MAPS") +
  xlab('Collapsed lof and consequence calls')

maps_sub %>% 
  plotly::ggplotly()

# Histogram
factors <- maps_table %>% arrange(N_variants)
maps_table$lof_csq_collapsed <- factor(maps_table$lof_csq_collapsed, levels = factors$lof_csq_collapsed)
ggplot(maps_table,
       aes(lof_csq_collapsed, N_variants)) +
  geom_bar(stat = 'identity') + 
  coord_flip() +
  ggthemes::theme_hc() +
  geom_text(aes(label=N_variants), position=position_dodge(width=0.9), hjust=-0.05) +
  ggtitle('Number of variants in each consequence')

# Recode into variant type
ht_syn_ps_recoded <- ht_syn_ps %>% 
  mutate(variant_type = case_when(
    (ref == 'G' & alt == 'A' & methylation_level == 2) | (ref == 'C' & alt == 'T' & methylation_level == 2) ~ 'CpG Transition',
    (ref == 'G' & alt == 'A' & methylation_level == 0) | (ref == 'C' & alt == 'T' & methylation_level == 0) ~ 'Non-CpG Transition',
    TRUE ~ 'Transversion'
  ))

# mu_snp vs ps 
mu_snp_vs_ps <- ht_syn_ps_recoded %>%
  ggplot(aes(x = mu_snp, y = ps, col = variant_type)) +
  geom_point() +
  ggthemes::theme_hc() +
  ggtitle('Singleton proportion vs. mutational rates - Synonymous Variants')

mu_snp_vs_ps %>% 
  plotly::ggplotly()


ht_syn_ps_recoded_cpg <- ht_syn_ps_recoded %>% 
  filter(variant_type == 'Transversion')

ht_syn_ps_recoded_non_cpg <- ht_syn_ps_recoded %>% 
  filter(variant_type == 'Non-CpG Transition')








