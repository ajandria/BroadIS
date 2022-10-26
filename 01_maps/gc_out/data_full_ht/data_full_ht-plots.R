
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Setup -------------------------------------------------------------------
maps_table <- read_csv('01_maps/gc_out/data_full_ht_26_Oct_22/02a_f_maps_table.csv')

ht_syn_ps <- read_csv('01_maps/gc_out/data_full_ht_26_Oct_22/ht_syn_ps.csv')

# Load --------------------------------------------------------------------
# MAPS bar plot
maps_full <- ggplot(maps_table,
       aes(x = consequence, y = maps, fill = consequence)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggthemes::theme_hc() +
  ggtitle('MAPS per variant functional class - full gnomAD.v3')

maps_full %>% 
  plotly::ggplotly()

# Recode into variant type
ht_syn_ps_recoded <- ht_syn_ps %>% 
  mutate(variant_type = 
           ifelse((ref == 'G' & alt == 'A' & methylation_level != 0) | (ref == 'C' & alt == 'T' & methylation_level != 0), 'CpG Transition',
                  ifelse((ref == 'G' & alt == 'A' & methylation_level == 0) | (ref == 'C' & alt == 'T' & methylation_level == 0), 'Non-CpG Transition', 'Transversion'))
           )

# mu_snp vs ps 
mu_snp_vs_ps <- ht_syn_ps %>%
  ggplot(aes(x = mu_snp, y = ps, col = context)) +
  geom_point() +
  ggthemes::theme_hc() +
  ggtitle('Singleton proportion vs. mutational rates - Synonymous Variants - full gnomAD.v3')

mu_snp_vs_ps %>% 
  plotly::ggplotly()

