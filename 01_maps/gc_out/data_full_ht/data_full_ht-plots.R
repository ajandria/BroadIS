
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Setup -------------------------------------------------------------------
maps_table <- read_csv('01_maps/gc_out/data_full_ht/02a_f_maps_table.csv')

ht_syn_ps <- read_csv('01_maps/gc_out/data_full_ht/ht_syn_ps.csv')

ht_inframe <- read_csv('01_maps/gc_out/data_full_ht/ht_f_inframe_del.csv')

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

# mu_snp vs ps 
mu_snp_vs_ps <- ht_syn_ps %>%
  ggplot(aes(x = mu_snp, y = ps, col = context)) +
  geom_point() +
  ggthemes::theme_hc() +
  ggtitle('Singleton proportion vs. mutational rates - Synonymous Variants - full gnomAD.v3')

mu_snp_vs_ps %>% 
  plotly::ggplotly()
