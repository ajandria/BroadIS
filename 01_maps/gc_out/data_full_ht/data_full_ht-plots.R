
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Setup -------------------------------------------------------------------
maps_table <- read_csv('01_maps/gc_out/data_full_ht_31_Oct_22_v2/02a_f_maps_table.csv') 

ht_syn_ps <- read_csv('01_maps/gc_out/data_full_ht_31_Oct_22_v2/ht_syn_ps.csv')

# Load --------------------------------------------------------------------
# MAPS bar plot
maps_full <- ggplot(maps_table,
       aes(x = consequence, y = maps, fill = consequence)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggthemes::theme_hc() +
  ggtitle('MAPS per variant functional class - based on variants with mu_snp defined and N_variants > 5000 full gnomAD.v3') +
  coord_flip()

maps_full %>% 
  plotly::ggplotly()

maps_sub <- ggplot(maps_table %>% filter(consequence %in% c('LC', 'HC', 'missense_variant', 'synonymous_variant')),
                    aes(x = consequence, y = maps, fill = consequence)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggthemes::theme_hc() +
  ggtitle('Subset of MAPS per variant functional class - based on variants with mu_snp defined and N_variants > 5000 full gnomAD.v3') +
  coord_flip()

maps_sub %>% 
  plotly::ggplotly()

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
  ggtitle('Singleton proportion vs. mutational rates - Synonymous Variants\nbased on variants with mu_snp defined full gnomAD.v3')

mu_snp_vs_ps %>% 
  plotly::ggplotly()


ht_syn_ps_recoded_cpg <- ht_syn_ps_recoded %>% 
  filter(variant_type == 'Transversion')

ht_syn_ps_recoded_non_cpg <- ht_syn_ps_recoded %>% 
  filter(variant_type == 'Non-CpG Transition')








