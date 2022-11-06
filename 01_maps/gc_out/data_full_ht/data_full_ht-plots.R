
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Setup -------------------------------------------------------------------
maps_table <- read_csv('01_maps/gc_out/data_full_ht_02_Nov_22_v4/02a_f_maps_table.csv') 
maps_table <- maps_table[order(maps_table$lof_csq_collapsed),]

ht_syn_ps <- read_csv('01_maps/gc_out/data_full_ht_02_Nov_22_v4/ht_syn_ps.csv')
ht_syn_ps <- ht_syn_ps[order(ht_syn_ps$context),]

# Load --------------------------------------------------------------------
# MAPS bar plot
color_table <- data.frame(
  consq = maps_table$lof_csq_collapsed,
  colors = c(
    '#AAAAAA', '#AAAAAA', '#AAAAAA', '#AAAAAA', '#9D1309', '#AAAAAA', '#AAAAAA',
    '#EE799F', '#FF6103', '#AAAAAA', '#AAAAAA', '#AAAAAA', '#AAAAAA', '#AAAAAA',
    '#AAAAAA', '#AAAAAA', '#AAAAAA', '#AAAAAA'
  )
)


(maps_full <- ggplot(maps_table,
       aes(x = lof_csq_collapsed, y = maps, fill = color_table$colors)) +
  geom_bar(stat = 'identity', width = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ggthemes::theme_hc() +
  ggtitle('MAPS across all variants') +
  coord_flip() +
  ylab("MAPS") +
  xlab('Collapsed lof and consequence calls') + 
  theme(legend.position='none') + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.46)) +
  theme(axis.title = element_text(face="bold"),
        axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        title = element_text(face="bold"))
 )

maps_full %>% 
  plotly::ggplotly()

# Sub on 4
maps_sub <- maps_table %>% 
  filter(lof_csq_collapsed %in% c('LC', 'HC', 'missense_variant', 'synonymous_variant'))
maps_sub$lof_csq_collapsed <- fct_rev(factor(maps_sub$lof_csq_collapsed, levels = maps_sub$lof_csq_collapsed))

(maps_sub_p <- ggplot(maps_sub,
                    aes(x = lof_csq_collapsed, y = maps)) +
                        #fill = c('#9D1309', '#EE799F', '#FF6103', '#AAAAAA'))) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', fill = c('grey', 'orange', 'darkred', 'darkred'),
               col = c('grey', 'orange', 'darkred', 'darkred'), binwidth = 0.015) +
  ylim(c(-0.01, 0.185)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() + theme(panel.border=element_blank(), legend.key = element_blank()) +
  ggtitle('MAPS in LoF & Missense & Synonymous') +
  ylab("MAPS") +
  xlab('Collapsed lof and consequence calls') + 
  theme(legend.position='none') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title = element_text(face="bold", size = 18),
        axis.text.x = element_text(face="bold", colour = 'black', size = 15),
        axis.text.y = element_text(face="bold", colour = 'black', size = 15),
        title = element_text(face="bold", size = 20))
)

maps_sub %>% 
  plotly::ggplotly()

# Histogram
factors <- maps_table %>% arrange(N_variants)
maps_table$lof_csq_collapsed <- factor(maps_table$lof_csq_collapsed, levels = factors$lof_csq_collapsed)

ggplot(maps_table,
       aes(lof_csq_collapsed, N_variants)) +
  geom_bar(stat = 'identity', aes(fill = '')) + 
  coord_flip() +
  geom_text(aes(label=N_variants), position=position_dodge(width=0.9), hjust=-0.05) +
  ggtitle('Number of variants in each consequence') +
  theme_bw() +
  ylim(c(0, max(maps_table$N_variants)+40000000))

# Recode into variant type
ht_syn_ps_recoded <- ht_syn_ps %>% 
  mutate(variant_type = case_when(
    (ref == 'G' & alt == 'A' & methylation_level == 2) | (ref == 'C' & alt == 'T' & methylation_level == 2) ~ 'CpG Transition',
    (ref == 'G' & alt == 'A' & methylation_level == 0) | (ref == 'C' & alt == 'T' & methylation_level == 0) ~ 'Non-CpG Transition',
    TRUE ~ 'Transversion'
  ))

ht_syn_ps_recoded <- ht_syn_ps_recoded[order(ht_syn_ps_recoded$variant_type),]
# mu_snp vs ps 
ht_syn_ps_recoded$context <- factor(ht_syn_ps_recoded$context, levels = unique(ht_syn_ps_recoded$context))

(mu_snp_vs_ps <- ht_syn_ps_recoded %>%
  ggplot(aes(x = mu_snp, y = ps, color = variant_type)) +
  geom_point(size = 5, alpha = 0.6) +
  scale_color_manual(values = c("CpG Transition" = '#2E9FFE',
                                "Non-CpG Transition" = "#458B00",
                                "Transversion" = "#EA4444")) +
  theme_classic() +
  ggtitle('Singleton proportion vs. mutational rates - Synonymous Variants')) +
  xlab('Mutation rate') +
  ylab('Singleton proportion') + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title = element_text(face="bold", size = 18),
        axis.text.x = element_text(face="bold", colour = 'black', size = 15),
        axis.text.y = element_text(face="bold", colour = 'black', size = 15),
        title = element_text(face="bold", size = 20),
        legend.text=element_text(size=18))+ 
  labs(color='Variant Type') 
)

mu_snp_vs_ps %>% 
  plotly::ggplotly()


ht_syn_ps_recoded_cpg <- ht_syn_ps_recoded %>% 
  filter(variant_type == 'Transversion')

ht_syn_ps_recoded_non_cpg <- ht_syn_ps_recoded %>% 
  filter(variant_type == 'Non-CpG Transition')








