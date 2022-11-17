
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Setup -------------------------------------------------------------------
maps_table <- read_csv("01_maps/gc_out/data_full_ht_02_Nov_22_v4/02a_f_maps_table.csv")
maps_table <- maps_table[order(maps_table$lof_csq_collapsed), ]

ht_syn_ps <- read_csv("01_maps/gc_out/data_full_ht_02_Nov_22_v4/ht_syn_ps.csv")
ht_syn_ps <- ht_syn_ps[order(ht_syn_ps$context), ]

# Load --------------------------------------------------------------------
# MAPS bar plot
color_table <- data.frame(
  consq = maps_table$lof_csq_collapsed,
  colors = c(
    "#AAAAAA", "#AAAAAA", "#AAAAAA", "#AAAAAA", "#9D1309", "#AAAAAA", "#AAAAAA",
    "#EE799F", "#FF6103", "#AAAAAA", "#AAAAAA", "#AAAAAA", "#AAAAAA", "#AAAAAA",
    "#AAAAAA", "#AAAAAA", "#AAAAAA", "#AAAAAA"
  )
)


(maps_full <- ggplot(
  maps_table,
  aes(x = lof_csq_collapsed, y = maps, fill = color_table$colors)
) +
  geom_bar(stat = "identity", width = 0.7) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggthemes::theme_hc() +
  ggtitle("MAPS across all variants") +
  coord_flip() +
  ylab("MAPS") +
  xlab("Collapsed lof and consequence calls") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.46)) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    title = element_text(face = "bold")
  )
)

maps_full %>%
  plotly::ggplotly()

# Sub on 4
maps_sub <- maps_table %>%
  filter(lof_csq_collapsed %in% c("LC", "HC", "missense_variant", "synonymous_variant"))
maps_sub$lof_csq_collapsed <- fct_rev(factor(maps_sub$lof_csq_collapsed, levels = maps_sub$lof_csq_collapsed))

(maps_sub_p <- ggplot(
  maps_sub,
  aes(x = lof_csq_collapsed, y = maps)
) +
  # fill = c('#9D1309', '#EE799F', '#FF6103', '#AAAAAA'))) +
  geom_dotplot(
    binaxis = "y", stackdir = "center", fill = c("grey", "orange", "darkred", "darkred"),
    col = c("grey", "orange", "darkred", "darkred"), binwidth = 0.015
  ) +
  ylim(c(-0.01, 0.185)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme_classic() +
  theme(panel.border = element_blank(), legend.key = element_blank()) +
  ggtitle("MAPS in LoF & Missense & Synonymous") +
  ylab("MAPS") +
  xlab("Collapsed lof and consequence calls") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    title = element_text(face = "bold", size = 20)
  )
)

maps_sub %>%
  plotly::ggplotly()

# Histogram
factors <- maps_table %>% arrange(N_variants)
maps_table$lof_csq_collapsed <- factor(maps_table$lof_csq_collapsed, levels = factors$lof_csq_collapsed)

ggplot(
  maps_table,
  aes(lof_csq_collapsed, N_variants)
) +
  geom_bar(stat = "identity", aes(fill = "")) +
  coord_flip() +
  geom_text(aes(label = N_variants), position = position_dodge(width = 0.9), hjust = -0.05) +
  ggtitle("Number of variants in each consequence") +
  theme_bw() +
  ylim(c(0, max(maps_table$N_variants) + 40000000))

# Recode into variant type
ht_syn_ps_recoded <- ht_syn_ps %>%
  mutate(variant_type = case_when(
    (ref == "G" & alt == "A" & methylation_level == 2) | (ref == "C" & alt == "T" & methylation_level == 2) ~ "CpG Transition",
    (ref == "G" & alt == "A" & methylation_level == 0) | (ref == "C" & alt == "T" & methylation_level == 0) ~ "Non-CpG Transition",
    TRUE ~ "Transversion"
  ))

ht_syn_ps_recoded <- ht_syn_ps_recoded[order(ht_syn_ps_recoded$variant_type), ]
# mu_snp vs ps
ht_syn_ps_recoded$context <- factor(ht_syn_ps_recoded$context, levels = unique(ht_syn_ps_recoded$context))

(mu_snp_vs_ps <- ht_syn_ps_recoded %>%
  ggplot(aes(x = mu_snp, y = ps, color = variant_type) +
    geom_point(size = 5, alpha = 0.6) +
    scale_color_manual(values = c(
      "CpG Transition" = "#2E9FFE",
      "Non-CpG Transition" = "#458B00",
      "Transversion" = "#EA4444"
    )) +
    theme_classic() +
    ggtitle("Singleton proportion vs. mutational rates - Synonymous Variants")) +
  xlab("Mutation rate") +
  ylab("Singleton proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18)
  ))

mu_snp_vs_ps %>%
  plotly::ggplotly()


ht_syn_ps_recoded_cpg <- ht_syn_ps_recoded %>%
  filter(variant_type == "Transversion")

ht_syn_ps_recoded_non_cpg <- ht_syn_ps_recoded %>%
  filter(variant_type == "Non-CpG Transition")



# per gene ----------------------------------------------------------------
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(ggplot2)

# Setup -------------------------------------------------------------------
maps_table <- read_csv("01_maps/gc_out/data_full_ht_08_Nov_22_v6/02a_f_maps_table.csv")

ht_syn_ps <- read_csv("01_maps/gc_out/data_full_ht_08_Nov_22_v6/ht_syn_ps.csv")

# ps vs.  mu_snp ----------------------------------------------------------
ht_syn_ps_recoded <- ht_syn_ps %>%
  mutate(variant_type = case_when(
    (ref == "G" & alt == "A" & methylation_level == 2) | (ref == "C" & alt == "T" & methylation_level == 2) ~ "CpG Transition",
    (ref == "G" & alt == "A" & methylation_level == 0) | (ref == "C" & alt == "T" & methylation_level == 0) ~ "Non-CpG Transition",
    TRUE ~ "Transversion"
  ))

(mu_snp_vs_ps <- ht_syn_ps_recoded %>%
  ggplot(aes(x = mu_snp, y = ps, color = variant_type)) +
  geom_point(size = 1, alpha = 0.6) +
  scale_color_manual(values = c(
    "CpG Transition" = "#2E9FFE",
    "Non-CpG Transition" = "#458B00",
    "Transversion" = "#EA4444"
  )) +
  theme_classic() +
  ggtitle("Singleton proportion vs. mutational rates - Synonymous Variants") +
  xlab("Mutation rate") +
  ylab("Singleton proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18)
  ))




ht_syn_ps_recoded_cpg <- ht_syn_ps_recoded %>%
  filter(variant_type == "Transversion")

ht_syn_ps_recoded_non_cpg <- ht_syn_ps_recoded %>%
  filter(variant_type == "Non-CpG Transition")



# Hist --------------------------------------------------------------------
maps_table <- maps_table %>% filter(lof_csq_collapsed == "synonymous_variant")
factors <- maps_table %>% arrange(N_variants)
maps_table$gene_ids <- factor(maps_table$gene_ids, levels = factors$gene_ids)

tiff("syn_per_gene.tiff", 1400, 1080, units = "px")
ggplot(
  maps_table,
  aes(gene_ids, N_variants)
) +
  geom_bar(stat = "identity", width = 0.5, color = "navy") +
  coord_flip() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank()
  )
# geom_text(aes(label=N_variants), position=position_dodge(width=0.9), hjust=-0.05) +
ggtitle("Number of synonymous variants in each gene") +
  theme_bw()
dev.off()

ecdf <- ecdf(maps_table$N_variants)
plot(ecdf)
ggplot(maps_table, aes(x = gene_ids, y = cumsum(N_variants))) +
  geom_line() +
  geom_point() +
  xlim(0, 500)




ggplot(
  maps_table,
  aes(N_variants, maps, col = lof_csq_collapsed)
) +
  geom_point() +
  theme_bw() +
  ggtitle("MAPS vs. N_variants | Gene's Synonymous Variants Only")

maps_table_ecdf <- maps_table %>%
  mutate(ecdf = ecdf(N_variants))

maps_table_ecdf2 <- maps_table_ecdf %>%
  filter(N_variants <= 500)

ggplot(maps_table_ecdf2, aes(N_variants, 1 - ecdf)) +
  geom_step()



#
#
#
#
# maps_per_gene_08_Nov_22_v1 ----------------------------------------------

# Dependencies
library(tidyverse)
library(ggplot2)

setwd("/Users/adrian/BroadIS/01_maps")

# Setup
maps_table <- read_csv("gc_out/maps_per_gene_08_Nov_22_v1/02a_f_maps_table.csv")
ht_syn_ps <- read_csv("gc_out/maps_per_gene_08_Nov_22_v1/ht_syn_ps.csv")

# Plots

# I)
ggplot(
  maps_table,
  aes(maps, N_variants, col = lof_csq_collapsed)
) +
  geom_point(size = 2) +
  ggthemes::theme_hc() +
  facet_grid(~lof_csq_collapsed) +
  ggtitle("MAPS vs. N_variants Per Gene\nSingleton proportion adjusted using all synonymous variants grouped into contexts by methylation level")  +
  xlab('MAPS') +
  ylab('N Variants') +
  scale_y_continuous(limits = c(0, max(maps_table$N_variants)+1000), 
                     breaks = seq(0, max(maps_table$N_variants)+1000, 2000)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18)
  )


# II)
# Recode into variant type
ht_syn_ps_recoded <- ht_syn_ps %>%
  mutate(variant_type = case_when(
    (ref == "G" & alt == "A" & methylation_level == 2) | (ref == "C" & alt == "T" & methylation_level == 2) ~ "CpG Transition",
    (ref == "G" & alt == "A" & methylation_level == 0) | (ref == "C" & alt == "T" & methylation_level == 0) ~ "Non-CpG Transition",
    TRUE ~ "Transversion"
  ))

ht_syn_ps_recoded <- ht_syn_ps_recoded[order(ht_syn_ps_recoded$variant_type), ]
# mu_snp vs ps
ht_syn_ps_recoded$context <- factor(ht_syn_ps_recoded$context, levels = unique(ht_syn_ps_recoded$context))

(mu_snp_vs_ps <- ht_syn_ps_recoded %>%
  ggplot(aes(x = mu_snp, y = ps, color = variant_type)) +
  geom_point(size = 5, alpha = 0.6) +
  scale_color_manual(values = c(
    "CpG Transition" = "#2E9FFE",
    "Non-CpG Transition" = "#458B00",
    "Transversion" = "#EA4444"
  )) +
  theme_classic() +
  ggtitle("Singleton proportion vs. mutational rates - Synonymous Variants") +
  xlab("Mutation rate") +
  ylab("Singleton proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18)
  ))

# III) ecdf on MAPS in synonymous
maps_table_synonymous <- maps_table %>%
  filter(lof_csq_collapsed == 'synonymous_variant')
# make ecdf of abs(maps)
ecdf_abs_maps <- ecdf(abs(maps_table_synonymous$maps))
# annotate the ecdf to table using the ecdf of abs(maps)
maps_table_synonymous_ecdf <- maps_table_synonymous %>% 
  mutate(ecdf_abs_maps = ecdf_abs_maps(abs(maps)))
# plot
ggplot(maps_table_synonymous_ecdf, 
       aes(abs(maps), 1 - ecdf_abs_maps)) +
  geom_step() +
  ggthemes::theme_hc() +
  xlab('abs(MAPS)') +
  ylab('Gene Ratio') +
  ggtitle('Cumulative Plot of abs(MAPS) Per Gene Based on Synonymous Variants')
  
# IV) MAPS missense vs. HC
maps_table_not_lc <- maps_table %>% 
  filter(lof_csq_collapsed != 'LC') %>%
  select(gene_ids, lof_csq_collapsed, maps) %>% 
  pivot_wider(
    names_from = lof_csq_collapsed,
    values_from = maps
  )

ggplot(maps_table_not_lc,
       aes(missense_variant, HC, alpha=1-abs(synonymous_variant))) +
  geom_point() +
  theme_bw() +
  ggtitle('MAPS of HC vs. Missense Variants') +
  xlab('Missense MAPS') +
  ylab('HC Maps') +
  labs(alpha = "1-|MAPS Synonymous|")

# V) MAPS HC vs. oe_lof
lof_metrics <- readr::read_delim('/Users/adrian/BroadIS/01_maps/data/gnomad.v2.1.1.lof_metrics.by_gene.txt')

maps_table_binned <- maps_table %>% 
  select(gene_ids, lof_csq_collapsed, maps) %>% 
  filter(lof_csq_collapsed == 'synonymous_variant') %>% 
  mutate(syn_maps_binned = case_when(
    (abs(maps) <= 0.001) ~ "(0, 0.001]",
    (abs(maps) <= 0.01) ~ "(0.001, 0.01]",
    (abs(maps) <= 0.05) ~ "(0.01, 0.05]",
    (abs(maps) <= 1) ~ "(0.05, 1]",
    TRUE ~ "ERROR"
  )) %>% 
  select(-c(maps)) %>% 
  select(gene_ids, syn_maps_binned) %>% 
  rename(gene_id = gene_ids)


maps_table_metrics <- maps_table %>% 
  select(gene_ids, N_variants, lof_csq_collapsed, maps) %>% 
  rename(gene_id = gene_ids) %>% 
  left_join(lof_metrics %>% select(gene_id, oe_lof_upper)) %>% 
  filter(lof_csq_collapsed == 'HC')

maps_table_metrics_final <- maps_table_metrics %>% 
  left_join(maps_table_binned)
  

ggplot(maps_table_metrics_final,
       aes(x = oe_lof_upper, y = maps, col = syn_maps_binned)) +
  geom_point() +
  ggthemes::theme_hc() +
  ggtitle('Per Gene oe_lof_upper vs. MAPS based on HC variants') +
  ylab('MAPS') +
  xlab('OE_LOF_UPPER') +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18)
  )

# VI) MAPS HC vs. oe_lof
lof_metrics <- readr::read_delim('/Users/adrian/BroadIS/01_maps/data/gnomad.v2.1.1.lof_metrics.by_gene.txt')

maps_table_binned_2 <- maps_table %>% 
  mutate(N_variants_binned = case_when(
    (N_variants <= 1) ~ "(0, 1]",
    (N_variants <= 2) ~ "(1, 2]",
    (N_variants <= 5) ~ "(2, 5]",
    (N_variants <= 10) ~ "(5, 10]",
    (N_variants <= 1e6) ~ "(10, 1e6)",
    TRUE ~ "ERROR"
  )) %>% 
  rename(gene_id = gene_ids) %>% 
  left_join(lof_metrics %>% select(gene_id, oe_lof_upper))

maps_table_binned_2$N_variants_binned <- 
  factor(maps_table_binned_2$N_variants_binned,
         levels = c("(0, 1]", "(1, 2]", "(2, 5]", "(5, 10]", "(10, 1e6)"))

ggplot(maps_table_binned_2,
       aes(x = oe_lof_upper, y = maps, col = N_variants_binned)) +
  geom_point() +
  ggthemes::theme_hc() +
  ggtitle('Per Gene oe_lof_upper vs. MAPS based on HC variants') +
  ylab('MAPS') +
  xlab('OE_LOF_UPPER') +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18)
  )


# VOII) ecdf on n_variants in hc
maps_table_hc <- maps_table %>%
  filter(lof_csq_collapsed == 'HC')
# make ecdf of abs(maps)
ecdf_N_Variants <- ecdf(maps_table_hc$N_variants)
# annotate the ecdf to table using the ecdf of abs(maps)
maps_table_hc_ecdf <- maps_table_hc %>% 
  mutate(ecdf_N_Variants = ecdf_N_Variants(N_variants))
# plot
maps_table_hc_ecdf_zoom <- maps_table_hc_ecdf %>% 
  filter(N_variants < 20)
ggplot(maps_table_hc_ecdf_zoom, 
       aes(N_variants, 1 - ecdf_N_Variants)) +
  geom_step() +
  ggthemes::theme_hc() +
  xlab('N_variants') +
  ylab('Gene Ratio') +
  ggtitle('Cumulative Plot of N_Variants < 20 Per Gene Based on HC Variants')

# VII) Meta-gene
maps_table_meta <- maps_table %>% 
  rename(gene_id = gene_ids) %>% 
  left_join(lof_metrics) %>% 
  select(gene_id, lof_csq_collapsed, oe_lof_upper_bin,
         N_singletons, expected_singletons, N_variants) %>% 
  group_by(gene_id, lof_csq_collapsed, oe_lof_upper_bin) %>% 
  summarise(N_singletons_bin = sum(N_singletons),
            expected_singletons_bin = sum(expected_singletons),
            N_variants_bin = sum(N_variants)) %>% 
  mutate(MAPS_bin = (N_singletons_bin - expected_singletons_bin)/N_variants_bin)

# Calculate 95% confidence interval for the mean

# a)
MAPS_bin_mean <- mean(maps_table_meta$MAPS_bin)
MAPS_bin_sd <- sd(maps_table_meta$MAPS_bin)
MAPS_bin_n_obs <- length(maps_table_meta$MAPS_bin)
MAPS_bin_standard_error <- MAPS_bin_sd/sqrt(MAPS_bin_n_obs)
t_a <- qt(p = 0.05/2, df = MAPS_bin_n_obs - 1, lower.tail = F)
error_margin <- t_a*MAPS_bin_standard_error
# Calculate the lower bound 
lower_bound <- MAPS_bin_mean - error_margin
# Calculate the upper bound
upper_bound <- MAPS_bin_mean + error_margin
print(c(lower_bound, upper_bound))

# b)
# Calculate the mean and standard error
model_MAPS_bin <- lm(MAPS_bin ~ 1, maps_table_meta)
# Find the confidence interval
confint(model_MAPS_bin, level=0.95)

# Plot MAPS vs. oe_lof_decile
maps_table_meta <- maps_table_meta %>% 
  mutate(N_variants_binned = case_when(
    (N_variants_bin <= 1) ~ "(0, 1]",
    (N_variants_bin <= 2) ~ "(1, 2]",
    (N_variants_bin <= 5) ~ "(2, 5]",
    (N_variants_bin <= 10) ~ "(5, 10]",
    (N_variants_bin <= 1e6) ~ "(10, 1e6)",
    TRUE ~ "ERROR"
  ))


maps_table_meta$N_variants_binned <- maps_table_meta$N_variants_binned <- 
  factor(maps_table_meta$N_variants_binned,
         levels = c("(0, 1]", "(1, 2]", "(2, 5]", "(5, 10]", "(10, 1e6)"))
maps_table_meta$oe_lof_upper_bin <- factor(maps_table_meta$oe_lof_upper_bin,
                                           levels = unique(sort(maps_table_meta$oe_lof_upper_bin)))

ggplot(maps_table_meta,
       aes(x = oe_lof_upper_bin, y = MAPS_bin, col = oe_lof_upper_bin)) +
  geom_boxplot(width = 0.7) +
  facet_wrap(~ lof_csq_collapsed) + 
  ggthemes::theme_hc() +
  ggtitle('MAPS vs. oe_lof_upper_bin -
- Decile bin of LOEUF for given transcript (lower values indicate more constrained)') +
  ylab('MAPS') +
  xlab('oe_lof_upper_bin') +
  #scale_y_continuous(limits = c(0,10),
  #                   breaks = c(1:10)) +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18))


# VIII) Exploratory
maps_table_meta_hc_lc <- maps_table_meta %>% 
  filter(lof_csq_collapsed == 'HC' | lof_csq_collapsed == 'LC') %>%
  select(gene_id, lof_csq_collapsed, MAPS_bin) %>% 
  pivot_wider(
    names_from = lof_csq_collapsed,
    values_from = MAPS_bin
  )


maps_table_meta_hc_lc %>%
    ggplot(aes(x = HC, y = LC)) +
    geom_point(size = 1) +
    theme_classic() +
    ggtitle("MAPS in HC vs. MAPS in LC For The Same Genes") +
    xlab("MAPS_HC") +
    ylab("MAPS_LC") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(
      axis.title = element_text(face = "bold", size = 18),
      axis.text.x = element_text(face = "bold", colour = "black", size = 15),
      axis.text.y = element_text(face = "bold", colour = "black", size = 15),
      title = element_text(face = "bold", size = 20),
      legend.text = element_text(size = 18)
    )

maps_table_meta_hc_missense <- maps_table_meta %>% 
  filter(lof_csq_collapsed == 'HC' | lof_csq_collapsed == 'missense_variant') %>%
  select(gene_id, lof_csq_collapsed, MAPS_bin) %>% 
  pivot_wider(
    names_from = lof_csq_collapsed,
    values_from = MAPS_bin
  )


maps_table_meta_hc_missense %>%
  ggplot(aes(x = HC, y = missense_variant)) +
  geom_point(size = 1) +
  theme_classic() +
  ggtitle("MAPS in HC vs. MAPS in Missense SNVs For The Same Genes") +
  xlab("MAPS_HC") +
  ylab("MAPS_MISSENSE") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(
    axis.title = element_text(face = "bold", size = 18),
    axis.text.x = element_text(face = "bold", colour = "black", size = 15),
    axis.text.y = element_text(face = "bold", colour = "black", size = 15),
    title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 18)
  )


