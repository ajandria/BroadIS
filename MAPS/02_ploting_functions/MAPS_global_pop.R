library(tidyverse)
library(ggplot2)


maps_global <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/maps_downsampling_multi_model.csv')
#maps_global$downsampling <- factor(maps_global$downsampling, levels = unique(maps_global$downsampling))

ggplot(maps_global %>% filter(population == 'global'),
       aes(downsampling, MAPS, color = lof_csq_collapsed)) +
  scale_color_manual(values = c("darkred", "darkred", "orange", "grey")) +
  geom_point() +
  facet_grid(~ lof_csq_collapsed) +
  theme_bw() +
  scale_y_continuous(limits = c(-0.1, 0.3),
                     breaks = seq(-0.1, 0.3, 0.05)) +
  scale_x_continuous(limits = c(0, 75000),
                     breaks = c(0, 20000, 40000, 60000, 80000))

maps_global_log_glob <- maps_global %>% 
  filter(population == 'global') %>% 
  mutate(downsampling = log10(downsampling))

ggplot(maps_global_log_glob,
       aes(downsampling, MAPS, color = lof_csq_collapsed)) +
  scale_color_manual(values = c("darkred", "red", "orange", "grey")) +
  geom_point() +
  #facet_grid(~ lof_csq_collapsed) +
  theme_bw() +
  scale_y_continuous(limits = c(-0.1, 0.3),
                     breaks = seq(-0.1, 0.3, 0.05)) +
  geom_line() +
  xlab("log10(downsampling)") +
  ggtitle('Downsampling plot of MAPS vs. global population')
  #scale_x_continuous(limits = c(0, 75000),
  #                   breaks = c(0, 20000, 40000, 60000, 80000))

# -----------------------------------------------------------------------------

maps_pop <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/maps_downsampling_multi_model.csv') %>% 
  filter(population != 'global')
maps_pop$population <- factor(maps_pop$population, levels = unique(maps_pop$population))

# I)

ggplot(maps_pop,
       aes(downsampling, MAPS, color = population)) +
  scale_color_manual(values = c(nfe = '#6AA5CD', 
                                fin = '#002F6C', 
                                mid = '#33CC33',  # Middle Eastern  
                                oth = '#ABB9B9',
                                ami = 'turquoise', # Amish
                                asj = 'coral',
                                afr = '#941494',
                                eas = '#108C44',
                                sas = '#FF9912', 
                                amr = '#ED1E24')) +
  geom_point() +
  facet_grid(~ lof_csq_collapsed) +
  theme_bw() +
  scale_y_continuous(limits = c(-0.1, 0.4),
                     breaks = seq(-0.1, 0.3, 0.05)) +
  scale_x_continuous(limits = c(0, 75000),
                     breaks = c(0, 20000, 40000, 60000, 80000))

# II)
ggplot(maps_pop %>% mutate(downsampling = log10(downsampling)),
       aes(downsampling, MAPS, group = population, color = population)) +
  scale_color_manual(values = c(nfe = '#6AA5CD', 
                                fin = '#002F6C', 
                                mid = '#33CC33',  # Middle Eastern  
                                oth = '#ABB9B9',
                                ami = 'turquoise', # Amish
                                asj = 'coral',
                                afr = '#941494',
                                eas = '#108C44',
                                sas = '#FF9912', 
                                amr = '#ED1E24')) +
  geom_point() +
  facet_grid(~ lof_csq_collapsed) +
  theme_bw() +
  scale_y_continuous(limits = c(-0.1, 0.4),
                     breaks = seq(-0.1, 0.3, 0.05)) +
  geom_line() +
  xlab("log10(downsampling)") +
  ggtitle('Downsampling plot of MAPS vs. per population')
  #scale_x_continuous(limits = c(0, 75000),
  #                   breaks = c(0, 20000, 40000, 60000, 80000))



# Debugging ---------------------------------------------------------------
# Haplo
haploinsufficient_1  <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/haploinsufficiency_severe_curated_2016.tsv',
                                 col_names = c('gene'))
haploinsufficient_2  <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/haploinsufficiency_moderate_curated_2016.tsv',
                                 col_names = c('gene'))
haploinsufficient_3  <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/haploinsufficiency_severe_curated_2016.tsv',
                                 col_names = c('gene'))

haploinsufficient_final <- union(union(haploinsufficient_1$gene,
                                       haploinsufficient_2$gene), haploinsufficient_3$gene) %>% 
  tibble() %>% 
  rename(gene = ".") %>% 
  mutate(haploinsufficient = 1)
# Metrics
metrics <- read_delim('/Users/adrian/BroadIS/MAPS/01_data/gnomad.v2.1.1.lof_metrics.by_gene.txt')
# Main
maps_global <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/maps_clean_global_Per_gene.csv') %>% 
  rename(gene_id = gene_ids) %>% 
  #filter(lof_csq_collapsed == 'HC') %>% 
  left_join(metrics) %>% 
  left_join(haploinsufficient_final) %>% 
  mutate(haploinsufficient = factor(ifelse(is.na(haploinsufficient), 0, haploinsufficient)))

# Density plots with semi-transparent fill
ggplot(maps_global, aes(MAPS, fill = haploinsufficient)) + 
  geom_density(alpha=.3) +
  theme_bw() + 
  ggtitle('Density of HC per gene MAPS')

ggplot(maps_global, aes(oe_lof_upper, fill = haploinsufficient)) + 
  geom_density(alpha=.3) +
  theme_bw() + 
  ggtitle('Density of HC per gene oe_lof_upper')

ggplot(maps_global, aes(x = MAPS, y = oe_lof_upper, col = haploinsufficient, alpha = N_variants)) +
  geom_point() +
  theme_bw() +
  ggtitle('Scatter plot of MAPS vs. LOEUF')

ggplot(maps_global, aes(MAPS, fill = haploinsufficient)) + 
  geom_histogram(alpha=.5, bins = 100) +
  theme_bw() + 
  ggtitle('Histogram of HC per gene MAPS')


maps_global %>% 
  filter(is.na(MAPS))


# Debugging ---------------------------------------------------------------
# AR
ar  <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/all_ar.tsv',
                                 col_names = c('gene')) %>% 
  mutate(AR = 1)
# Metrics
metrics <- read_delim('/Users/adrian/BroadIS/MAPS/01_data/gnomad.v2.1.1.lof_metrics.by_gene.txt')
# Main
maps_global <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/maps_clean_global_Per_gene.csv') %>% 
  rename(gene_id = gene_ids) %>% 
  filter(lof_csq_collapsed == 'HC') %>% 
  left_join(metrics) %>% 
  left_join(ar) %>% 
  mutate(AR = factor(ifelse(is.na(AR), 0, AR)))

# Density plots with semi-transparent fill
ar1 <- ggplot(maps_global, aes(MAPS, fill = AR)) + 
  geom_density(alpha=.3) +
  theme_bw() + 
  ggtitle('AR Gene List | Density of HC per gene MAPS')

ar2 <- ggplot(maps_global, aes(oe_lof_upper, fill = AR)) + 
  geom_density(alpha=.3) +
  theme_bw() + 
  ggtitle('AR Gene List | Density of HC per gene oe_lof_upper')

ggplot(maps_global %>% arrange(AR), aes(x = MAPS, y = oe_lof_upper, col = AR, alpha = N_variants)) +
  ggsci::scale_color_lancet() +
  geom_point() +
  theme_bw() +
  ggtitle('AR Gene List | Scatter plot of MAPS vs. LOEUF by N_variants')

ggplot(maps_global, aes(MAPS, fill = AR)) + 
  geom_histogram(alpha=.5, bins = 100) +
  theme_bw() + 
  ggtitle('Histogram of HC per gene MAPS')



# Debugging ---------------------------------------------------------------
# AD
ad  <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/all_ad.tsv',
                col_names = c('gene')) %>% 
  mutate(AD = 1)
# Metrics
metrics <- read_delim('/Users/adrian/BroadIS/MAPS/01_data/gnomad.v2.1.1.lof_metrics.by_gene.txt')
# Main
maps_global <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/maps_clean_global_Per_gene.csv') %>% 
  rename(gene_id = gene_ids) %>% 
  filter(lof_csq_collapsed == 'HC') %>% 
  left_join(metrics) %>% 
  left_join(ad) %>% 
  mutate(AD = factor(ifelse(is.na(AD), 0, AD)))

# Density plots with semi-transparent fill
ad1 <- ggplot(maps_global, aes(MAPS, fill = AD)) + 
  geom_density(alpha=.3) +
  theme_bw() + 
  ggtitle('AD Gene List | Density of HC per gene MAPS')

ad2 <- ggplot(maps_global, aes(oe_lof_upper, fill = AD)) + 
  geom_density(alpha=.3) +
  theme_bw() + 
  ggtitle('AD Gene List | Density of HC per gene oe_lof_upper')

ggplot(maps_global %>% arrange(AD), aes(x = MAPS, y = oe_lof_upper, col = AD, alpha = N_variants)) +
  ggsci::scale_color_jco() +
  geom_point() +
  theme_bw() +
  ggtitle('AD Gene List | Scatter plot of MAPS vs. LOEUF by N_variants')

ggplot(maps_global, aes(MAPS, fill = AD)) + 
  geom_histogram(alpha=.5, bins = 100) +
  theme_bw() + 
  ggtitle('Histogram of HC per gene MAPS')










  