library(tidyverse)
library(ggplot2)


maps_global <- read_csv('/Users/adrian/BroadIS/01_maps/data/maps_downsampling_multi_model.csv')
#maps_global$downsampling <- factor(maps_global$downsampling, levels = unique(maps_global$downsampling))

ggplot(maps_global,
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


# -------------------------------------------------------------------------

maps_pop <- read_csv('/Users/adrian/BroadIS/01_maps/data/maps_downsampling_multi_model.csv') %>% 
  filter(population != 'global')
maps_pop$population <- factor(maps_pop$population, levels = unique(maps_pop$population))

# I)

ggplot(maps_pop,
       aes(downsampling, MAPS, color = population)) +
  geom_point() +
  facet_grid(~ lof_csq_collapsed) +
  theme_bw() +
  scale_y_continuous(limits = c(-0.1, 0.3),
                     breaks = seq(-0.1, 0.3, 0.05)) +
  scale_x_continuous(limits = c(0, 75000),
                     breaks = c(0, 20000, 40000, 60000, 80000))

# II)

ggplot(maps_pop %>% mutate(downsampling = log10(downsampling)),
       aes(downsampling, MAPS, color = population)) +
  scale_color_manual(values = c('black', 'brown', 'green', 'blue', 'gold', 'red', 'magenta', 'purple',
                                'grey', 'orange')) +
  geom_point() +
  facet_grid(~ lof_csq_collapsed) +
  theme_bw() +
  scale_y_continuous(limits = c(-0.1, 0.3),
                     breaks = seq(-0.1, 0.3, 0.05)) +
  geom_line() +
  xlab("log10(downsampling)") +
  ggtitle('Downsampling plot of MAPS vs. per population')
  #scale_x_continuous(limits = c(0, 75000),
  #                   breaks = c(0, 20000, 40000, 60000, 80000))





