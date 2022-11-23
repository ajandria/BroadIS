library(tidyverse)
library(ggplot2)


maps_global <- read_csv('/Users/adrian/BroadIS/01_maps/data/MAPS_global_downsampling.csv')
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


maps_pop <- read_csv('/Users/adrian/BroadIS/01_maps/data/MAPS_population_downsampling.csv') %>% 
  filter(population != 'global')
maps_pop$population <- factor(maps_pop$population, levels = unique(maps_pop$population))

ggplot(maps_pop,
       aes(downsampling, MAPS, color = population)) +
  geom_point() +
  facet_grid(~ lof_csq_collapsed) +
  theme_bw() +
  scale_y_continuous(limits = c(-0.1, 0.3),
                     breaks = seq(-0.1, 0.3, 0.05)) +
  scale_x_continuous(limits = c(0, 75000),
                     breaks = c(0, 20000, 40000, 60000, 80000))
