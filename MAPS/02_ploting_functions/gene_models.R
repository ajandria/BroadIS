
library(tidyverse)
library(ggplot2)

maps_global <- read_csv('/Users/adrian/BroadIS/01_maps/data/maps_downsampling_multi_model.csv')

maps_global_log_glob <- maps_global %>% 
  filter(population == 'global')
























