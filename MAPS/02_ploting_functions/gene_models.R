
# Dependencies ------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(pROC)
library(caret)
library(ggplot2)

# Setup -------------------------------------------------------------------
maps_global <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/maps_clean_global_Per_gene.csv') %>% 
  rename(gene_id = gene_ids)

metrics <- read_delim('/Users/adrian/BroadIS/MAPS/01_data/gnomad.v2.1.1.lof_metrics.by_gene.txt')

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
  
olfactory <- read_csv('/Users/adrian/BroadIS/MAPS/01_data/olfactory_receptors.tsv',
                                 col_names = c('gene'))

# Load --------------------------------------------------------------------
# Join the tables
maps_haplo <- maps_global %>% 
  left_join(metrics) %>% 
  left_join(haploinsufficient_final) %>% 
  mutate(haploinsufficient = factor(ifelse(is.na(haploinsufficient), 0, haploinsufficient)))

# Logit for just MAPS
# Logit for just oe_lof_upper
haploinsufficient_MAPS_logit <- 
  glm(haploinsufficient ~ MAPS, 
      data = maps_haplo %>% filter(lof_csq_collapsed == 'HC'), 
      family = "binomial")

knitr::kable(broom::tidy(haploinsufficient_MAPS_logit))

# Logit for just oe_lof_upper
haploinsufficient_lof_logit <- 
  glm(haploinsufficient ~ oe_lof_upper, 
      data = maps_haplo %>% filter(lof_csq_collapsed == 'HC'), 
      family = "binomial")

knitr::kable(broom::tidy(haploinsufficient_lof_logit))

# Add effect of MAPS
haploinsufficient_lof_maps_logit <- 
  glm(haploinsufficient ~ oe_lof_upper + MAPS, 
      data = maps_haplo %>% filter(lof_csq_collapsed == 'HC'), 
      family = "binomial")

knitr::kable(broom::tidy(haploinsufficient_lof_maps_logit))









