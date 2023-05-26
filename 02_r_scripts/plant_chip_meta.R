# a meta analysis of Arabidopsis clock ChIP data and overlap with cold TF network

library(tidyverse)
library(janitor)
library(MetaCycle)
devtools::install_github("vqv/ggbiplot")
library(ggbiplot)
library(ggrepel)
library(circlize)
library(UpSetR)

library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)

setwd("/Users/Allan/Documents/plant_ChIP_meta/")

# gather and join together all the WebPlotDigitizer files for each cluster group
# WebPlotDigitizer https://apps.automeris.io/wpd/
clusters_aggregated <- list.files(path = './00_raw_data/cluster_image_analysis_aggregate', 
                                  pattern = '*.csv', full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  purrr::reduce(full_join, by = 'id') %>% 
  remove_empty(which = 'cols') %>% 
  relocate(c(cluster_2, cluster_3, cluster_4, cluster_5, cluster_6, cluster_7, 
             cluster_8, cluster_9), .before = cluster_10) %>% 
  write_csv('./01_tidy_data/cluster_profiles_aggregated.csv')

# split by individual days
clusters_aggregated_day1 <- clusters_aggregated[row.names(clusters_aggregated) %in% 1:9, ]
clusters_aggregated_day2 <- clusters_aggregated[row.names(clusters_aggregated) %in% 9:17, ]
clusters_aggregated_day5 <- clusters_aggregated[row.names(clusters_aggregated) %in% 18:26, ]

# transpose so that columns are sample numbers (s1 to s9) and rows are cluster groups, and rename columns
clusters_aggregated_day1_df <- data.frame(t(clusters_aggregated_day1[, -1])) %>% 
  rename_with(~ paste0("s", 1:9)) %>% 
  mutate(clusters = paste0("cluster_", 1:74)) %>% 
  relocate(clusters) %>% 
  write_csv('./01_tidy_data/clusters_aggregated_day1.csv')

clusters_aggregated_day2_df <- data.frame(t(clusters_aggregated_day2[, -1])) %>% 
  rename_with(~ paste0("s", 9:17)) %>% 
  mutate(clusters = paste0("cluster_", 1:74)) %>% 
  relocate(clusters) %>% 
  write_csv('./01_tidy_data/clusters_aggregated_day2.csv')

clusters_aggregated_day5_df <- data.frame(t(clusters_aggregated_day5[, -1])) %>% 
  rename_with(~ paste0("s", 18:26)) %>% 
  mutate(clusters = paste0("cluster_", 1:74)) %>% 
  relocate(clusters) %>% 
  write_csv('./01_tidy_data/clusters_aggregated_day5.csv')

# use meta2d from MetaCycle package to detect rhythmic signals from time-series datasets with multiple methods
# https://cran.r-project.org/web/packages/MetaCycle/MetaCycle.pdf
# https://cran.r-project.org/web/packages/MetaCycle/vignettes/implementation.html
# output files are in specified outdir
meta2d(infile = './01_tidy_data/clusters_aggregated_day1.csv', filestyle = 'csv', outdir = './00_raw_data/cluster_days_output_day1', 
       timepoints = seq(0, 24, by = 3))

meta2d(infile = './01_tidy_data/clusters_aggregated_day2.csv', filestyle = 'csv', outdir = './00_raw_data/cluster_days_output_day2', 
       timepoints = seq(0, 24, by = 3))

meta2d(infile = './01_tidy_data/clusters_aggregated_day5.csv', filestyle = 'csv', outdir = './00_raw_data/cluster_days_output_day5', 
       timepoints = seq(0, 24, by = 3))

# read back in the meta2d output file with selected columns
day1_output <- read_csv('./00_raw_data/cluster_days_output_day1/meta2d_clusters_aggregated_day1.csv') %>% 
  select(CycID, meta2d_pvalue, meta2d_period, meta2d_phase, meta2d_AMP) %>% 
  dplyr::rename(cluster = CycID, d1_meta2d_pvalue = meta2d_pvalue, d1_meta2d_period = meta2d_period, d1_meta2d_phase = meta2d_phase, d1_meta2d_AMP = meta2d_AMP)

day2_output <- read_csv('./00_raw_data/cluster_days_output_day2/meta2d_clusters_aggregated_day2.csv') %>% 
  select(CycID, meta2d_pvalue, meta2d_period, meta2d_phase, meta2d_AMP) %>% 
  dplyr::rename(cluster = CycID, d2_meta2d_pvalue = meta2d_pvalue, d2_meta2d_period = meta2d_period, d2_meta2d_phase = meta2d_phase, d2_meta2d_AMP = meta2d_AMP)

day5_output <- read_csv('./00_raw_data/cluster_days_output_day5/meta2d_clusters_aggregated_day5.csv') %>% 
  select(CycID, meta2d_pvalue, meta2d_period, meta2d_phase, meta2d_AMP) %>% 
  dplyr::rename(cluster = CycID, d5_meta2d_pvalue = meta2d_pvalue, d5_meta2d_period = meta2d_period, d5_meta2d_phase = meta2d_phase, d5_meta2d_AMP = meta2d_AMP)

# compare d1 vs d2 and set rules for classifying amplitude and rhythm changes between days
# _150 means 1.5 fold up or down difference in amplitude
day1_vs_day2_150 <- full_join(day1_output, day2_output, by = "cluster") %>% 
  mutate(AMP_change = d2_meta2d_AMP/d1_meta2d_AMP *100) %>% 
  mutate(pval_flag = case_when(d1_meta2d_pvalue > 0.05 & d2_meta2d_pvalue > 0.05 ~ 'd1_nr_d2_nr', 
                               d1_meta2d_pvalue <= 0.05 & d2_meta2d_pvalue > 0.05 ~ 'd1_r_d2_nr',
                               d1_meta2d_pvalue > 0.05 & d2_meta2d_pvalue <= 0.05 ~ 'd1_nr_d2_r',
                               d1_meta2d_pvalue <= 0.05 & d2_meta2d_pvalue <= 0.05 ~ 'd1_r_d2_r',
                               TRUE ~ 'rhythmic')) %>% 
  mutate(amp_flag = case_when(AMP_change <= 66.6  & AMP_change > 33.3 ~ 'lose_medium',
                              AMP_change <= 33.3 ~ 'lose_high',
                              AMP_change >= 150 & AMP_change < 300 ~ 'gain_medium',
                              AMP_change >= 300 ~ 'gain_high',
                              TRUE ~ 'other')) %>% 
  mutate(amp_rhythm_flag = case_when(amp_flag == 'gain_high' & pval_flag == 'd1_r_d2_r' ~ 'gain_amp_high_rhy_stay',
                                     amp_flag == 'gain_medium' & pval_flag == 'd1_r_d2_r' ~ 'gain_amp_med_rhy_stay',
                                     amp_flag == 'gain_high' & pval_flag == 'd1_nr_d2_r' ~ 'gain_amp_high_rhy_gain',
                                     amp_flag == 'gain_medium' & pval_flag == 'd1_nr_d2_r' ~ 'gain_amp_med_rhy_gain',
                                     amp_flag == 'gain_high' & pval_flag == 'd1_r_d2_nr' ~ 'gain_amp_high_rhy_lose',
                                     amp_flag == 'gain_medium' & pval_flag == 'd1_r_d2_nr' ~ 'gain_amp_med_rhy_lose',
                                     amp_flag == 'gain_high' & pval_flag == 'd1_nr_d2_nr' ~ 'gain_amp_high_rhy_none',
                                     amp_flag == 'gain_medium' & pval_flag == 'd1_nr_d2_nr' ~ 'gain_amp_med_rhy_none',
                                     amp_flag == 'lose_high' & pval_flag == 'd1_r_d2_r' ~ 'lose_amp_high_rhy_stay',
                                     amp_flag == 'lose_medium' & pval_flag == 'd1_r_d2_r' ~ 'lose_amp_med_rhy_stay',
                                     amp_flag == 'lose_high' & pval_flag == 'd1_nr_d2_r' ~ 'lose_amp_high_rhy_gain',
                                     amp_flag == 'lose_medium' & pval_flag == 'd1_nr_d2_r' ~ 'lose_amp_med_rhy_gain',
                                     amp_flag == 'lose_high' & pval_flag == 'd1_r_d2_nr' ~ 'lose_amp_high_rhy_lose',
                                     amp_flag == 'lose_medium' & pval_flag == 'd1_r_d2_nr' ~ 'lose_amp_med_rhy_lose',
                                     amp_flag == 'lose_high' & pval_flag == 'd1_nr_d2_nr' ~ 'lose_amp_high_rhy_none',
                                     amp_flag == 'lose_medium' & pval_flag == 'd1_nr_d2_nr' ~ 'lose_amp_med_rhy_none',
                                     TRUE ~ 'other')) %>% 
  mutate(cluster_id = as.character(paste(1:74))) %>% 
  relocate(cluster_id)

write_csv(day1_vs_day2_150, './01_tidy_data/day1_vs_day2_150.csv')

# compare d1 vs d5 and set rules for classifying amplitude and rhythm changes between days
# _150 means 1.5 fold up or down difference in amplitude
day1_vs_day5_150 <- full_join(day1_output, day5_output, by = "cluster") %>% 
  mutate(AMP_change = d5_meta2d_AMP/d1_meta2d_AMP *100) %>% 
  mutate(pval_flag = case_when(d1_meta2d_pvalue > 0.05 & d5_meta2d_pvalue > 0.05 ~ 'd1_nr_d5_nr', 
                               d1_meta2d_pvalue <= 0.05 & d5_meta2d_pvalue > 0.05 ~ 'd1_r_d5_nr',
                               d1_meta2d_pvalue > 0.05 & d5_meta2d_pvalue <= 0.05 ~ 'd1_nr_d5_r',
                               d1_meta2d_pvalue <= 0.05 & d5_meta2d_pvalue <= 0.05 ~ 'd1_r_d5_r',
                               TRUE ~ 'rhythmic')) %>% 
  mutate(amp_flag = case_when(AMP_change <= 66.6  & AMP_change > 33.3 ~ 'lose_medium',
                              AMP_change <= 33.3 ~ 'lose_high',
                              AMP_change >= 150 & AMP_change < 300 ~ 'gain_medium',
                              AMP_change >= 300 ~ 'gain_high',
                              TRUE ~ 'other')) %>% 
  mutate(amp_rhythm_flag = case_when(amp_flag == 'gain_high' & pval_flag == 'd1_r_d5_r' ~ 'gain_amp_high_rhy_stay',
                                     amp_flag == 'gain_medium' & pval_flag == 'd1_r_d5_r' ~ 'gain_amp_med_rhy_stay',
                                     amp_flag == 'gain_high' & pval_flag == 'd1_nr_d5_r' ~ 'gain_amp_high_rhy_gain',
                                     amp_flag == 'gain_medium' & pval_flag == 'd1_nr_d5_r' ~ 'gain_amp_med_rhy_gain',
                                     amp_flag == 'gain_high' & pval_flag == 'd1_r_d5_nr' ~ 'gain_amp_high_rhy_lose',
                                     amp_flag == 'gain_medium' & pval_flag == 'd1_r_d5_nr' ~ 'gain_amp_med_rhy_lose',
                                     amp_flag == 'gain_high' & pval_flag == 'd1_nr_d5_nr' ~ 'gain_amp_high_rhy_none',
                                     amp_flag == 'gain_medium' & pval_flag == 'd1_nr_d5_nr' ~ 'gain_amp_med_rhy_none',
                                     amp_flag == 'lose_high' & pval_flag == 'd1_r_d5_r' ~ 'lose_amp_high_rhy_stay',
                                     amp_flag == 'lose_medium' & pval_flag == 'd1_r_d5_r' ~ 'lose_amp_med_rhy_stay',
                                     amp_flag == 'lose_high' & pval_flag == 'd1_nr_d5_r' ~ 'lose_amp_high_rhy_gain',
                                     amp_flag == 'lose_medium' & pval_flag == 'd1_nr_d5_r' ~ 'lose_amp_med_rhy_gain',
                                     amp_flag == 'lose_high' & pval_flag == 'd1_r_d5_nr' ~ 'lose_amp_high_rhy_lose',
                                     amp_flag == 'lose_medium' & pval_flag == 'd1_r_d5_nr' ~ 'lose_amp_med_rhy_lose',
                                     amp_flag == 'lose_high' & pval_flag == 'd1_nr_d5_nr' ~ 'lose_amp_high_rhy_none',
                                     amp_flag == 'lose_medium' & pval_flag == 'd1_nr_d5_nr' ~ 'lose_amp_med_rhy_none',
                                     TRUE ~ 'other')) %>% 
  mutate(cluster_id = as.character(paste(1:74))) %>% 
  relocate(cluster_id)

write_csv(day1_vs_day5_150, './01_tidy_data/day1_vs_day5_150.csv')

# Scatter plot of the MetaCycle outputs for Amplitude coloured by the amp_flag grouping
plot_day1_vs_day2_150_AMP <- day1_vs_day2_150 %>% 
  filter(pval_flag == 'd1_nr_d2_r' | pval_flag == 'd1_r_d2_nr' | pval_flag == 'd1_r_d2_r') %>%
  ggplot(aes(x = d1_meta2d_AMP, y = d2_meta2d_AMP, colour = amp_flag, shape = pval_flag)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  scale_x_continuous(limits = c(0, 1.8), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  geom_point(size = 2.5) +
  ggpubr::theme_pubr() +
  theme(legend.position = "right", legend.key = element_blank()) +
  scale_colour_brewer(palette = "Set1", labels = c("gain - high", "gain - medium", "lose - high", "lose - medium", "other")) +
  geom_text_repel(aes(label = cluster_id), show.legend = FALSE) +
  labs(colour = "Amplitude", shape = "Rhythm",
       y= "day 2 Amplitude", x = "day 1 Amplitude") +
  ggtitle("Day 1 (20C steady state) vs Day 2 (to 4C transient-cooling)",
          subtitle = "meta2d Amplitude") +
  theme(plot.title = element_text(color = "grey30"),
        plot.subtitle = element_text(color = "grey30")
  )

ggsave('./03_plots/plot_day1_vs_day2_150_AMP.png', dpi = 300, height = 6, width = 6, units = 'in')

plot_day1_vs_day5_150_AMP <- day1_vs_day5_150 %>%
  filter(pval_flag == 'd1_nr_d5_r' | pval_flag == 'd1_r_d5_nr' | pval_flag == 'd1_r_d5_r') %>%
  ggplot(aes(x = d1_meta2d_AMP, y = d5_meta2d_AMP, colour = amp_flag, shape = pval_flag)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  scale_x_continuous(limits = c(0, 1.8), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  geom_point(size = 2.5) +
  ggpubr::theme_pubr() +
  theme(legend.position = "right", legend.key = element_blank()) +
  scale_colour_brewer(palette = "Set1", labels = c("gain - high", "gain - medium", "lose - high", "lose - medium", "other")) +
  geom_text_repel(aes(label = cluster_id), show.legend = FALSE) +
  labs(colour = "Amplitude", shape = "Rhythm",
       y= "day 5 Amplitude", x = "day 1 Amplitude") +
  ggtitle("Day 1 (20C steady state) vs Day 5 (4C steady state)",
          subtitle = "meta2d Amplitude") +
  theme(plot.title = element_text(color = "grey30"),
        plot.subtitle = element_text(color = "grey30")
  )

ggsave('./03_plots/plot_day1_vs_day5_150_AMP.png', dpi = 300, height = 6, width = 6, units = 'in')


get_clusters <- function(df, filter_col, amp_flag_id){
  
  clusters <- df %>% 
    dplyr::filter({{filter_col}} == {{amp_flag_id}}) %>% 
    dplyr::select(cluster_id) %>% 
    dplyr::rename(cluster = cluster_id)
  
  return(clusters)
}

# day1 vs day 2
amp_gain_high_clusters_d1_d2 <- get_clusters(day1_vs_day2_150, amp_flag, 'gain_high')
amp_gain_medium_clusters_d1_d2 <- get_clusters(day1_vs_day2_150, amp_flag, 'gain_medium')
amp_lose_high_clusters_d1_d2 <- get_clusters(day1_vs_day2_150, amp_flag, 'lose_high')
amp_lose_medium_clusters_d1_d2 <- get_clusters(day1_vs_day2_150, amp_flag, 'lose_medium')
amp_other_clusters_d1_d2 <- get_clusters(day1_vs_day2_150, amp_flag, 'other')

# day1 vs day 5
amp_gain_high_clusters_d1_d5 <- get_clusters(day1_vs_day5_150, amp_flag, 'gain_high')
amp_gain_medium_clusters_d1_d5 <- get_clusters(day1_vs_day5_150, amp_flag, 'gain_medium')
amp_lose_high_clusters_d1_d5 <- get_clusters(day1_vs_day5_150, amp_flag, 'lose_high')
amp_lose_medium_clusters_d1_d5 <- get_clusters(day1_vs_day5_150, amp_flag, 'lose_medium')
amp_other_clusters_d1_d5 <- get_clusters(day1_vs_day5_150, amp_flag, 'other')

# An analysis of established and published Arabidopsis clock ChIP targets in TF network
# read the TF network cluster
# The TF network is 7302 genes over 75 clusters (clusters 0-74)

TF <- read_csv("./00_raw_data/TF Network Cluster Nov2018.csv") %>%
  pivot_longer(cols = starts_with ('cluster'),
               names_to = "cluster", 
               values_to = "gene_ID",
               values_drop_na = TRUE) %>% 
  filter(grepl('AT', gene_ID)) %>%
  mutate(cluster = str_sub(cluster, 9, -1)) %>% 
  write_csv("./01_tidy_data/TF Network Cluster Nov2018 pivot longer.csv")

# LHY dataset
# read in the Adams LHY paper dataset and skip first 2 lines
# The LHY paper is Adams et al. (2018) New Phytologist 220(3); 897
# supplemental data set (Table S2) 
adams <- read_csv("./00_raw_data/nph15415-sup-0002-tables2.csv",skip=2) %>% 
  dplyr::select(gene_ID = 1) %>%
  filter(!is.na(gene_ID)) %>%
  distinct(gene_ID) %>% 
  write_csv("./01_tidy_data/LHY_targets.csv")

# CCA1 Nagel et al. dataset
# read in the Nagel CCA1 paper dataset and skip first 2 lines
# The CCA1 paper is Nagel et al. (2015) PNAS 112(34); E4802
# supplemental data set (Table S1) 
nagel <- read_csv("./00_raw_data/pnas.1513609112.sd01.csv",skip=2) %>% 
  dplyr::select(gene_ID = 10) %>% 
  mutate(gene_ID = str_sub(gene_ID, end = 9)) %>% 
  distinct(gene_ID) %>% 
  write_csv("./01_tidy_data/CCA1_nagel_targets.csv")

# CCA1 Kamioka et al. dataset
# read in the Kamioka CCA1 paper dataset and skip first 2 lines
# The CCA1 paper is Kamioka et al. (2016) Plant Cell 28(3); 696
# supplemental data set (Table S1C)
kamioka <- read_csv("./00_raw_data/TPC2015-00737-RAR3_Supplemental_Data_Set_1C.csv",skip=3) %>% 
  dplyr::select(gene_ID = 10) %>% 
  distinct(gene_ID) %>% 
  write_csv("./01_tidy_data/CCA1_kamioka_targets.csv")

# merge the nagel and kamioka CCA1 datasets
# use inner_join from dplyr
# 249 obs.
kamioka_nagel_merge <- inner_join(nagel, kamioka, by = "gene_ID")

# TOC1 dataset
# read in the Huang TOC1 paper dataset
# The TOC1 paper is Huang et al. (2012) Science 336:75
# supplemental data set (Table S1)
huang <- read_csv("./00_raw_data/Huang TOC1 CHiP TableS1.csv") %>% 
  dplyr::select(gene_ID = 14) %>% 
  distinct(gene_ID) %>% 
  write_csv("./01_tidy_data/TOC1_huang_targets.csv")

# PRR5 dataset
# read in the Nakamichi PRR5 paper dataset
# The PRR5 paper is Nakamichi et al. (2012) PNAS 109:17123
# supplemental data set (Table S3) 
nakamichi <- read_csv("./00_raw_data/Dataset S3 Nakamichi et al PRR5 binding targets PNAS 2012.csv", skip=2) %>% 
  dplyr::select(gene_ID = 3) %>% 
  distinct(gene_ID) %>%
  write_csv("./01_tidy_data/PRR5_nakamichi_targets.csv")

# PRR7 dataset
# read in the Liu PRR7 paper dataset
# The PRR7 paper is Liu et al. (2013) The Plant Journal 76:101
# supplemental data set (Table S1)
liu <- read_csv("./00_raw_data/Dataset S1 Liu et al PRR7 edit.csv") %>% 
  dplyr::select(gene_ID = 17) %>% 
  distinct(gene_ID) %>%
  write_csv("./01_tidy_data/PRR7_liu_targets.csv")

# LUX dataset
# read in the Ezer EC paper for the LUX dataset (LUX_17 tab)
# The Evening Complex (EC) paper is Ezer et al. (2017) Nature Plants 3: article 17087
# supplemental data set (LUX_17 tab of Table S6) 
ezer_LUX <- read_csv("./00_raw_data/Ezer et al nplants Suppl Table S6.csv") %>% 
  dplyr::select(gene_ID = 1) %>% 
  distinct(gene_ID) %>% 
  write_csv("./01_tidy_data/LUX_ezer_targets.csv")

# ELF3 dataset
# read in the Ezer EC paper for the ELF3 dataset (ELF3_22 tab)
# The Evening Complex (EC) paper is Ezer et al. (2017) Nature Plants 3: article 17087
# supplemental data set (ELF3_22 tab of Table S6) 
ezer_ELF3 <- read_csv("./00_raw_data/ELF3_22 Ezer Table S6.csv") %>% 
  dplyr::select(gene_ID = 1) %>% 
  distinct(gene_ID) %>%
  write_csv("./01_tidy_data/ELF3_ezer_targets.csv")

# ELF4 dataset
# read in the Ezer EC paper for the ELF4 dataset (ELF4_22 tab)
# The Evening Complex (EC) paper is Ezer et al. (2017) Nature Plants 3: article 17087
# supplemental data set (ELF4_22 tab of Table S6) 
ezer_ELF4 <- read_csv("./00_raw_data/ELF4_22 Ezer Table S6.csv") %>% 
  dplyr::select(gene_ID = 1) %>% 
  distinct(gene_ID) %>%
  write_csv("./01_tidy_data/ELF4_ezer_targets.csv")


merge_TF_clock <- function(df, clock_id){
  
  merge <- inner_join(TF, df, by = 'gene_ID') %>%
    arrange(nchar(cluster), cluster) %>% 
    mutate(clock = {{clock_id}})
  
  return(merge)
}

TF_adams_merge <- merge_TF_clock(adams, 'LHY')
TF_nagel_merge <- merge_TF_clock(nagel, 'CCA1')
TF_kamioka_merge <- merge_TF_clock(kamioka, 'CCA1')
TF_kamioka_nagel_merge <- merge_TF_clock(kamioka_nagel_merge, 'CCA1')
TF_huang_merge <- merge_TF_clock(huang, 'TOC1')
TF_nakamichi_merge <- merge_TF_clock(nakamichi, 'PRR5')
TF_liu_merge <- merge_TF_clock(liu, 'PRR7')
TF_ezer_LUX_merge <- merge_TF_clock(ezer_LUX, 'LUX')
TF_ezer_ELF3_merge <- merge_TF_clock(ezer_ELF3, 'ELF3')
TF_ezer_ELF4_merge <- merge_TF_clock(ezer_ELF4, 'ELF4')


clock_clusters <- function(df_clock, df_clusters, label){
  
  merge_clock_cluster_types <- df_clock %>% 
    inner_join(df_clusters, by = 'cluster') %>% 
    mutate(type = {{label}})
  
  return(merge_clock_cluster_types)
  
}

# d1-d2----
# gain-high d1d2
LHY_gain_high_d1d2 <- clock_clusters(TF_adams_merge, amp_gain_high_clusters_d1_d2, 'gain_high_d1_d2')
CCA1_nagel_gain_high_d1d2 <- clock_clusters(TF_nagel_merge, amp_gain_high_clusters_d1_d2, 'gain_high_d1_d2')
CCA1_kamioka_gain_high_d1d2 <- clock_clusters(TF_kamioka_merge, amp_gain_high_clusters_d1_d2, 'gain_high_d1_d2')
CCA1_nagel_kamioka_gain_high_d1d2 <- clock_clusters(TF_kamioka_nagel_merge, amp_gain_high_clusters_d1_d2, 'gain_high_d1_d2')
TOC1_gain_high_d1d2 <- clock_clusters(TF_huang_merge, amp_gain_high_clusters_d1_d2, 'gain_high_d1_d2')
PRR5_gain_high_d1d2 <- clock_clusters(TF_nakamichi_merge, amp_gain_high_clusters_d1_d2, 'gain_high_d1_d2')
PRR7_gain_high_d1d2 <- clock_clusters(TF_liu_merge, amp_gain_high_clusters_d1_d2, 'gain_high_d1_d2')
LUX_gain_high_d1d2 <- clock_clusters(TF_ezer_LUX_merge, amp_gain_high_clusters_d1_d2, 'gain_high_d1_d2')
ELF3_gain_high_d1d2 <- clock_clusters(TF_ezer_ELF3_merge, amp_gain_high_clusters_d1_d2, 'gain_high_d1_d2')
ELF4_gain_high_d1d2 <- clock_clusters(TF_ezer_ELF4_merge, amp_gain_high_clusters_d1_d2, 'gain_high_d1_d2')

# gain-medium d1d2
LHY_gain_medium_d1d2 <- clock_clusters(TF_adams_merge, amp_gain_medium_clusters_d1_d2, 'gain_medium_d1_d2')
CCA1_nagel_gain_medium_d1d2 <- clock_clusters(TF_nagel_merge, amp_gain_medium_clusters_d1_d2, 'gain_medium_d1_d2')
CCA1_kamioka_gain_medium_d1d2 <- clock_clusters(TF_kamioka_merge, amp_gain_medium_clusters_d1_d2, 'gain_medium_d1_d2')
CCA1_nagel_kamioka_gain_medium_d1d2 <- clock_clusters(TF_kamioka_nagel_merge, amp_gain_medium_clusters_d1_d2, 'gain_medium_d1_d2')
TOC1_gain_medium_d1d2 <- clock_clusters(TF_huang_merge, amp_gain_medium_clusters_d1_d2, 'gain_medium_d1_d2')
PRR5_gain_medium_d1d2 <- clock_clusters(TF_nakamichi_merge, amp_gain_medium_clusters_d1_d2, 'gain_medium_d1_d2')
PRR7_gain_medium_d1d2 <- clock_clusters(TF_liu_merge, amp_gain_medium_clusters_d1_d2, 'gain_medium_d1_d2')
LUX_gain_medium_d1d2 <- clock_clusters(TF_ezer_LUX_merge, amp_gain_medium_clusters_d1_d2, 'gain_medium_d1_d2')
ELF3_gain_medium_d1d2 <- clock_clusters(TF_ezer_ELF3_merge, amp_gain_medium_clusters_d1_d2, 'gain_medium_d1_d2')
ELF4_gain_medium_d1d2 <- clock_clusters(TF_ezer_ELF4_merge, amp_gain_medium_clusters_d1_d2, 'gain_medium_d1_d2')

# lose-high d1d2
LHY_lose_high_d1d2 <- clock_clusters(TF_adams_merge, amp_lose_high_clusters_d1_d2, 'lose_high_d1_d2')
CCA1_nagel_lose_high_d1d2 <- clock_clusters(TF_nagel_merge, amp_lose_high_clusters_d1_d2, 'lose_high_d1_d2')
CCA1_kamioka_lose_high_d1d2 <- clock_clusters(TF_kamioka_merge, amp_lose_high_clusters_d1_d2, 'lose_high_d1_d2')
CCA1_nagel_kamioka_lose_high_d1d2 <- clock_clusters(TF_kamioka_nagel_merge, amp_lose_high_clusters_d1_d2, 'lose_high_d1_d2')
TOC1_lose_high_d1d2 <- clock_clusters(TF_huang_merge, amp_lose_high_clusters_d1_d2, 'lose_high_d1_d2')
PRR5_lose_high_d1d2 <- clock_clusters(TF_nakamichi_merge, amp_lose_high_clusters_d1_d2, 'lose_high_d1_d2')
PRR7_lose_high_d1d2 <- clock_clusters(TF_liu_merge, amp_lose_high_clusters_d1_d2, 'lose_high_d1_d2')
LUX_lose_high_d1d2 <- clock_clusters(TF_ezer_LUX_merge, amp_lose_high_clusters_d1_d2, 'lose_high_d1_d2')
ELF3_lose_high_d1d2 <- clock_clusters(TF_ezer_ELF3_merge, amp_lose_high_clusters_d1_d2, 'lose_high_d1_d2')
ELF4_lose_high_d1d2 <- clock_clusters(TF_ezer_ELF4_merge, amp_lose_high_clusters_d1_d2, 'lose_high_d1_d2')

# lose-medium d1d2
LHY_lose_medium_d1d2 <- clock_clusters(TF_adams_merge, amp_lose_medium_clusters_d1_d2, 'lose_medium_d1_d2')
CCA1_nagel_lose_medium_d1d2 <- clock_clusters(TF_nagel_merge, amp_lose_medium_clusters_d1_d2, 'lose_medium_d1_d2')
CCA1_kamioka_lose_medium_d1d2 <- clock_clusters(TF_kamioka_merge, amp_lose_medium_clusters_d1_d2, 'lose_medium_d1_d2')
CCA1_nagel_kamioka_lose_medium_d1d2 <- clock_clusters(TF_kamioka_nagel_merge, amp_lose_medium_clusters_d1_d2, 'lose_medium_d1_d2')
TOC1_lose_medium_d1d2 <- clock_clusters(TF_huang_merge, amp_lose_medium_clusters_d1_d2, 'lose_medium_d1_d2')
PRR5_lose_medium_d1d2 <- clock_clusters(TF_nakamichi_merge, amp_lose_medium_clusters_d1_d2, 'lose_medium_d1_d2')
PRR7_lose_medium_d1d2 <- clock_clusters(TF_liu_merge, amp_lose_medium_clusters_d1_d2, 'lose_medium_d1_d2')
LUX_lose_medium_d1d2 <- clock_clusters(TF_ezer_LUX_merge, amp_lose_medium_clusters_d1_d2, 'lose_medium_d1_d2')
ELF3_lose_medium_d1d2 <- clock_clusters(TF_ezer_ELF3_merge, amp_lose_medium_clusters_d1_d2, 'lose_medium_d1_d2')
ELF4_lose_medium_d1d2 <- clock_clusters(TF_ezer_ELF4_merge, amp_lose_medium_clusters_d1_d2, 'lose_medium_d1_d2')

# other d1d2
LHY_other_d1d2 <- clock_clusters(TF_adams_merge, amp_other_clusters_d1_d2, 'other_d1_d2')
CCA1_nagel_other_d1d2 <- clock_clusters(TF_nagel_merge, amp_other_clusters_d1_d2, 'other_d1_d2')
CCA1_kamioka_other_d1d2 <- clock_clusters(TF_kamioka_merge, amp_other_clusters_d1_d2, 'other_d1_d2')
CCA1_nagel_kamioka_other_d1d2 <- clock_clusters(TF_kamioka_nagel_merge, amp_other_clusters_d1_d2, 'other_d1_d2')
TOC1_other_d1d2 <- clock_clusters(TF_huang_merge, amp_other_clusters_d1_d2, 'other_d1_d2')
PRR5_other_d1d2 <- clock_clusters(TF_nakamichi_merge, amp_other_clusters_d1_d2, 'other_d1_d2')
PRR7_other_d1d2 <- clock_clusters(TF_liu_merge, amp_other_clusters_d1_d2, 'other_d1_d2')
LUX_other_d1d2 <- clock_clusters(TF_ezer_LUX_merge, amp_other_clusters_d1_d2, 'other_d1_d2')
ELF3_other_d1d2 <- clock_clusters(TF_ezer_ELF3_merge, amp_other_clusters_d1_d2, 'other_d1_d2')
ELF4_other_d1d2 <- clock_clusters(TF_ezer_ELF4_merge, amp_other_clusters_d1_d2, 'other_d1_d2')

# d1-d5----
# gain-high d1d5
LHY_gain_high_d1d5 <- clock_clusters(TF_adams_merge, amp_gain_high_clusters_d1_d5, 'gain_high_d1_d5')
CCA1_nagel_gain_high_d1d5 <- clock_clusters(TF_nagel_merge, amp_gain_high_clusters_d1_d5, 'gain_high_d1_d5')
CCA1_kamioka_gain_high_d1d5 <- clock_clusters(TF_kamioka_merge, amp_gain_high_clusters_d1_d5, 'gain_high_d1_d5')
CCA1_nagel_kamioka_gain_high_d1d5 <- clock_clusters(TF_kamioka_nagel_merge, amp_gain_high_clusters_d1_d5, 'gain_high_d1_d5')
TOC1_gain_high_d1d5 <- clock_clusters(TF_huang_merge, amp_gain_high_clusters_d1_d5, 'gain_high_d1_d5')
PRR5_gain_high_d1d5 <- clock_clusters(TF_nakamichi_merge, amp_gain_high_clusters_d1_d5, 'gain_high_d1_d5')
PRR7_gain_high_d1d5 <- clock_clusters(TF_liu_merge, amp_gain_high_clusters_d1_d5, 'gain_high_d1_d5')
LUX_gain_high_d1d5 <- clock_clusters(TF_ezer_LUX_merge, amp_gain_high_clusters_d1_d5, 'gain_high_d1_d5')
ELF3_gain_high_d1d5 <- clock_clusters(TF_ezer_ELF3_merge, amp_gain_high_clusters_d1_d5, 'gain_high_d1_d5')
ELF4_gain_high_d1d5 <- clock_clusters(TF_ezer_ELF4_merge, amp_gain_high_clusters_d1_d5, 'gain_high_d1_d5')

# gain-medium d1d5
LHY_gain_medium_d1d5 <- clock_clusters(TF_adams_merge, amp_gain_medium_clusters_d1_d5, 'gain_medium_d1_d5')
CCA1_nagel_gain_medium_d1d5 <- clock_clusters(TF_nagel_merge, amp_gain_medium_clusters_d1_d5, 'gain_medium_d1_d5')
CCA1_kamioka_gain_medium_d1d5 <- clock_clusters(TF_kamioka_merge, amp_gain_medium_clusters_d1_d5, 'gain_medium_d1_d5')
CCA1_nagel_kamioka_gain_medium_d1d5 <- clock_clusters(TF_kamioka_nagel_merge, amp_gain_medium_clusters_d1_d5, 'gain_medium_d1_d5')
TOC1_gain_medium_d1d5 <- clock_clusters(TF_huang_merge, amp_gain_medium_clusters_d1_d5, 'gain_medium_d1_d5')
PRR5_gain_medium_d1d5 <- clock_clusters(TF_nakamichi_merge, amp_gain_medium_clusters_d1_d5, 'gain_medium_d1_d5')
PRR7_gain_medium_d1d5 <- clock_clusters(TF_liu_merge, amp_gain_medium_clusters_d1_d5, 'gain_medium_d1_d5')
LUX_gain_medium_d1d5 <- clock_clusters(TF_ezer_LUX_merge, amp_gain_medium_clusters_d1_d5, 'gain_medium_d1_d5')
ELF3_gain_medium_d1d5 <- clock_clusters(TF_ezer_ELF3_merge, amp_gain_medium_clusters_d1_d5, 'gain_medium_d1_d5')
ELF4_gain_medium_d1d5 <- clock_clusters(TF_ezer_ELF4_merge, amp_gain_medium_clusters_d1_d5, 'gain_medium_d1_d5')

# lose-high d1d5
LHY_lose_high_d1d5 <- clock_clusters(TF_adams_merge, amp_lose_high_clusters_d1_d5, 'lose_high_d1_d5')
CCA1_nagel_lose_high_d1d5 <- clock_clusters(TF_nagel_merge, amp_lose_high_clusters_d1_d5, 'lose_high_d1_d5')
CCA1_kamioka_lose_high_d1d5 <- clock_clusters(TF_kamioka_merge, amp_lose_high_clusters_d1_d5, 'lose_high_d1_d5')
CCA1_nagel_kamioka_lose_high_d1d5 <- clock_clusters(TF_kamioka_nagel_merge, amp_lose_high_clusters_d1_d5, 'lose_high_d1_d5')
TOC1_lose_high_d1d5 <- clock_clusters(TF_huang_merge, amp_lose_high_clusters_d1_d5, 'lose_high_d1_d5')
PRR5_lose_high_d1d5 <- clock_clusters(TF_nakamichi_merge, amp_lose_high_clusters_d1_d5, 'lose_high_d1_d5')
PRR7_lose_high_d1d5 <- clock_clusters(TF_liu_merge, amp_lose_high_clusters_d1_d5, 'lose_high_d1_d5')
LUX_lose_high_d1d5 <- clock_clusters(TF_ezer_LUX_merge, amp_lose_high_clusters_d1_d5, 'lose_high_d1_d5')
ELF3_lose_high_d1d5 <- clock_clusters(TF_ezer_ELF3_merge, amp_lose_high_clusters_d1_d5, 'lose_high_d1_d5')
ELF4_lose_high_d1d5 <- clock_clusters(TF_ezer_ELF4_merge, amp_lose_high_clusters_d1_d5, 'lose_high_d1_d5')

# lose-medium d1d5
LHY_lose_medium_d1d5 <- clock_clusters(TF_adams_merge, amp_lose_medium_clusters_d1_d5, 'lose_medium_d1_d5')
CCA1_nagel_lose_medium_d1d5 <- clock_clusters(TF_nagel_merge, amp_lose_medium_clusters_d1_d5, 'lose_medium_d1_d5')
CCA1_kamioka_lose_medium_d1d5 <- clock_clusters(TF_kamioka_merge, amp_lose_medium_clusters_d1_d5, 'lose_medium_d1_d5')
CCA1_nagel_kamioka_lose_medium_d1d5 <- clock_clusters(TF_kamioka_nagel_merge, amp_lose_medium_clusters_d1_d5, 'lose_medium_d1_d5')
TOC1_lose_medium_d1d5 <- clock_clusters(TF_huang_merge, amp_lose_medium_clusters_d1_d5, 'lose_medium_d1_d5')
PRR5_lose_medium_d1d5 <- clock_clusters(TF_nakamichi_merge, amp_lose_medium_clusters_d1_d5, 'lose_medium_d1_d5')
PRR7_lose_medium_d1d5 <- clock_clusters(TF_liu_merge, amp_lose_medium_clusters_d1_d5, 'lose_medium_d1_d5')
LUX_lose_medium_d1d5 <- clock_clusters(TF_ezer_LUX_merge, amp_lose_medium_clusters_d1_d5, 'lose_medium_d1_d5')
ELF3_lose_medium_d1d5 <- clock_clusters(TF_ezer_ELF3_merge, amp_lose_medium_clusters_d1_d5, 'lose_medium_d1_d5')
ELF4_lose_medium_d1d5 <- clock_clusters(TF_ezer_ELF4_merge, amp_lose_medium_clusters_d1_d5, 'lose_medium_d1_d5')

# other d1d5
LHY_other_d1d5 <- clock_clusters(TF_adams_merge, amp_other_clusters_d1_d5, 'other_d1_d5')
CCA1_nagel_other_d1d5 <- clock_clusters(TF_nagel_merge, amp_other_clusters_d1_d5, 'other_d1_d5')
CCA1_kamioka_other_d1d5 <- clock_clusters(TF_kamioka_merge, amp_other_clusters_d1_d5, 'other_d1_d5')
CCA1_nagel_kamioka_other_d1d5 <- clock_clusters(TF_kamioka_nagel_merge, amp_other_clusters_d1_d5, 'other_d1_d5')
TOC1_other_d1d5 <- clock_clusters(TF_huang_merge, amp_other_clusters_d1_d5, 'other_d1_d5')
PRR5_other_d1d5 <- clock_clusters(TF_nakamichi_merge, amp_other_clusters_d1_d5, 'other_d1_d5')
PRR7_other_d1d5 <- clock_clusters(TF_liu_merge, amp_other_clusters_d1_d5, 'other_d1_d5')
LUX_other_d1d5 <- clock_clusters(TF_ezer_LUX_merge, amp_other_clusters_d1_d5, 'other_d1_d5')
ELF3_other_d1d5 <- clock_clusters(TF_ezer_ELF3_merge, amp_other_clusters_d1_d5, 'other_d1_d5')
ELF4_other_d1d5 <- clock_clusters(TF_ezer_ELF4_merge, amp_other_clusters_d1_d5, 'other_d1_d5')

# bind all clock targets----
# * d1d2----
# append all the LHY targets:
# LHY amp_gain_high (45 obs) + LHY amp_gain_medium (100 obs) + LHY amp_lose_high (52 obs) + LHY amp_lose_medium (20 obs) + LHY amp_other (110 obs): total equals 327 obs
LHY_bind_d1d2 <- bind_rows(LHY_gain_high_d1d2, LHY_gain_medium_d1d2, LHY_lose_high_d1d2, LHY_lose_medium_d1d2, LHY_other_d1d2)

names(LHY_bind_d1d2)[3] <- "LHY"

write_csv(LHY_bind_d1d2, './01_tidy_data/LHY_targets_d1d2.csv')

# CCA1-Nagel amp_gain_high (69 obs) + CCA1-Nagel amp_gain_medium (151 obs) + CCA1-Nagel amp_lose_high (97 obs) + CCA1-Nagel amp_lose_medium (33 obs) + CCA1-Nagel amp_other (183 obs): total equals 533 obs
CCA1_nagel_bind_d1d2 <- bind_rows(CCA1_nagel_gain_high_d1d2, CCA1_nagel_gain_medium_d1d2, CCA1_nagel_lose_high_d1d2, CCA1_nagel_lose_medium_d1d2, CCA1_nagel_other_d1d2)

names(CCA1_nagel_bind_d1d2)[3] <- "CCA1 Nagel"

write_csv(CCA1_nagel_bind_d1d2, './01_tidy_data/CCA1_nagel_bind_d1d2.csv')

# CCA1-Kamioka amp_gain_high (37 obs) + CCA1-Kamioka amp_gain_medium (59 obs) + CCA1-Kamioka amp_lose_high (24 obs) + CCA1-Kamioka amp_lose_medium (14 obs) + CCA1-Kamioka amp_other (68 obs): total equals 202 obs
CCA1_kamioka_bind_d1d2 <- bind_rows(CCA1_kamioka_gain_high_d1d2, CCA1_kamioka_gain_medium_d1d2, CCA1_kamioka_lose_high_d1d2, CCA1_kamioka_lose_medium_d1d2, CCA1_kamioka_other_d1d2)

names(CCA1_kamioka_bind_d1d2)[3] <- "CCA1 Kamioka"

write_csv(CCA1_kamioka_bind_d1d2, './01_tidy_data/CCA1_kamioka_bind_d1d2.csv')

# CCA1-Nagel-Kamioka amp_gain_high (22 obs) + CCA1-Nagel-Kamioka amp_gain_medium (44 obs) + CCA1-Nagel-Kamioka amp_lose_high (10 obs) + CCA1-Nagel-Kamioka amp_lose_medium (7 obs) + CCA1-Nagel-Kamioka amp_other (50 obs): total equals 133 obs
CCA1_nagel_kamioka_bind_d1d2 <- bind_rows(CCA1_nagel_kamioka_gain_high_d1d2, CCA1_nagel_kamioka_gain_medium_d1d2, CCA1_nagel_kamioka_lose_high_d1d2, CCA1_nagel_kamioka_lose_medium_d1d2, CCA1_nagel_kamioka_other_d1d2)

names(CCA1_nagel_kamioka_bind_d1d2)[3] <- "CCA1 Nagel-Kamioka"

write_csv(CCA1_nagel_kamioka_bind_d1d2, './01_tidy_data/CCA1_nagel_kamioka_bind_d1d2.csv')

# TOC1 amp_gain_high (79 obs) + TOC1 amp_gain_medium (73 obs) + TOC1 amp_lose_high (45 obs) + TOC1 amp_lose_medium (30 obs) + TOC1 amp_other (65 obs): total equals 292 obs
TOC1_bind_d1d2 <- bind_rows(TOC1_gain_high_d1d2, TOC1_gain_medium_d1d2, TOC1_lose_high_d1d2, TOC1_lose_medium_d1d2, TOC1_other_d1d2)

names(TOC1_bind_d1d2)[3] <- "TOC1"

write_csv(TOC1_bind_d1d2, './01_tidy_data/TOC1_bind_d1d2.csv')

# PRR5 amp_gain_high (30 obs) + PRR5 amp_gain_medium (14 obs) + PRR5 amp_lose_high (2 obs) + PRR5 amp_lose_medium (4 obs) + PRR5 amp_other (6 obs): total equals 56 obs
PRR5_bind_d1d2 <- bind_rows(PRR5_gain_high_d1d2, PRR5_gain_medium_d1d2, PRR5_lose_high_d1d2, PRR5_lose_medium_d1d2, PRR5_other_d1d2)

names(PRR5_bind_d1d2)[3] <- "PRR5"

write_csv(PRR5_bind_d1d2, './01_tidy_data/PRR5_bind_d1d2.csv')

# PRR7 amp_gain_high (29 obs) + PRR7 amp_gain_medium (10 obs) + PRR7 amp_lose_high (11 obs) + PRR7 amp_lose_medium (2 obs) + PRR7 amp_other (14 obs): total equals 66 obs
PRR7_bind_d1d2 <- bind_rows(PRR7_gain_high_d1d2, PRR7_gain_medium_d1d2, PRR7_lose_high_d1d2, PRR7_lose_medium_d1d2, PRR7_other_d1d2)

names(PRR7_bind_d1d2)[3] <- "PRR7"

write_csv(PRR7_bind_d1d2, './01_tidy_data/PRR7_bind_d1d2.csv')

# LUX amp_gain_high (58 obs) + LUX amp_gain_medium (86 obs) + LUX amp_lose_high (103 obs) + LUX amp_lose_medium (42 obs) + LUX amp_other (103 obs): total equals 392 obs
LUX_bind_d1d2 <- bind_rows(LUX_gain_high_d1d2, LUX_gain_medium_d1d2, LUX_lose_high_d1d2, LUX_lose_medium_d1d2, LUX_other_d1d2)

names(LUX_bind_d1d2)[3] <- "LUX"

write_csv(LUX_bind_d1d2, './01_tidy_data/LUX_bind_d1d2.csv')

# ELF3 amp_gain_high (19 obs) + ELF3 amp_gain_medium (19 obs) + ELF3 amp_lose_high (44 obs) + ELF3 amp_lose_medium (11 obs) + ELF3 amp_other (27 obs): total equals 120 obs
ELF3_bind_d1d2 <- bind_rows(ELF3_gain_high_d1d2, ELF3_gain_medium_d1d2, ELF3_lose_high_d1d2, ELF3_lose_medium_d1d2, ELF3_other_d1d2)

names(ELF3_bind_d1d2)[3] <- "ELF3"

write_csv(ELF3_bind_d1d2, './01_tidy_data/ELF3_bind_d1d2.csv')

# ELF4 amp_gain_high (8 obs) + ELF4 amp_gain_medium (3 obs) + ELF4 amp_lose_high (12 obs) + ELF4 amp_lose_medium (2 obs) + ELF4 amp_other (5 obs): total equals 30 obs
ELF4_bind_d1d2 <- bind_rows(ELF4_gain_high_d1d2, ELF4_gain_medium_d1d2, ELF4_lose_high_d1d2, ELF4_lose_medium_d1d2, ELF4_other_d1d2)

names(ELF4_bind_d1d2)[3] <- "ELF4"

write_csv(ELF4_bind_d1d2, './01_tidy_data/ELF4_bind_d1d2.csv')

# * d1d5----
# append all the LHY targets:
# LHY amp_gain_high (0 obs) + LHY amp_gain_medium (16 obs) + LHY amp_lose_high (87 obs) + LHY amp_lose_medium (53 obs) + LHY amp_other (144 obs): total equals 300 obs
LHY_bind_d1d5 <- bind_rows(LHY_gain_high_d1d5, LHY_gain_medium_d1d5, LHY_lose_high_d1d5, LHY_lose_medium_d1d5, LHY_other_d1d5)

names(LHY_bind_d1d5)[3] <- "LHY"

write_csv(LHY_bind_d1d5, './01_tidy_data/LHY_bind_d1d5.csv')

# CCA1-Nagel amp_gain_high (0 obs) + CCA1-Nagel amp_gain_medium (26 obs) + CCA1-Nagel amp_lose_high (166 obs) + CCA1-Nagel amp_lose_medium (82 obs) + CCA1-Nagel amp_other (210 obs): total equals 484 obs
CCA1_nagel_bind_d1d5 <- bind_rows(CCA1_nagel_gain_high_d1d5, CCA1_nagel_gain_medium_d1d5, CCA1_nagel_lose_high_d1d5, CCA1_nagel_lose_medium_d1d5, CCA1_nagel_other_d1d5)

names(CCA1_nagel_bind_d1d5)[3] <- "CCA1 Nagel"

write_csv(CCA1_nagel_bind_d1d5, './01_tidy_data/CCA1_nagel_bind_d1d5.csv')

# CCA1-Kamioka amp_gain_high (1 obs) + CCA1-Kamioka amp_gain_medium (6 obs) + CCA1-Kamioka amp_lose_high (53 obs) + CCA1-Kamioka amp_lose_medium (31 obs) + CCA1-Kamioka amp_other (88 obs): total equals 179 obs
CCA1_kamioka_bind_d1d5 <- bind_rows(CCA1_kamioka_gain_high_d1d5, CCA1_kamioka_gain_medium_d1d5, CCA1_kamioka_lose_high_d1d5, CCA1_kamioka_lose_medium_d1d5, CCA1_kamioka_other_d1d5)

names(CCA1_kamioka_bind_d1d5)[3] <- "CCA1 Kamioka"

write_csv(CCA1_kamioka_bind_d1d5, './01_tidy_data/CCA1_kamioka_bind_d1d5.csv')

# CCA1-Nagel-Kamioka amp_gain_high (0 obs) + CCA1-Nagel-Kamioka amp_gain_medium (4 obs) + CCA1-Nagel-Kamioka amp_lose_high (21 obs) + CCA1-Nagel-Kamioka amp_lose_medium (28 obs) + CCA1-Nagel-Kamioka amp_other (64 obs): total equals 117 obs
CCA1_nagel_kamioka_bind_d1d5 <- bind_rows(CCA1_nagel_kamioka_gain_high_d1d5, CCA1_nagel_kamioka_gain_medium_d1d5, CCA1_nagel_kamioka_lose_high_d1d5, CCA1_nagel_kamioka_lose_medium_d1d5, CCA1_nagel_kamioka_other_d1d5)

names(CCA1_nagel_kamioka_bind_d1d5)[3] <- "CCA1 Nagel-Kamioka"

write_csv(CCA1_nagel_kamioka_bind_d1d5, './01_tidy_data/CCA1_nagel_kamioka_bind_d1d5.csv')

# TOC1 amp_gain_high (0 obs) + TOC1 amp_gain_medium (13 obs) + TOC1 amp_lose_high (95 obs) + TOC1 amp_lose_medium (30 obs) + TOC1 amp_other (110 obs): total equals 248 obs
TOC1_bind_d1d5 <- bind_rows(TOC1_gain_high_d1d5, TOC1_gain_medium_d1d5, TOC1_lose_high_d1d5, TOC1_lose_medium_d1d5, TOC1_other_d1d5)

names(TOC1_bind_d1d5)[3] <- "TOC1"

write_csv(TOC1_bind_d1d5, './01_tidy_data/TOC1_bind_d1d5.csv')

# PRR5 amp_gain_high (0 obs) + PRR5 amp_gain_medium (1 obs) + PRR5 amp_lose_high (6 obs) + PRR5 amp_lose_medium (5 obs) + PRR5 amp_other (31 obs): total equals 43 obs
PRR5_bind_d1d5 <- bind_rows(PRR5_gain_high_d1d5, PRR5_gain_medium_d1d5, PRR5_lose_high_d1d5, PRR5_lose_medium_d1d5, PRR5_other_d1d5)

names(PRR5_bind_d1d5)[3] <- "PRR5"

write_csv(PRR5_bind_d1d5, './01_tidy_data/PRR5_bind_d1d5.csv')

# PRR7 amp_gain_high (0 obs) + PRR7 amp_gain_medium (3 obs) + PRR7 amp_lose_high (18 obs) + PRR7 amp_lose_medium (6 obs) + PRR7 amp_other (33 obs): total equals 60 obs
PRR7_bind_d1d5 <- bind_rows(PRR7_gain_high_d1d5, PRR7_gain_medium_d1d5, PRR7_lose_high_d1d5, PRR7_lose_medium_d1d5, PRR7_other_d1d5)

names(PRR7_bind_d1d5)[3] <- "PRR7"

write_csv(PRR7_bind_d1d5, './01_tidy_data/PRR7_bind_d1d5.csv')

# LUX amp_gain_high (3 obs) + LUX amp_gain_medium (21 obs) + LUX amp_lose_high (181 obs) + LUX amp_lose_medium (43 obs) + LUX amp_other (116 obs): total equals 364 obs
LUX_bind_d1d5 <- bind_rows(LUX_gain_high_d1d5, LUX_gain_medium_d1d5, LUX_lose_high_d1d5, LUX_lose_medium_d1d5, LUX_other_d1d5)

names(LUX_bind_d1d5)[3] <- "LUX"

write_csv(LUX_bind_d1d5, './01_tidy_data/LUX_bind_d1d5.csv')

# ELF3 amp_gain_high (0 obs) + ELF3 amp_gain_medium (7 obs) + ELF3 amp_lose_high (59 obs) + ELF3 amp_lose_medium (16 obs) + ELF3 amp_other (31 obs): total equals 113 obs
ELF3_bind_d1d5 <- bind_rows(ELF3_gain_high_d1d5, ELF3_gain_medium_d1d5, ELF3_lose_high_d1d5, ELF3_lose_medium_d1d5, ELF3_other_d1d5)

names(ELF3_bind_d1d5)[3] <- "ELF3"

write_csv(ELF3_bind_d1d5, './01_tidy_data/ELF3_bind_d1d5.csv')

# ELF4 amp_gain_high (0 obs) + ELF4 amp_gain_medium (0 obs) + ELF4 amp_lose_high (14 obs) + ELF4 amp_lose_medium (3 obs) + ELF4 amp_other (10 obs): total equals 27 obs
ELF4_bind_d1d5 <- bind_rows(ELF4_gain_high_d1d5, ELF4_gain_medium_d1d5, ELF4_lose_high_d1d5, ELF4_lose_medium_d1d5, ELF4_other_d1d5)

names(ELF4_bind_d1d5)[3] <- "ELF4"

write_csv(ELF4_bind_d1d5, './01_tidy_data/ELF4_bind_d1d5.csv')

summarise_targets <- function(df, col_str, clock_id){
  
  summary <- df %>% 
    group_by({{col_str}}) %>% 
    dplyr::summarise(n=n()) %>% 
    mutate(freq = (n/sum(n) *100)) %>% 
    mutate(clock = {{clock_id}})
  
  return(summary)
  
}

# d1-d2----
LHY_d1d2_summary <- summarise_targets(LHY_bind_d1d2, type, 'LHY')
CCA1_nagel_d1d2_summary <- summarise_targets(CCA1_nagel_bind_d1d2, type, 'CCA1 Nagel')
CCA1_kamioka_d1d2_summary <- summarise_targets(CCA1_kamioka_bind_d1d2, type, 'CCA1 Kamioka')
CCA1_nagel_kamioka_d1d2_summary <- summarise_targets(CCA1_nagel_kamioka_bind_d1d2, type, 'CCA1')
TOC1_d1d2_summary <- summarise_targets(TOC1_bind_d1d2, type, 'TOC1')
PRR5_d1d2_summary <- summarise_targets(PRR5_bind_d1d2, type, 'PRR5')
PRR7_d1d2_summary <- summarise_targets(PRR7_bind_d1d2, type, 'PRR7')
LUX_d1d2_summary <- summarise_targets(LUX_bind_d1d2, type, 'LUX')
ELF3_d1d2_summary <- summarise_targets(ELF3_bind_d1d2, type, 'ELF3')
ELF4_d1d2_summary <- summarise_targets(ELF4_bind_d1d2, type, 'ELF4')

clock_d1_d2 <- bind_rows(LHY_d1d2_summary, CCA1_nagel_kamioka_d1d2_summary, TOC1_d1d2_summary,
                         PRR5_d1d2_summary, PRR7_d1d2_summary, LUX_d1d2_summary, ELF3_d1d2_summary,
                         ELF4_d1d2_summary) 

# d1-d5----
LHY_d1d5_summary <- summarise_targets(LHY_bind_d1d5, type, 'LHY')
CCA1_nagel_d1d5_summary <- summarise_targets(CCA1_nagel_bind_d1d5, type, 'CCA1 Nagel')
CCA1_kamioka_d1d5_summary <- summarise_targets(CCA1_kamioka_bind_d1d5, type, 'CCA1 Kamioka')
CCA1_nagel_kamioka_d1d5_summary <- summarise_targets(CCA1_nagel_kamioka_bind_d1d5, type, 'CCA1')
TOC1_d1d5_summary <- summarise_targets(TOC1_bind_d1d5, type, 'TOC1')
PRR5_d1d5_summary <- summarise_targets(PRR5_bind_d1d5, type, 'PRR5')
PRR7_d1d5_summary <- summarise_targets(PRR7_bind_d1d5, type, 'PRR7')
LUX_d1d5_summary <- summarise_targets(LUX_bind_d1d5, type, 'LUX')
ELF3_d1d5_summary <- summarise_targets(ELF3_bind_d1d5, type, 'ELF3')
ELF4_d1d5_summary <- summarise_targets(ELF4_bind_d1d5, type, 'ELF4')

clock_d1_d5 <- bind_rows(LHY_d1d5_summary, CCA1_nagel_kamioka_d1d5_summary, TOC1_d1d5_summary,
                         PRR5_d1d5_summary, PRR7_d1d5_summary, LUX_d1d5_summary, ELF3_d1d5_summary,
                         ELF4_d1d5_summary)

plot_clock_d1d2 <- clock_d1_d2 %>% 
  mutate(type = factor(type, levels = c('gain_high_d1_d2', 'gain_medium_d1_d2', 'other_d1_d2', 'lose_medium_d1_d2', 'lose_high_d1_d2')),
         clock = factor(clock, levels = c('CCA1', 'LHY', 'TOC1', 'PRR5', 'PRR7', 'LUX', 'ELF3', 'ELF4'))) %>% 
  ggplot(aes(x = clock, y = freq, fill = type)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_fill_manual(values = c('#31a354', '#74c476', '#cccccc', '#fb6a4a', '#de2d26'), labels = c('gain - high', 'gain - medium', 'other', 'lose - medium', 'lose - high')) +
  geom_bar(position = 'stack', stat = 'identity') +
  ggpubr::theme_pubr() +
  theme(legend.position = "right") +
  theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1)) +
  labs(fill = 'Amplitude', y = 'Proportion (%)', x = '') +
  ggtitle("Amplitude patterns for clock targets",
          subtitle = "Compare day 1 (20C steady state) with day 2 (to 4C transient-cooling)") 

plot_clock_d1d2

ggsave('./03_plots/plot_clock_d1d2.png', dpi = 300, height = 6, width = 6, units = 'in')

plot_clock_d1d5 <- clock_d1_d5 %>% 
  mutate(type = factor(type, levels = c('gain_high_d1_d5', 'gain_medium_d1_d5', 'other_d1_d5', 'lose_medium_d1_d5', 'lose_high_d1_d5')),
         clock = factor(clock, levels = c('CCA1', 'LHY', 'TOC1', 'PRR5', 'PRR7', 'LUX', 'ELF3', 'ELF4'))) %>% 
  ggplot(aes(x = clock, y = freq, fill = type)) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80, 100)) +
  scale_fill_manual(values = c('#31a354', '#74c476', '#cccccc', '#fb6a4a', '#de2d26'), labels = c('gain - high', 'gain - medium', 'other', 'lose - medium', 'lose - high')) +
  geom_bar(position = 'stack', stat = 'identity') +
  ggpubr::theme_pubr() +
  theme(legend.position = "right") +
  theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1)) +
  labs(fill = 'Amplitude', y = 'Proportion (%)', x = '') +
  ggtitle("Amplitude patterns for clock targets",
          subtitle = "Compare day 1 (20C steady state) with Day 5 (4C steady state)") 

plot_clock_d1d5

ggsave('./03_plots/plot_clock_d1d5.png', dpi = 300, height = 6, width = 6, units = 'in')

# circular barplot----
# d1d2----

gain_amp_high_d1d2 <- amp_gain_high_clusters_d1_d2 %>% 
  mutate(group = 'gain high')

gain_amp_med_d1d2 <- amp_gain_medium_clusters_d1_d2 %>% 
  mutate(group = 'gain medium')

lose_amp_high_d1d2 <- amp_lose_high_clusters_d1_d2 %>% 
  mutate(group = 'lose high')

lose_amp_med_d1d2 <- amp_lose_medium_clusters_d1_d2 %>% 
  mutate(group = 'lose medium')

other_amp_d1d2 <- amp_other_clusters_d1_d2 %>% 
  mutate(group = 'other')


# list all clusters with their amplitude profile description
clusters_d1d2 <- bind_rows(gain_amp_high_d1d2,
                           gain_amp_med_d1d2,
                           lose_amp_high_d1d2,
                           lose_amp_med_d1d2,
                           other_amp_d1d2) %>% 
  mutate_at(1, as.numeric) %>%
  arrange(cluster)

get_CCGs_clusters <- function(df1, col_str1, col_str2, df2, clock_id, type1, type2, type3, type4){
  
  tidy_clock_summary <- df1 %>% 
    group_by({{col_str1}}, {{col_str2}}) %>% 
    dplyr::summarise(n=n()) %>%
    mutate_at(1, as.numeric) %>%
    arrange(cluster) %>%
    ungroup() %>% 
    mutate(group = case_when(type == {{type1}} ~ 'gain high',
                             type == {{type2}} ~ 'gain medium',
                             type == {{type3}} ~ 'lose high',
                             type == {{type4}} ~ 'lose medium',
                             TRUE ~ 'other')) %>%
    select(-type)
  
  complete_clusters <- df2 %>%
    left_join(tidy_clock_summary) %>% 
    mutate_at(3, ~replace_na(.,0)) %>% 
    dplyr::rename({{clock_id}} := n)
  
  return(complete_clusters)
  
}

# d1-d2----
LHY_CCG_cl_d1d2 <- get_CCGs_clusters(LHY_bind_d1d2, cluster, type, clusters_d1d2, LHY, 'gain_high_d1_d2', 'gain_medium_d1_d2', 'lose_high_d1_d2', 'lose_medium_d1_d2')
CCA1_CCG_cl_d1d2 <- get_CCGs_clusters(CCA1_nagel_kamioka_bind_d1d2, cluster, type, clusters_d1d2, CCA1, 'gain_high_d1_d2', 'gain_medium_d1_d2', 'lose_high_d1_d2', 'lose_medium_d1_d2')
TOC1_CCG_cl_d1d2 <- get_CCGs_clusters(TOC1_bind_d1d2, cluster, type, clusters_d1d2, TOC1, 'gain_high_d1_d2', 'gain_medium_d1_d2', 'lose_high_d1_d2', 'lose_medium_d1_d2')
PRR5_CCG_cl_d1d2 <- get_CCGs_clusters(PRR5_bind_d1d2, cluster, type, clusters_d1d2, PRR5, 'gain_high_d1_d2', 'gain_medium_d1_d2', 'lose_high_d1_d2', 'lose_medium_d1_d2')
PRR7_CCG_cl_d1d2 <- get_CCGs_clusters(PRR7_bind_d1d2, cluster, type, clusters_d1d2, PRR7, 'gain_high_d1_d2', 'gain_medium_d1_d2', 'lose_high_d1_d2', 'lose_medium_d1_d2')
LUX_CCG_cl_d1d2 <- get_CCGs_clusters(LUX_bind_d1d2, cluster, type, clusters_d1d2, LUX, 'gain_high_d1_d2', 'gain_medium_d1_d2', 'lose_high_d1_d2', 'lose_medium_d1_d2')
ELF3_CCG_cl_d1d2 <- get_CCGs_clusters(ELF3_bind_d1d2, cluster, type, clusters_d1d2, ELF3, 'gain_high_d1_d2', 'gain_medium_d1_d2', 'lose_high_d1_d2', 'lose_medium_d1_d2')
ELF4_CCG_cl_d1d2 <- get_CCGs_clusters(ELF4_bind_d1d2, cluster, type, clusters_d1d2, ELF4, 'gain_high_d1_d2', 'gain_medium_d1_d2', 'lose_high_d1_d2', 'lose_medium_d1_d2')

clock_clusters_d1d2 <- purrr::reduce(list(LHY_CCG_cl_d1d2, CCA1_CCG_cl_d1d2, TOC1_CCG_cl_d1d2,
                                          PRR5_CCG_cl_d1d2, PRR7_CCG_cl_d1d2, LUX_CCG_cl_d1d2,
                                          ELF3_CCG_cl_d1d2, ELF4_CCG_cl_d1d2), dplyr::left_join) %>% 
  rowwise(cluster) %>%
  dplyr::mutate(sum = sum(c_across(LHY:ELF4))) %>% 
  relocate(sum, .before = LHY) %>% 
  mutate(group = factor(group, levels = c('gain high', 'gain medium', 'lose high', 'lose medium', 'other')))

circbar_d1d2 <- clock_clusters_d1d2 %>% 
  pivot_longer(cols = LHY:ELF4, values_to = 'number')

#https://www.r-graph-gallery.com/295-basic-circular-barplot.html
#https://www.r-graph-gallery.com/296-add-labels-to-circular-barplot
#https://www.r-graph-gallery.com/299-circular-stacked-barplot.html

# Set a number of 'empty bar' to add at the end of each group
circbar_d1d2_empty_bar <- 2

circbar_d1d2_nOBsType <- nlevels(as.factor(circbar_d1d2$name))

circbar_d1d2_to_add <- data.frame( matrix(NA, circbar_d1d2_empty_bar*nlevels(circbar_d1d2$group)*circbar_d1d2_nOBsType, ncol(circbar_d1d2)) )

colnames(circbar_d1d2_to_add) <- colnames(circbar_d1d2)

circbar_d1d2_to_add$group <- rep(levels(circbar_d1d2$group), each=circbar_d1d2_empty_bar*circbar_d1d2_nOBsType )

circbar_d1d2 <- rbind(circbar_d1d2, circbar_d1d2_to_add)

circbar_d1d2_gc <- circbar_d1d2 %>% arrange(group, cluster)

circbar_d1d2_gs <- circbar_d1d2_gc %>% arrange(group, sum)

circbar_d1d2_gs$id <- rep( seq(1, nrow(circbar_d1d2_gs)/circbar_d1d2_nOBsType) , each=circbar_d1d2_nOBsType)

circbar_d1d2_gs$name <- factor(circbar_d1d2_gs$name, levels = c("CCA1", "LHY", "TOC1", "PRR5", "PRR7", "LUX", "ELF3", "ELF4"))

# Get the name and the y position of each label
label_data_circbar_d1d2<- circbar_d1d2_gs %>% dplyr::group_by(id, cluster) %>% dplyr::summarize(tot=sum(number))

number_of_bar_circbar_d1d2 <- nrow(label_data_circbar_d1d2)

angle_circbar_d1d2 <- 90 - 360 * (label_data_circbar_d1d2$id-0.5) /number_of_bar_circbar_d1d2

label_data_circbar_d1d2$hjust <- ifelse(angle_circbar_d1d2 < -90, 1, 0)

label_data_circbar_d1d2$angle <- ifelse(angle_circbar_d1d2 < -90, angle_circbar_d1d2+180, angle_circbar_d1d2)

# prepare a data frame for base lines
base_data_circbar_d1d2 <- circbar_d1d2_gs %>% dplyr::group_by(group) %>% dplyr::summarize(start=min(id), end=max(id) - circbar_d1d2_empty_bar) %>% dplyr::rowwise() %>% dplyr::mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data_circbar_d1d2 <- base_data_circbar_d1d2

grid_data_circbar_d1d2$end <- grid_data_circbar_d1d2$end[ c( nrow(grid_data_circbar_d1d2), 1:nrow(grid_data_circbar_d1d2)-1)] + 1

grid_data_circbar_d1d2$start <- grid_data_circbar_d1d2$start - 1

grid_data_circbar_d1d2 <- grid_data_circbar_d1d2[-1,]
circbar_d1d2_plot

# Make the plot
circbar_d1d2_plot <- ggplot(circbar_d1d2_gs) +
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=number, fill=name), stat="identity", alpha=0.75) +
  scale_fill_manual (values=c("#1a9850", "#a6d96a", "#4575b4", "#fdae61", "#f46d43", "#542788", "#636363", "#cccccc")) +
  
  # Add a val=150/100/50/0 lines
  geom_segment(data= grid_data_circbar_d1d2, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey30", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data= grid_data_circbar_d1d2, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey30", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data= grid_data_circbar_d1d2, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey30", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data= grid_data_circbar_d1d2, aes(x = end, y = 150, xend = start, yend = 150), colour = "grey30", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 150/100/50/0 line
  ggplot2::annotate("text", x = rep(max(circbar_d1d2_gs$id),4), y = c(0, 50, 100, 150), label = c("0", "50", "100", "150") , color="grey30", size=5 , angle=0, fontface="bold", hjust=1) +
  ylim(-150,max(20+label_data_circbar_d1d2$tot, na.rm=T)) +
  theme_minimal() +
  theme(legend.position=c(0.25,0.8),
        legend.text = element_text(color = "black", size = 9),
        legend.title= element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm")) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data_circbar_d1d2, aes(x=id, y=tot+5, label=cluster, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4, angle= label_data_circbar_d1d2$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data_circbar_d1d2, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data_circbar_d1d2, aes(x = title, y = -15, label=group), hjust=c(1,1,0.5,0,0), vjust=c(0,0,-1,0,1), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

ggsave(circbar_d1d2_plot, file="./03_plots/circbar_d1d2_plot.png", width=8, height=8, units="in",dpi=200)

# d1-d5----

gain_amp_high_d1d5 <- amp_gain_high_clusters_d1_d5 %>% 
  mutate(group = 'gain high')

gain_amp_med_d1d5 <- amp_gain_medium_clusters_d1_d5 %>% 
  mutate(group = 'gain medium')

lose_amp_high_d1d5 <- amp_lose_high_clusters_d1_d5 %>% 
  mutate(group = 'lose high')

lose_amp_med_d1d5 <- amp_lose_medium_clusters_d1_d5 %>% 
  mutate(group = 'lose medium')

other_amp_d1d5 <- amp_other_clusters_d1_d5 %>% 
  mutate(group = 'other')

clusters_d1d5 <- bind_rows(gain_amp_high_d1d5,
                           gain_amp_med_d1d5,
                           lose_amp_high_d1d5,
                           lose_amp_med_d1d5,
                           other_amp_d1d5) %>% 
  mutate_at(1, as.numeric) %>%
  arrange(cluster)

LHY_CCG_cl_d1d5 <- get_CCGs_clusters(LHY_bind_d1d5, cluster, type, clusters_d1d5, LHY, 'gain_high_d1_d5', 'gain_medium_d1_d5', 'lose_high_d1_d5', 'lose_medium_d1_d5')
CCA1_CCG_cl_d1d5 <- get_CCGs_clusters(CCA1_nagel_kamioka_bind_d1d5, cluster, type, clusters_d1d5, CCA1, 'gain_high_d1_d5', 'gain_medium_d1_d5', 'lose_high_d1_d5', 'lose_medium_d1_d5')
TOC1_CCG_cl_d1d5 <- get_CCGs_clusters(TOC1_bind_d1d5, cluster, type, clusters_d1d5, TOC1, 'gain_high_d1_d5', 'gain_medium_d1_d5', 'lose_high_d1_d5', 'lose_medium_d1_d5')
PRR5_CCG_cl_d1d5 <- get_CCGs_clusters(PRR5_bind_d1d5, cluster, type, clusters_d1d5, PRR5, 'gain_high_d1_d5', 'gain_medium_d1_d5', 'lose_high_d1_d5', 'lose_medium_d1_d5')
PRR7_CCG_cl_d1d5 <- get_CCGs_clusters(PRR7_bind_d1d5, cluster, type, clusters_d1d5, PRR7, 'gain_high_d1_d5', 'gain_medium_d1_d5', 'lose_high_d1_d5', 'lose_medium_d1_d5')
LUX_CCG_cl_d1d5 <- get_CCGs_clusters(LUX_bind_d1d5, cluster, type, clusters_d1d5, LUX, 'gain_high_d1_d5', 'gain_medium_d1_d5', 'lose_high_d1_d5', 'lose_medium_d1_d5')
ELF3_CCG_cl_d1d5 <- get_CCGs_clusters(ELF3_bind_d1d5, cluster, type, clusters_d1d5, ELF3, 'gain_high_d1_d5', 'gain_medium_d1_d5', 'lose_high_d1_d5', 'lose_medium_d1_d5')
ELF4_CCG_cl_d1d5 <- get_CCGs_clusters(ELF4_bind_d1d5, cluster, type, clusters_d1d5, ELF4, 'gain_high_d1_d5', 'gain_medium_d1_d5', 'lose_high_d1_d5', 'lose_medium_d1_d5')

clock_clusters_d1d5 <- purrr::reduce(list(LHY_CCG_cl_d1d5, CCA1_CCG_cl_d1d5, TOC1_CCG_cl_d1d5,
                                          PRR5_CCG_cl_d1d5, PRR7_CCG_cl_d1d5, LUX_CCG_cl_d1d5,
                                          ELF3_CCG_cl_d1d5, ELF4_CCG_cl_d1d5), dplyr::left_join) %>% 
  rowwise(cluster) %>%
  dplyr::mutate(sum = sum(c_across(LHY:ELF4))) %>% 
  relocate(sum, .before = LHY) %>% 
  mutate(group = factor(group, levels = c('gain high', 'gain medium', 'lose high', 'lose medium', 'other')))

circbar_d1d5 <- clock_clusters_d1d5 %>% 
  pivot_longer(cols = LHY:ELF4, values_to = 'number') 

# Set a number of 'empty bar' to add at the end of each group
circbar_d1d5_empty_bar <- 2

circbar_d1d5_nOBsType <- nlevels(as.factor(circbar_d1d5$name))

circbar_d1d5_to_add <- data.frame( matrix(NA, circbar_d1d5_empty_bar*nlevels(circbar_d1d5$group)*circbar_d1d5_nOBsType, ncol(circbar_d1d5)) )

colnames(circbar_d1d5_to_add) <- colnames(circbar_d1d5)

circbar_d1d5_to_add$group <- rep(levels(circbar_d1d5$group), each=circbar_d1d5_empty_bar*circbar_d1d5_nOBsType )

circbar_d1d5 <- rbind(circbar_d1d5, circbar_d1d5_to_add)

circbar_d1d5_gc <- circbar_d1d5 %>% arrange(group, cluster)

circbar_d1d5_gs <- circbar_d1d5_gc %>% arrange(group, sum)

circbar_d1d5_gs$id <- rep( seq(1, nrow(circbar_d1d5_gs)/circbar_d1d5_nOBsType) , each=circbar_d1d5_nOBsType)

circbar_d1d5_gs$name <- factor(circbar_d1d5_gs$name, levels = c("CCA1", "LHY", "TOC1", "PRR5", "PRR7", "LUX", "ELF3", "ELF4"))

# Get the name and the y position of each label
label_data_circbar_d1d5<- circbar_d1d5_gs %>% dplyr::group_by(id, cluster) %>% dplyr::summarize(tot=sum(number))

number_of_bar_circbar_d1d5 <- nrow(label_data_circbar_d1d5)

angle_circbar_d1d5 <- 90 - 360 * (label_data_circbar_d1d5$id-0.5) /number_of_bar_circbar_d1d5

label_data_circbar_d1d5$hjust <- ifelse(angle_circbar_d1d5 < -90, 1, 0)

label_data_circbar_d1d5$angle <- ifelse(angle_circbar_d1d5 < -90, angle_circbar_d1d5+180, angle_circbar_d1d5)

# prepare a data frame for base lines
base_data_circbar_d1d5 <- circbar_d1d5_gs %>% dplyr::group_by(group) %>% dplyr::summarize(start=min(id), end=max(id) - circbar_d1d5_empty_bar) %>% dplyr::rowwise() %>% dplyr::mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data_circbar_d1d5 <- base_data_circbar_d1d5

grid_data_circbar_d1d5$end <- grid_data_circbar_d1d5$end[ c( nrow(grid_data_circbar_d1d5), 1:nrow(grid_data_circbar_d1d5)-1)] + 1

grid_data_circbar_d1d5$start <- grid_data_circbar_d1d5$start - 1

grid_data_circbar_d1d5 <- grid_data_circbar_d1d5[-1,]

circbar_d1d5_plot

# Make the plot
circbar_d1d5_plot <- ggplot(circbar_d1d5_gs) +
  
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=number, fill=name), stat="identity", alpha=0.75) +
  scale_fill_manual (values=c("#1a9850", "#a6d96a", "#4575b4", "#fdae61", "#f46d43", "#542788", "#636363", "#cccccc")) +
  
  # Add a val=150/100/50/0 lines
  geom_segment(data= grid_data_circbar_d1d5, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey30", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data= grid_data_circbar_d1d5, aes(x = end, y = 50, xend = start, yend = 50), colour = "grey30", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data= grid_data_circbar_d1d5, aes(x = end, y = 100, xend = start, yend = 100), colour = "grey30", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data= grid_data_circbar_d1d5, aes(x = end, y = 150, xend = start, yend = 150), colour = "grey30", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 150/100/50/0 line
  ggplot2::annotate("text", x = rep(max(circbar_d1d5_gs$id),4), y = c(0, 50, 100, 150), label = c("0", "50", "100", "150") , color="grey30", size=5 , angle=0, fontface="bold", hjust=1) +
  ylim(-150,max(20+label_data_circbar_d1d5$tot, na.rm=T)) +
  theme_minimal() +
  theme(legend.position=c(0.25,0.8),
        legend.text = element_text(color = "black", size = 9),
        legend.title= element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-1,4), "cm")) +
  coord_polar() +
  
  # Add labels on top of each bar
  geom_text(data=label_data_circbar_d1d5, aes(x=id, y=tot+5, label=cluster, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=4, angle= label_data_circbar_d1d5$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data_circbar_d1d5, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data_circbar_d1d5, aes(x = title, y = -15, label=group), hjust=c(0.5,1,1,0,0), vjust=c(0.5,0,-1,0,1), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

ggsave(circbar_d1d5_plot, file="./03_plots/circbar_d1d5_plot.png", width=8, height=8, units="in",dpi=200)

# Odds Ratios----

# 7302 gene loci grouped in 75 clusters (0-74); error in row 5772
TF_network_clusters <- read_csv('./00_raw_data/TF Network Cluster Nov2018 long format.csv') %>%
  filter(!row_number() %in% 5772)

# count of gene loci in each TF cluster
TF_network_clusters_count <-  TF_network_clusters %>% 
  dplyr::group_by(cluster) %>%
  dplyr::summarise(cluster_number = n()) %>% 
  filter(!cluster == 0)

LHY_CCG_cl_d1d2_fishers <- LHY_CCG_cl_d1d2 %>% 
  dplyr::select(-group) %>% 
  inner_join(TF_network_clusters_count) %>% 
  mutate(fishers_col2 = sum(LHY) - LHY,
         fishers_col4 = sum(cluster_number) - cluster_number) %>% 
  dplyr::select(LHY, cluster_number, fishers_col2, fishers_col4) %>% 
  dplyr::rename(fishers_col1 = LHY,
                fishers_col3 = cluster_number) %>% 
  relocate(fishers_col2, .after = fishers_col1)

fishers_prep<- function(df1, col_str1, col_str2, clock_id, col_str3, col_str4){
  
  prep1 <- df1 %>% 
    dplyr::select(-group) %>% 
    inner_join(TF_network_clusters_count) %>%
    dplyr::mutate({{col_str1}} := sum({{clock_id}}) - {{clock_id}},
           {{col_str2}} := sum(cluster_number) - cluster_number)
  
  prep2 <- prep1 %>% 
    dplyr::select({{clock_id}}, cluster_number, {{col_str1}}, {{col_str2}}) %>% 
    dplyr::rename({{col_str3}} := {{clock_id}},
                  {{col_str4}} := cluster_number) %>% 
    relocate({{col_str1}}, .after = {{col_str3}})
  
  return(prep2)
  
  
}

#d1-d2----

LHY_d1d2_fishers_prep <- fishers_prep(LHY_CCG_cl_d1d2, fishers_col2, fishers_col4, LHY, fishers_col1, fishers_col3)
CCA1_d1d2_fishers_prep <- fishers_prep(CCA1_CCG_cl_d1d2, fishers_col2, fishers_col4, CCA1, fishers_col1, fishers_col3)
TOC1_d1d2_fishers_prep <- fishers_prep(TOC1_CCG_cl_d1d2, fishers_col2, fishers_col4, TOC1, fishers_col1, fishers_col3)
PRR5_d1d2_fishers_prep <- fishers_prep(PRR5_CCG_cl_d1d2, fishers_col2, fishers_col4, PRR5, fishers_col1, fishers_col3)
PRR7_d1d2_fishers_prep <- fishers_prep(PRR7_CCG_cl_d1d2, fishers_col2, fishers_col4, PRR7, fishers_col1, fishers_col3)
LUX_d1d2_fishers_prep <- fishers_prep(LUX_CCG_cl_d1d2, fishers_col2, fishers_col4, LUX, fishers_col1, fishers_col3)
ELF3_d1d2_fishers_prep <- fishers_prep(ELF3_CCG_cl_d1d2, fishers_col2, fishers_col4, ELF3, fishers_col1, fishers_col3)
ELF4_d1d2_fishers_prep <- fishers_prep(ELF4_CCG_cl_d1d2, fishers_col2, fishers_col4, ELF4, fishers_col1, fishers_col3)


#d1-d5----

LHY_d1d5_fishers_prep <- fishers_prep(LHY_CCG_cl_d1d5, fishers_col2, fishers_col4, LHY, fishers_col1, fishers_col3)
CCA1_d1d5_fishers_prep <- fishers_prep(CCA1_CCG_cl_d1d5, fishers_col2, fishers_col4, CCA1, fishers_col1, fishers_col3)
TOC1_d1d5_fishers_prep <- fishers_prep(TOC1_CCG_cl_d1d5, fishers_col2, fishers_col4, TOC1, fishers_col1, fishers_col3)
PRR5_d1d5_fishers_prep <- fishers_prep(PRR5_CCG_cl_d1d5, fishers_col2, fishers_col4, PRR5, fishers_col1, fishers_col3)
PRR7_d1d5_fishers_prep <- fishers_prep(PRR7_CCG_cl_d1d5, fishers_col2, fishers_col4, PRR7, fishers_col1, fishers_col3)
LUX_d1d5_fishers_prep <- fishers_prep(LUX_CCG_cl_d1d5, fishers_col2, fishers_col4, LUX, fishers_col1, fishers_col3)
ELF3_d1d5_fishers_prep <- fishers_prep(ELF3_CCG_cl_d1d5, fishers_col2, fishers_col4, ELF3, fishers_col1, fishers_col3)
ELF4_d1d5_fishers_prep <- fishers_prep(ELF4_CCG_cl_d1d5, fishers_col2, fishers_col4, ELF4, fishers_col1, fishers_col3)

get_odds <- function(df){
  
  odds<- df %>%
    data.frame(apply(., 1, function(x) fisher.test(matrix(x, nr=2), alternative="greater")$estimate))
  
  return(odds)
    
}

#d1-d2----
LHY_d1d2_fishers <- get_odds(LHY_d1d2_fishers_prep)
colnames(LHY_d1d2_fishers)[5] <- "LHY"
LHY_d1d2_fishers <- LHY_d1d2_fishers %>% 
  dplyr::select(LHY) %>% 
  round(2)

CCA1_d1d2_fishers <- get_odds(CCA1_d1d2_fishers_prep)
colnames(CCA1_d1d2_fishers)[5] <- "CCA1"
CCA1_d1d2_fishers <- CCA1_d1d2_fishers %>% 
  dplyr::select(CCA1) %>% 
  round(2)

TOC1_d1d2_fishers <- get_odds(TOC1_d1d2_fishers_prep)
colnames(TOC1_d1d2_fishers)[5] <- "TOC1"
TOC1_d1d2_fishers <- TOC1_d1d2_fishers %>% 
  dplyr::select(TOC1) %>% 
  round(2)

PRR5_d1d2_fishers <- get_odds(PRR5_d1d2_fishers_prep)
colnames(PRR5_d1d2_fishers)[5] <- "PRR5"
PRR5_d1d2_fishers <- PRR5_d1d2_fishers %>% 
  dplyr::select(PRR5) %>% 
  round(2)

PRR7_d1d2_fishers <- get_odds(PRR7_d1d2_fishers_prep)
colnames(PRR7_d1d2_fishers)[5] <- "PRR7"
PRR7_d1d2_fishers <- PRR7_d1d2_fishers %>% 
  dplyr::select(PRR7) %>% 
  round(2)

LUX_d1d2_fishers <- get_odds(LUX_d1d2_fishers_prep)
colnames(LUX_d1d2_fishers)[5] <- "LUX"
LUX_d1d2_fishers <- LUX_d1d2_fishers %>% 
  dplyr::select(LUX) %>% 
  round(2)

ELF3_d1d2_fishers <- get_odds(ELF3_d1d2_fishers_prep)
colnames(ELF3_d1d2_fishers)[5] <- "ELF3"
ELF3_d1d2_fishers <- ELF3_d1d2_fishers %>% 
  dplyr::select(ELF3) %>% 
  round(2)

ELF4_d1d2_fishers <- get_odds(ELF4_d1d2_fishers_prep)
colnames(ELF4_d1d2_fishers)[5] <- "ELF4"
ELF4_d1d2_fishers <- ELF4_d1d2_fishers %>% 
  dplyr::select(ELF4) %>% 
  round(2)

odds_d1d2 <- bind_cols(clusters_d1d2,
                       LHY_d1d2_fishers, CCA1_d1d2_fishers, TOC1_d1d2_fishers, PRR5_d1d2_fishers,
                       PRR7_d1d2_fishers, LUX_d1d2_fishers, ELF3_d1d2_fishers, ELF4_d1d2_fishers) %>% 
  left_join(TF_network_clusters_count) %>% 
  dplyr::rename(cluster_size = cluster_number) %>% 
  relocate(cluster_size, .before = LHY)

#d1-d5----
LHY_d1d5_fishers <- get_odds(LHY_d1d5_fishers_prep)
colnames(LHY_d1d5_fishers)[5] <- "LHY"
LHY_d1d5_fishers <- LHY_d1d5_fishers %>% 
  dplyr::select(LHY) %>% 
  round(2)

CCA1_d1d5_fishers <- get_odds(CCA1_d1d5_fishers_prep)
colnames(CCA1_d1d5_fishers)[5] <- "CCA1"
CCA1_d1d5_fishers <- CCA1_d1d5_fishers %>% 
  dplyr::select(CCA1) %>% 
  round(2)

TOC1_d1d5_fishers <- get_odds(TOC1_d1d5_fishers_prep)
colnames(TOC1_d1d5_fishers)[5] <- "TOC1"
TOC1_d1d5_fishers <- TOC1_d1d5_fishers %>% 
  dplyr::select(TOC1) %>% 
  round(2)

PRR5_d1d5_fishers <- get_odds(PRR5_d1d5_fishers_prep)
colnames(PRR5_d1d5_fishers)[5] <- "PRR5"
PRR5_d1d5_fishers <- PRR5_d1d5_fishers %>% 
  dplyr::select(PRR5) %>% 
  round(2)

PRR7_d1d5_fishers <- get_odds(PRR7_d1d5_fishers_prep)
colnames(PRR7_d1d5_fishers)[5] <- "PRR7"
PRR7_d1d5_fishers <- PRR7_d1d5_fishers %>% 
  dplyr::select(PRR7) %>% 
  round(2)

LUX_d1d5_fishers <- get_odds(LUX_d1d5_fishers_prep)
colnames(LUX_d1d5_fishers)[5] <- "LUX"
LUX_d1d5_fishers <- LUX_d1d5_fishers %>% 
  dplyr::select(LUX) %>% 
  round(2)

ELF3_d1d5_fishers <- get_odds(ELF3_d1d5_fishers_prep)
colnames(ELF3_d1d5_fishers)[5] <- "ELF3"
ELF3_d1d5_fishers <- ELF3_d1d5_fishers %>% 
  dplyr::select(ELF3) %>% 
  round(2)

ELF4_d1d5_fishers <- get_odds(ELF4_d1d5_fishers_prep)
colnames(ELF4_d1d5_fishers)[5] <- "ELF4"
ELF4_d1d5_fishers <- ELF4_d1d5_fishers %>% 
  dplyr::select(ELF4) %>% 
  round(2)

odds_d1d5 <- bind_cols(clusters_d1d5,
                       LHY_d1d5_fishers, CCA1_d1d5_fishers, TOC1_d1d5_fishers, PRR5_d1d5_fishers,
                       PRR7_d1d5_fishers, LUX_d1d5_fishers, ELF3_d1d5_fishers, ELF4_d1d5_fishers) %>% 
  left_join(TF_network_clusters_count) %>% 
  dplyr::rename(cluster_size = cluster_number) %>% 
  relocate(cluster_size, .before = LHY)

#heatmaps----
#d1-d2----
#* full----
odds_d1d2_matrix <- odds_d1d2 %>%
  dplyr::select(4:11)
  
odds_d1d2_matrix<-as.matrix(odds_d1d2_matrix)

rownames(odds_d1d2_matrix) <- 1:74

col_fun = colorRamp2(c(0,15,30), c('white','blue','red'))

col_fun2 = colorRamp2(c(0,15,35), c('#FFFFFFFF','#9ecae1','#3182bd'))

col_fun(seq(-5,5))
col_fun2(seq(-5,5))

row_ha_hmap_d1d2_full = rowAnnotation(context = odds_d1d2$group, size = anno_barplot(odds_d1d2$cluster_size), col = list(context = c("gain high" = "#31a354", "gain medium" = "#74c476", "other" = "#cccccc", "lose medium" = "#fb6a4a", "lose high" = "#de2d26")), show_legend = FALSE)

m_d1d2 = Heatmap(odds_d1d2_matrix,
                 name='Odds Ratio', 
                 heatmap_legend_param = list(legend_direction = "horizontal", 
                                      title_gp = gpar(fontsize = 10), 
                                      legend_width = unit(5, "cm"), at = c(0, 5, 10, 15, 20, 25, 30), 
                                      labels = c(0, 5, 10, 15, 20, 25, 30), 
                                      title = "Odds Ratio", 
                                      legend_height = unit(4, "cm"), 
                                      title_position = "topleft", border="gray40"), 
                 column_names_rot = 45, 
                 column_title = NULL, 
                 rect_gp = gpar(col = "gray40", lwd = 1), 
                 col=col_fun2,
                 right_annotation = row_ha_hmap_d1d2_full,
                 row_title=NULL, 
                 row_dend_width = unit(4, "cm"), 
                 row_names_gp = gpar(fontsize = 8), 
                 row_gap = unit(2, "mm"), 
                 column_km=6, 
                 column_km_repeats = 100,
                 column_dend_height = unit(2, "cm"),
                 column_gap = unit(5, "mm"), 
                 left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c(8,8,8,8)), 
                                                           labels = c("g1", "g2", "g3", "g4"),
                                                           labels_gp = gpar(col = "black", fontsize = 10))), 
                 row_km=4, 
                 row_km_repeats = 100, 
                 cell_fun = function(j, i, x, y, width, height, fill) {if(odds_d1d2_matrix[i, j] > 5) grid.text(sprintf("%.1f", odds_d1d2_matrix[i, j]), x, y, gp = gpar(fontsize =8, col='black'))})


lgd_list=list(Legend(labels=c("Gain High", "Gain Medium", "Other", "Lose Medium", "Lose High"), title="cluster context", type = "points", pch=15, title_gp = gpar(fontsize = 10), size=unit(5,"mm"), border="black", legend_gp=gpar(col=c('#31a354', '#74c476', '#cccccc', '#fb6a4a', '#de2d26'))))

draw(m_d1d2, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list)

decorate_annotation("size", {grid.text("cluster size", y = unit(1, "npc") + unit(2, "mm"), just = "bottom", gp = gpar(fontsize = 10))})

# https://jokergoo.github.io/2020/05/11/set-cell-width/height-in-the-heatmap/

m = draw(m_d1d2, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list)
w = ComplexHeatmap:::width(m)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(m)
h = convertY(h, "inch", valueOnly = TRUE)
c(w, h)

# dev.off()
# plot(rnorm(50), rnorm(50))
# dev.cur()
# dev.off(2)

#* trimmed----
odds_d1d2_trimmed <- odds_d1d2 %>%
  filter_at(vars(LHY:ELF4), any_vars(.>5)) 
  
odds_d1d2_trimmed_matrix <- odds_d1d2_trimmed %>%
  dplyr::select(4:11)
  
odds_d1d2_trimmed_matrix <- as.matrix(odds_d1d2_trimmed_matrix)

rownames(odds_d1d2_trimmed_matrix) <- c(9, 10, 11, 20, 22, 24, 25, 29, 31, 34, 38, 42, 44, 52, 56, 58, 65, 67, 71, 72)

row_ha_hmap_d1d2_trimmed = rowAnnotation(context = odds_d1d2_trimmed$group, size = anno_barplot(odds_d1d2_trimmed$cluster_size), col = list(context = c("gain high" = "#31a354", "gain medium" = "#74c476", "other" = "#cccccc", "lose medium" = "#fb6a4a", "lose high" = "#de2d26")), show_legend = FALSE)
draw(row_ha_hmap_d1d2_trimmed)
m_d1d2_trimmed = Heatmap(odds_d1d2_trimmed_matrix,
                         name='Odds Ratio',
                         heatmap_legend_param = list(legend_direction = "horizontal", 
                                                     title_gp = gpar(fontsize = 10), 
                                                     legend_width = unit(5, "cm"), at = c(0, 5, 10, 15, 20, 25, 30, 35), 
                                                     labels = c(0, 5, 10, 15, 20, 25, 30, 35), 
                                                     title = "Odds Ratio",
                                                     legend_height = unit(4, "cm"), 
                                                     title_position = "topleft", border="gray40"), 
                         column_names_rot = 45, 
                         column_title = NULL, 
                         rect_gp = gpar(col = "gray40", lwd = 1), 
                         col=col_fun2,
                         right_annotation = row_ha_hmap_d1d2_trimmed,
                         row_title=NULL, 
                         row_dend_width = unit(4, "cm"),
                         row_names_gp = gpar(fontsize = 10), 
                         row_gap = unit(2, "mm"), 
                         column_km=6, 
                         column_km_repeats = 100,
                         column_dend_height = unit(2, "cm"),
                         column_gap = unit(5, "mm"), 
                         left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c(8,8,8,8)), 
                                                                  labels = c("g1", "g2", "g3", "g4"),
                                                                  labels_gp = gpar(col = "black", fontsize = 10))), 
                         row_km=4, 
                         row_km_repeats = 100,
                         cell_fun = function(j, i, x, y, width, height, fill) {if(odds_d1d2_trimmed_matrix[i, j] > 5) grid.text(sprintf("%.1f", odds_d1d2_trimmed_matrix[i, j]), x, y, gp = gpar(fontsize =11, col='black', fontface = 'bold'))})


lgd_list=list(Legend(labels=c("Gain High", "Gain Medium", "Other", "Lose Medium", "Lose High"), title="cluster context", type = "points", pch=15, title_gp = gpar(fontsize = 10), size=unit(5,"mm"), border="black", legend_gp=gpar(col=c('#31a354', '#74c476', '#cccccc', '#fb6a4a', '#de2d26'))))

draw(m_d1d2_trimmed, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list)

decorate_annotation("size", {grid.text("cluster size", y = unit(1, "npc") + unit(2, "mm"), just = "bottom", gp = gpar(fontsize = 10))})

# https://jokergoo.github.io/2020/05/11/set-cell-width/height-in-the-heatmap/

m = draw(m_d1d2_trimmed, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list)
w = ComplexHeatmap:::width(m)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(m)
h = convertY(h, "inch", valueOnly = TRUE)
c(w, h)

#d1-d5----

odds_d1d5_matrix <- odds_d1d5 %>%
  dplyr::select(4:11)

odds_d1d5_matrix<-as.matrix(odds_d1d5_matrix)

rownames(odds_d1d5_matrix) <- 1:74

row_ha_hmap_d1d5_full = rowAnnotation(context = odds_d1d5$group, size = anno_barplot(odds_d1d5$cluster_size), col = list(context = c("gain high" = "#31a354", "gain medium" = "#74c476", "other" = "#cccccc", "lose medium" = "#fb6a4a", "lose high" = "#de2d26")), show_legend = FALSE)

m_d1d5 = Heatmap(odds_d1d5_matrix,
                 name='Odds Ratio', 
                 heatmap_legend_param = list(legend_direction = "horizontal", 
                                             title_gp = gpar(fontsize = 10), 
                                             legend_width = unit(5, "cm"), at = c(0, 5, 10, 15, 20, 25, 30), 
                                             labels = c(0, 5, 10, 15, 20, 25, 30), 
                                             title = "Odds Ratio", 
                                             legend_height = unit(4, "cm"), 
                                             title_position = "topleft", border="gray40"), 
                 column_names_rot = 45, 
                 column_title = NULL, 
                 rect_gp = gpar(col = "gray40", lwd = 1), 
                 col=col_fun2,
                 right_annotation = row_ha_hmap_d1d5_full,
                 row_title=NULL, 
                 row_dend_width = unit(4, "cm"), 
                 row_names_gp = gpar(fontsize = 8), 
                 row_gap = unit(2, "mm"), 
                 column_km=6, 
                 column_km_repeats = 100,
                 column_dend_height = unit(2, "cm"),
                 column_gap = unit(5, "mm"), 
                 left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c(8,8,8,8)), 
                                                                  labels = c("g1", "g2", "g3", "g4"),
                                                                  labels_gp = gpar(col = "black", fontsize = 10))), 
                 row_km=4, 
                 row_km_repeats = 100, 
                 cell_fun = function(j, i, x, y, width, height, fill) {if(odds_d1d5_matrix[i, j] > 5) grid.text(sprintf("%.1f", odds_d1d5_matrix[i, j]), x, y, gp = gpar(fontsize =8, col='black'))})

draw(m_d1d5, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list)

decorate_annotation("size", {grid.text("cluster size", y = unit(1, "npc") + unit(2, "mm"), just = "bottom", gp = gpar(fontsize = 10))})

# https://jokergoo.github.io/2020/05/11/set-cell-width/height-in-the-heatmap/

m = draw(m_d1d5, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list)
w = ComplexHeatmap:::width(m)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(m)
h = convertY(h, "inch", valueOnly = TRUE)
c(w, h)

#* trimmed----
odds_d1d5_trimmed <- odds_d1d5 %>%
  filter_at(vars(LHY:ELF4), any_vars(.>5)) 

odds_d1d5_trimmed_matrix <- odds_d1d5_trimmed %>%
  dplyr::select(4:11)

odds_d1d5_trimmed_matrix <- as.matrix(odds_d1d5_trimmed_matrix)

rownames(odds_d1d5_trimmed_matrix) <- c(9, 10, 11, 20, 22, 24, 25, 29, 31, 34, 38, 42, 44, 52, 56, 58, 65, 67, 71, 72)

row_ha_hmap_d1d5_trimmed = rowAnnotation(context = odds_d1d5_trimmed$group, size = anno_barplot(odds_d1d5_trimmed$cluster_size), col = list(context = c("gain high" = "#31a354", "gain medium" = "#74c476", "other" = "#cccccc", "lose medium" = "#fb6a4a", "lose high" = "#de2d26")), show_legend = FALSE)

m_d1d5_trimmed = Heatmap(odds_d1d5_trimmed_matrix,
                         name='Odds Ratio',
                         heatmap_legend_param = list(legend_direction = "horizontal", 
                                                     title_gp = gpar(fontsize = 10), 
                                                     legend_width = unit(5, "cm"), at = c(0, 5, 10, 15, 20, 25, 30, 35), 
                                                     labels = c(0, 5, 10, 15, 20, 25, 30, 35), 
                                                     title = "Odds Ratio",
                                                     legend_height = unit(4, "cm"), 
                                                     title_position = "topleft", border="gray40"), 
                         column_names_rot = 45, 
                         column_title = NULL, 
                         rect_gp = gpar(col = "gray40", lwd = 1), 
                         col=col_fun2,
                         right_annotation = row_ha_hmap_d1d5_trimmed,
                         row_title=NULL, 
                         row_dend_width = unit(4, "cm"),
                         row_names_gp = gpar(fontsize = 10), 
                         row_gap = unit(2, "mm"), 
                         column_km=6, 
                         column_km_repeats = 100,
                         column_dend_height = unit(2, "cm"),
                         column_gap = unit(5, "mm"), 
                         left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c(8,8,8,8)), 
                                                                          labels = c("g1", "g2", "g3", "g4"),
                                                                          labels_gp = gpar(col = "black", fontsize = 10))), 
                         row_km=4, 
                         row_km_repeats = 100,
                         cell_fun = function(j, i, x, y, width, height, fill) {if(odds_d1d5_trimmed_matrix[i, j] > 5) grid.text(sprintf("%.1f", odds_d1d5_trimmed_matrix[i, j]), x, y, gp = gpar(fontsize =11, col='black', fontface = 'bold'))})


lgd_list=list(Legend(labels=c("Gain High", "Gain Medium", "Other", "Lose Medium", "Lose High"), title="cluster context", type = "points", pch=15, title_gp = gpar(fontsize = 10), size=unit(5,"mm"), border="black", legend_gp=gpar(col=c('#31a354', '#74c476', '#cccccc', '#fb6a4a', '#de2d26'))))

draw(m_d1d5_trimmed, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list)

decorate_annotation("size", {grid.text("cluster size", y = unit(1, "npc") + unit(2, "mm"), just = "bottom", gp = gpar(fontsize = 10))})

# https://jokergoo.github.io/2020/05/11/set-cell-width/height-in-the-heatmap/

m = draw(m_d1d5_trimmed, heatmap_legend_side = "bottom", annotation_legend_list = lgd_list)
w = ComplexHeatmap:::width(m)
w = convertX(w, "inch", valueOnly = TRUE)
h = ComplexHeatmap:::height(m)
h = convertY(h, "inch", valueOnly = TRUE)
c(w, h)

# UpSet Plots----

# full----

myGeneSets <- list(TF_network_LHY = TF_adams_merge$gene_ID,
                   TF_network_CCA1 = TF_kamioka_nagel_merge$gene_ID,
                   TF_network_TOC1 = TF_huang_merge$gene_ID,
                   TF_network_PRR5 = TF_nakamichi_merge$gene_ID,
                   TF_network_PRR7 = TF_liu_merge$gene_ID,
                   TF_network_LUX = TF_ezer_LUX_merge$gene_ID,
                   TF_network_ELF3 = TF_ezer_ELF3_merge$gene_ID,
                   TF_network_ELF4 = TF_ezer_ELF4_merge$gene_ID)

# fromList: a function to convert a list of named vectors to a data frame compatible with UpSetR
sets <- fromList(myGeneSets)

UpSet<-upset(sets, nsets=8, number.angles = 30, order.by = "freq", matrix.color='grey40', point.size = 2.5, sets.x.label = "Clock ChIP targets", mainbar.y.label = "Gene Set Intersections", sets.bar.color = c("#542788", "#a6d96a", "#4575b4", "#1a9850", "#636363", "#f46d43", "#fdae61", "#cccccc"), text.scale = c(1.3, 1.3, 1, 1, 1, 0.75))

UpSet

png("Upset LHY CCA1 TOC1 LUX.png", width = 6, height = 6, units = 'in', res = 300)

# for all 8 clock
# png("Upset LHY CCA1 TOC1 LUX PRR5 PRR7 ELF3 ELF4.png", width = 8, height = 6, units = 'in', res = 300)

UpSet

dev.off()

# trimmed----

TF_adams_merge_trimmed <- TF_adams_merge %>% 
  filter(cluster %in% c(9, 10, 11, 20, 22, 24, 25, 29, 31, 34, 38, 42, 44, 52, 56, 58, 65, 67, 71, 72))

TF_kamioka_nagel_merge_trimmed <- TF_kamioka_nagel_merge %>%
  filter(cluster %in% c(9, 10, 11, 20, 22, 24, 25, 29, 31, 34, 38, 42, 44, 52, 56, 58, 65, 67, 71, 72))

TF_huang_merge_trimmed <- TF_huang_merge %>%
  filter(cluster %in% c(9, 10, 11, 20, 22, 24, 25, 29, 31, 34, 38, 42, 44, 52, 56, 58, 65, 67, 71, 72))

TF_nakamichi_merge_trimmed <- TF_nakamichi_merge %>%
  filter(cluster %in% c(9, 10, 11, 20, 22, 24, 25, 29, 31, 34, 38, 42, 44, 52, 56, 58, 65, 67, 71, 72))

TF_liu_merge_trimmed <- TF_liu_merge %>%
  filter(cluster %in% c(9, 10, 11, 20, 22, 24, 25, 29, 31, 34, 38, 42, 44, 52, 56, 58, 65, 67, 71, 72))

TF_ezer_LUX_merge_trimmed <- TF_ezer_LUX_merge %>%
  filter(cluster %in% c(9, 10, 11, 20, 22, 24, 25, 29, 31, 34, 38, 42, 44, 52, 56, 58, 65, 67, 71, 72))

TF_ezer_ELF3_merge_trimmed <- TF_ezer_ELF3_merge %>%
  filter(cluster %in% c(9, 10, 11, 20, 22, 24, 25, 29, 31, 34, 38, 42, 44, 52, 56, 58, 65, 67, 71, 72))

TF_ezer_ELF4_merge_trimmed <- TF_ezer_ELF4_merge %>%
  filter(cluster %in% c(9, 10, 11, 20, 22, 24, 25, 29, 31, 34, 38, 42, 44, 52, 56, 58, 65, 67, 71, 72))
  
myGeneSets_trimmed <- list(LHY = TF_adams_merge_trimmed$gene_ID,
                           CCA1 = TF_kamioka_nagel_merge_trimmed$gene_ID,
                           TOC1 = TF_huang_merge_trimmed$gene_ID,
                           PRR5 = TF_nakamichi_merge_trimmed$gene_ID,
                           PRR7 = TF_liu_merge_trimmed$gene_ID,
                           LUX = TF_ezer_LUX_merge_trimmed$gene_ID,
                           ELF3 = TF_ezer_ELF3_merge_trimmed$gene_ID,
                           ELF4 = TF_ezer_ELF4_merge_trimmed$gene_ID) 
  
sets_trimmed <- fromList(myGeneSets_trimmed)

UpSet_trimmed<-upset(sets_trimmed, nsets=8, number.angles = 30, order.by = "freq", matrix.color='grey40', point.size = 2.5, sets.x.label = "Clock ChIP targets", mainbar.y.label = "Gene Set Intersections", sets.bar.color = c("#542788", "#a6d96a", "#4575b4", "#1a9850", "#636363", "#f46d43", "#fdae61", "#cccccc"), text.scale = c(1.3, 1.3, 1, 1, 1, 0.9))

png("./03_plots/UpSet_trimmed.png", width = 8, height = 6, units = 'in', res = 300)

UpSet_trimmed

dev.off()

# UpSet column identities----

# col1----
# LHY targets alone i.e. not targets of CCA1, TOC1, LUX, ELF3, PRR7, PRR5 and ELF4
# firstly an anti_join of CCA1 with LHY
TF_LHY_TOC1 <- anti_join(TF_adams_merge_trimmed, TF_huang_merge_trimmed, by="gene_ID")

# take this and rule out LUX targets
TF_LHY_TOC1_notLUX <- anti_join(TF_LHY_TOC1, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_LHY_TOC1_notLUX_notCCA1 <- anti_join(TF_LHY_TOC1_notLUX, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_LHY_TOC1_notLUX_notCCA1_notELF3 <- anti_join(TF_LHY_TOC1_notLUX_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR5 <- anti_join(TF_LHY_TOC1_notLUX_notCCA1_notELF3, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR5_notPRR7 <- anti_join(TF_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR5, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR5_notPRR7_notELF4 <- anti_join(TF_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR5_notPRR7, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 1 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col1 <- TF_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR5_notPRR7_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col1$clock <- "LHY"

UpSet_trimmed_col1 <- UpSet_trimmed_col1 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col1 <- UpSet_trimmed_col1 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col2----
# TOC1 targets alone i.e. not targets of CCA1, LHY, LUX, ELF3, PRR7, PRR5 and ELF4
# firstly an anti_join of CCA1 with LHY
TF_TOC1_LHY <- anti_join(TF_huang_merge_trimmed, TF_adams_merge_trimmed, by="gene_ID")

# take this and rule out LUX targets
TF_TOC1_LHY_notLUX <- anti_join(TF_TOC1_LHY, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_TOC1_LHY_notLUX_notCCA1 <- anti_join(TF_TOC1_LHY_notLUX, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_TOC1_LHY_notLUX_notCCA1_notELF3 <- anti_join(TF_TOC1_LHY_notLUX_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_TOC1_LHY_notLUX_notCCA1_notELF3_notPRR5 <- anti_join(TF_TOC1_LHY_notLUX_notCCA1_notELF3, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_TOC1_LHY_notLUX_notCCA1_notELF3_notPRR5_notPRR7 <- anti_join(TF_TOC1_LHY_notLUX_notCCA1_notELF3_notPRR5, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_TOC1_LHY_notLUX_notCCA1_notELF3_notPRR5_notPRR7_notELF4 <- anti_join(TF_TOC1_LHY_notLUX_notCCA1_notELF3_notPRR5_notPRR7, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 2 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col2 <- TF_TOC1_LHY_notLUX_notCCA1_notELF3_notPRR5_notPRR7_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col2$clock <- "TOC1"

UpSet_trimmed_col2 <- UpSet_trimmed_col2 %>% 
  left_join(TOC1_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col2 <- UpSet_trimmed_col2 %>% 
  left_join(TOC1_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col3----
# LUX targets alone i.e. not targets of CCA1, LHY, TOC1, ELF3, PRR7, PRR5 and ELF4
# firstly an anti_join of CCA1 with LHY
TF_LUX_LHY <- anti_join(TF_ezer_LUX_merge_trimmed, TF_adams_merge_trimmed, by="gene_ID")

# take this and rule out TOC1 targets
TF_LUX_LHY_notTOC1 <- anti_join(TF_LUX_LHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_LUX_LHY_notTOC1_notCCA1 <- anti_join(TF_LUX_LHY_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_LUX_LHY_notTOC1_notCCA1_notELF3 <- anti_join(TF_LUX_LHY_notTOC1_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR5 <- anti_join(TF_LUX_LHY_notTOC1_notCCA1_notELF3, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR5_notPRR7 <- anti_join(TF_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR5, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR5_notPRR7_notELF4 <- anti_join(TF_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR5_notPRR7, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 3 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col3 <- TF_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR5_notPRR7_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col3$clock <- "LUX"

UpSet_trimmed_col3 <- UpSet_trimmed_col3 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col3 <- UpSet_trimmed_col3 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col4----
# CCA1 and LHY common targets - not targets of LUX, TOC1, ELF3, PRR7, PRR5 and ELF4
# firstly an inner_join of CCA1 with LHY
TF_common_LHY_CCA1 <- inner_join(TF_adams_merge_trimmed[,1:2], TF_kamioka_nagel_merge_trimmed[,2], by="gene_ID")

# take this and rule out TOC1 targets
TF_common_LHY_CCA1_notTOC1 <- anti_join(TF_common_LHY_CCA1, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out LUX targets
TF_common_LHY_CCA1_notTOC1_notLUX <- anti_join(TF_common_LHY_CCA1_notTOC1, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LHY_CCA1_notTOC1_notLUX_notELF3 <- anti_join(TF_common_LHY_CCA1_notTOC1_notLUX, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LHY_CCA1_notTOC1_notLUX_notELF3_notPRR5 <- anti_join(TF_common_LHY_CCA1_notTOC1_notLUX_notELF3, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LHY_CCA1_notTOC1_notLUX_notELF3_notPRR5_notPRR7 <- anti_join(TF_common_LHY_CCA1_notTOC1_notLUX_notELF3_notPRR5, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LHY_CCA1_notTOC1_notLUX_notELF3_notPRR5_notPRR7_notELF4 <- anti_join(TF_common_LHY_CCA1_notTOC1_notLUX_notELF3_notPRR5_notPRR7, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 4 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col4 <- TF_common_LHY_CCA1_notTOC1_notLUX_notELF3_notPRR5_notPRR7_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col4$clock <- "LHY and CCA1"

UpSet_trimmed_col4 <- UpSet_trimmed_col4 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col4 <- UpSet_trimmed_col4 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col5----
# LUX and ELF3 common targets - not targets of LHY, CCA1, TOC1, PRR7, PRR5 and ELF4
# firstly an inner_join of LUX with ELF3
TF_common_LUX_ELF3 <- inner_join(TF_ezer_LUX_merge_trimmed[,1:2], TF_ezer_ELF3_merge_trimmed[,2], by="gene_ID")

# take this and rule out LHY targets
TF_common_LUX_ELF3_notLHY <- anti_join(TF_common_LUX_ELF3, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_ELF3_notLHY_notCCA1 <- anti_join(TF_common_LUX_ELF3_notLHY, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_common_LUX_ELF3_notLHY_notCCA1_notTOC1 <- anti_join(TF_common_LUX_ELF3_notLHY_notCCA1, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_ELF3_notLHY_notCCA1_notTOC1_notPRR5 <- anti_join(TF_common_LUX_ELF3_notLHY_notCCA1_notTOC1, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_ELF3_notLHY_notCCA1_notTOC1_notPRR5_notPRR7 <- anti_join(TF_common_LUX_ELF3_notLHY_notCCA1_notTOC1_notPRR5, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_ELF3_notLHY_notCCA1_notTOC1_notPRR5_notPRR7_notELF4 <- anti_join(TF_common_LUX_ELF3_notLHY_notCCA1_notTOC1_notPRR5_notPRR7, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 5 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col5 <- TF_common_LUX_ELF3_notLHY_notCCA1_notTOC1_notPRR5_notPRR7_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col5$clock <- "LUX and ELF3"

UpSet_trimmed_col5 <- UpSet_trimmed_col5 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col5 <- UpSet_trimmed_col5 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col6----
# PRR5 targets alone - not targets of LUX, LHY, TOC1, CCA1, ELF3, PRR7 and ELF4
# firstly an anti_join of PRR5 with LUX
TF_PRR5_LUX <- anti_join(TF_nakamichi_merge_trimmed, TF_ezer_LUX_merge_trimmed, by="gene_ID")

# take this and rule out LHY targets
TF_PRR5_LUX_notLHY <- anti_join(TF_PRR5_LUX, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_PRR5_LUX_notLHY_notTOC1 <- anti_join(TF_PRR5_LUX_notLHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_PRR5_LUX_notLHY_notTOC1_notCCA1 <- anti_join(TF_PRR5_LUX_notLHY_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_PRR5_LUX_notLHY_notTOC1_notCCA1_notELF3 <- anti_join(TF_PRR5_LUX_notLHY_notTOC1_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_PRR5_LUX_notLHY_notTOC1_notCCA1_notELF3_notPRR7 <- anti_join(TF_PRR5_LUX_notLHY_notTOC1_notCCA1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_PRR5_LUX_notLHY_notTOC1_notCCA1_notELF3_notPRR7_notELF4 <- anti_join(TF_PRR5_LUX_notLHY_notTOC1_notCCA1_notELF3_notPRR7, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 6 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col6 <- TF_PRR5_LUX_notLHY_notTOC1_notCCA1_notELF3_notPRR7_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col6$clock <- "PRR5"

UpSet_trimmed_col6 <- UpSet_trimmed_col6 %>% 
  left_join(PRR5_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col6 <- UpSet_trimmed_col6 %>% 
  left_join(PRR5_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col7----
# PRR7 targets alone - not targets of LUX, LHY, TOC1, CCA1, ELF3, PRR5 and ELF4
# firstly an anti_join of PRR7 with LUX
TF_PRR7_LUX <- anti_join(TF_liu_merge_trimmed, TF_ezer_LUX_merge_trimmed, by="gene_ID")

# take this and rule out LHY targets
TF_PRR7_LUX_notLHY <- anti_join(TF_PRR7_LUX, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_PRR7_LUX_notLHY_notTOC1 <- anti_join(TF_PRR7_LUX_notLHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_PRR7_LUX_notLHY_notTOC1_notCCA1 <- anti_join(TF_PRR7_LUX_notLHY_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_PRR7_LUX_notLHY_notTOC1_notCCA1_notELF3 <- anti_join(TF_PRR7_LUX_notLHY_notTOC1_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_PRR7_LUX_notLHY_notTOC1_notCCA1_notELF3_notPRR5 <- anti_join(TF_PRR7_LUX_notLHY_notTOC1_notCCA1_notELF3, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_PRR7_LUX_notLHY_notTOC1_notCCA1_notELF3_notPRR5_notELF4 <- anti_join(TF_PRR7_LUX_notLHY_notTOC1_notCCA1_notELF3_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 7 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col7 <- TF_PRR7_LUX_notLHY_notTOC1_notCCA1_notELF3_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col7$clock <- "PRR7"

UpSet_trimmed_col7 <- UpSet_trimmed_col7 %>% 
  left_join(PRR7_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col7 <- UpSet_trimmed_col7 %>% 
  left_join(PRR7_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col8----
# LUX and LHY common targets - not targets of TOC1, CCA1, ELF3, PRR7, PRR5 and ELF4
# firstly an inner_join of PRR7 with LUX
TF_common_LUX_LHY <- inner_join(TF_ezer_LUX_merge_trimmed[,1:2], TF_adams_merge_trimmed[,2], by="gene_ID")

# take this and rule out TOC1 targets
TF_common_LUX_LHY_notTOC1 <- anti_join(TF_common_LUX_LHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_LHY_notTOC1_notCCA1 <- anti_join(TF_common_LUX_LHY_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LUX_LHY_notTOC1_notCCA1_notELF3 <- anti_join(TF_common_LUX_LHY_notTOC1_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR7 <- anti_join(TF_common_LUX_LHY_notTOC1_notCCA1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR7_notPRR5 <- anti_join(TF_common_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 8 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col8 <- TF_common_LUX_LHY_notTOC1_notCCA1_notELF3_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col8$clock <- "LUX and LHY"

UpSet_trimmed_col8 <- UpSet_trimmed_col8 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col8 <- UpSet_trimmed_col8 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col9----
# LUX and LHY and CCA1 common targets - not targets of TOC1, ELF3, PRR7, PRR5 and ELF4
# firstly take TF_common_LUX_LHY and then inner_join with CCA1
TF_common_LUX_LHY_CCA1 <- inner_join(TF_common_LUX_LHY[,1:2], TF_kamioka_nagel_merge_trimmed[,2], by="gene_ID")

# take this and rule out TOC1 targets
TF_common_LUX_LHY_CCA1_notTOC1 <- anti_join(TF_common_LUX_LHY_CCA1, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LUX_LHY_CCA1_notTOC1_notELF3 <- anti_join(TF_common_LUX_LHY_CCA1_notTOC1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_LHY_CCA1_notTOC1_notELF3_notPRR7 <- anti_join(TF_common_LUX_LHY_CCA1_notTOC1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_LHY_CCA1_notTOC1_notELF3_notPRR7_notPRR5 <- anti_join(TF_common_LUX_LHY_CCA1_notTOC1_notELF3_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_LHY_CCA1_notTOC1_notELF3_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LUX_LHY_CCA1_notTOC1_notELF3_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 9 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col9 <- TF_common_LUX_LHY_CCA1_notTOC1_notELF3_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col9$clock <- "LUX and LHY and CCA1"

UpSet_trimmed_col9 <- UpSet_trimmed_col9 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col9 <- UpSet_trimmed_col9 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col10----
# LUX and TOC1 common targets - not targets of LHY, CCA1, ELF3, PRR7, PRR5 and ELF4
# firstly take an inner_join of LUX with TOC1
TF_common_LUX_TOC1 <- inner_join(TF_ezer_LUX_merge_trimmed[,1:2], TF_huang_merge_trimmed[,2], by="gene_ID")

# take this and rule out LHY targets
TF_common_LUX_TOC1_notLHY <- anti_join(TF_common_LUX_TOC1, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_TOC1_notLHY_notCCA1 <- anti_join(TF_common_LUX_TOC1_notLHY, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LUX_TOC1_notLHY_notCCA1_notELF3 <- anti_join(TF_common_LUX_TOC1_notLHY_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_TOC1_notLHY_notCCA1_notELF3_notPRR7 <- anti_join(TF_common_LUX_TOC1_notLHY_notCCA1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_TOC1_notLHY_notCCA1_notELF3_notPRR7_notPRR5 <- anti_join(TF_common_LUX_TOC1_notLHY_notCCA1_notELF3_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_TOC1_notLHY_notCCA1_notELF3_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LUX_TOC1_notLHY_notCCA1_notELF3_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 10 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col10 <- TF_common_LUX_TOC1_notLHY_notCCA1_notELF3_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col10$clock <- "LUX and TOC1"

UpSet_trimmed_col10 <- UpSet_trimmed_col10 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col10 <- UpSet_trimmed_col10 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)



