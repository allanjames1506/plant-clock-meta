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
  
sets_trimmed <- fromList(myGeneSets_trimmed) %>% 
  write_csv('./00_raw_data/sets_trimmed.csv')

UpSet_trimmed<-upset(sets_trimmed, 
                     nsets=8, 
                     number.angles = 30, 
                     order.by = "freq", 
                     matrix.color='grey40', 
                     point.size = 2.5, 
                     sets.x.label = "Clock ChIP targets", 
                     mainbar.y.label = "Gene Set Intersections", 
                     sets.bar.color = c("#542788", "#a6d96a", "#4575b4", "#1a9850", "#636363", "#f46d43", "#fdae61", "#cccccc"), 
                     text.scale = c(1.3, 1.3, 1, 1, 1, 0.9),
                     base_annotations=list(
                       'Intersection size'=intersection_size(
                         counts=FALSE,
                         mapping=aes(fill=mpaa)
                       )))

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

# col11----
# CCA1 targets alone - not targets of LUX, LHY, TOC1, ELF3, PRR7, PRR5 and ELF4
# firstly take an anti_join of CCA1 with LUX
TF_CCA1_LUX <- anti_join(TF_kamioka_nagel_merge_trimmed, TF_ezer_LUX_merge_trimmed, by="gene_ID")

# take this and rule out LHY targets
TF_CCA1_LUX_notLHY <- anti_join(TF_CCA1_LUX, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_CCA1_LUX_notLHY_notTOC1 <- anti_join(TF_CCA1_LUX_notLHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_CCA1_LUX_notLHY_notTOC1_notELF3 <- anti_join(TF_CCA1_LUX_notLHY_notTOC1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_CCA1_LUX_notLHY_notTOC1_notELF3_notPRR7 <- anti_join(TF_CCA1_LUX_notLHY_notTOC1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_CCA1_LUX_notLHY_notTOC1_notELF3_notPRR7_notPRR5 <- anti_join(TF_CCA1_LUX_notLHY_notTOC1_notELF3_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_CCA1_LUX_notLHY_notTOC1_notELF3_notPRR7_notPRR5_notELF4 <- anti_join(TF_CCA1_LUX_notLHY_notTOC1_notELF3_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 11 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col11 <- TF_CCA1_LUX_notLHY_notTOC1_notELF3_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col11$clock <- "CCA1"

UpSet_trimmed_col11 <- UpSet_trimmed_col11 %>% 
  left_join(CCA1_nagel_kamioka_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col11 <- UpSet_trimmed_col11 %>% 
  left_join(CCA1_nagel_kamioka_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col12----
# LUX and ELF3 and ELF4 common targets - not targets of LHY, TOC1, CCA1, PRR7, and PRR5 
# firstly take TF_common_LUX_ELF3 and do an inner_join with ELF4
TF_common_LUX_ELF3_ELF4 <- inner_join(TF_common_LUX_ELF3[,1:2], TF_ezer_ELF4_merge_trimmed[,2], by="gene_ID")

# take this and rule out LHY targets
TF_common_LUX_ELF3_ELF4_notLHY <- anti_join(TF_common_LUX_ELF3_ELF4, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_common_LUX_ELF3_ELF4_notLHY_notTOC1 <- anti_join(TF_common_LUX_ELF3_ELF4_notLHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_ELF3_ELF4_notLHY_notTOC1_notCCA1 <- anti_join(TF_common_LUX_ELF3_ELF4_notLHY_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_ELF3_ELF4_notLHY_notTOC1_notCCA1_notPRR7 <- anti_join(TF_common_LUX_ELF3_ELF4_notLHY_notTOC1_notCCA1, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_ELF3_ELF4_notLHY_notTOC1_notCCA1_notPRR7_notPRR5 <- anti_join(TF_common_LUX_ELF3_ELF4_notLHY_notTOC1_notCCA1_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 12 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col12 <- TF_common_LUX_ELF3_ELF4_notLHY_notTOC1_notCCA1_notPRR7_notPRR5[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col12$clock <- "LUX and ELF3 and ELF4"

UpSet_trimmed_col12 <- UpSet_trimmed_col12 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col12 <- UpSet_trimmed_col12 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col13----
# LHY and TOC1 common targets - not targets of LUX, CCA1, ELF3, PRR7, PRR5 and ELF4 
# firstly take LHY and do an inner_join with TOC1
TF_common_LHY_TOC1 <- inner_join(TF_adams_merge_trimmed[,1:2], TF_huang_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_common_LHY_TOC1_notLUX <- anti_join(TF_common_LHY_TOC1, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LHY_TOC1_notLUX_notCCA1 <- anti_join(TF_common_LHY_TOC1_notLUX, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LHY_TOC1_notLUX_notCCA1_notELF3 <- anti_join(TF_common_LHY_TOC1_notLUX_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR7 <- anti_join(TF_common_LHY_TOC1_notLUX_notCCA1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR7_notPRR5 <- anti_join(TF_common_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')


# tidy up dataset
# select first two columns
# Column 13 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col13 <- TF_common_LHY_TOC1_notLUX_notCCA1_notELF3_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col13$clock <- "LHY and TOC1"

UpSet_trimmed_col13 <- UpSet_trimmed_col13 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col13 <- UpSet_trimmed_col13 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col14----
# LHY and TOC1 and CCA1 common targets - not targets of LUX, ELF3, PRR7, PRR5 and ELF4 
# firstly take TF_common_LHY_TOC1 and do an inner_join with CCA1
TF_common_LHY_TOC1_CCA1 <- inner_join(TF_common_LHY_TOC1[,1:2], TF_kamioka_nagel_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_common_LHY_TOC1_CCA1_notLUX <- anti_join(TF_common_LHY_TOC1_CCA1, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LHY_TOC1_CCA1_notLUX_notELF3 <- anti_join(TF_common_LHY_TOC1_CCA1_notLUX, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LHY_TOC1_CCA1_notLUX_notELF3_notPRR7 <- anti_join(TF_common_LHY_TOC1_CCA1_notLUX_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LHY_TOC1_CCA1_notLUX_notELF3_notPRR7_notPRR5 <- anti_join(TF_common_LHY_TOC1_CCA1_notLUX_notELF3_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LHY_TOC1_CCA1_notLUX_notELF3_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LHY_TOC1_CCA1_notLUX_notELF3_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')


# tidy up dataset
# select first two columns
# Column 14 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col14 <- TF_common_LHY_TOC1_CCA1_notLUX_notELF3_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col14$clock <- "LHY and TOC1 and CCA1"

UpSet_trimmed_col14 <- UpSet_trimmed_col14 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col14 <- UpSet_trimmed_col14 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col15----
# LUX and TOC1 and PRR5 common targets - not targets of LHY, CCA1, ELF3, PRR7 and ELF4 
# firstly take TF_common_LUX_TOC1 and do an inner_join with PRR5
TF_common_LUX_TOC1_PRR5 <- inner_join(TF_common_LUX_TOC1[,1:2], TF_nakamichi_merge_trimmed[,2], by="gene_ID")

# take this and rule out LHY targets
TF_common_LUX_TOC1_PRR5_notLHY <- anti_join(TF_common_LUX_TOC1_PRR5, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_TOC1_PRR5_notLHY_notCCA1 <- anti_join(TF_common_LUX_TOC1_PRR5_notLHY, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LUX_TOC1_PRR5_notLHY_notCCA1_notELF3 <- anti_join(TF_common_LUX_TOC1_PRR5_notLHY_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_TOC1_PRR5_notLHY_notCCA1_notELF3_notPRR7 <- anti_join(TF_common_LUX_TOC1_PRR5_notLHY_notCCA1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_TOC1_PRR5_notLHY_notCCA1_notELF3_notPRR7_notELF4 <- anti_join(TF_common_LUX_TOC1_PRR5_notLHY_notCCA1_notELF3_notPRR7, TF_ezer_ELF4_merge_trimmed, by='gene_ID')


# tidy up dataset
# select first two columns
# Column 15 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col15 <- TF_common_LUX_TOC1_PRR5_notLHY_notCCA1_notELF3_notPRR7_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col15$clock <- "LUX and TOC1 and PRR5"

UpSet_trimmed_col15 <- UpSet_trimmed_col15 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col15 <- UpSet_trimmed_col15 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col16----
# LUX and LHY and ELF3 common targets - not targets of TOC1, CCA1, PRR7, PRR5 and ELF4 
# firstly take TF_common_LUX_LHY and do an inner_join with ELF3
TF_common_LUX_LHY_ELF3 <- inner_join(TF_common_LUX_LHY[,1:2], TF_ezer_ELF3_merge_trimmed[,2], by="gene_ID")

# take this and rule out TOC1 targets
TF_common_LUX_LHY_ELF3_notTOC1 <- anti_join(TF_common_LUX_LHY_ELF3, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_LHY_ELF3_notTOC1_notCCA1 <- anti_join(TF_common_LUX_LHY_ELF3_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_LHY_ELF3_notTOC1_notCCA1_notPRR7 <- anti_join(TF_common_LUX_LHY_ELF3_notTOC1_notCCA1, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_LHY_ELF3_notTOC1_notCCA1_notPRR7_notPRR5 <- anti_join(TF_common_LUX_LHY_ELF3_notTOC1_notCCA1_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_LHY_ELF3_notTOC1_notCCA1_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LUX_LHY_ELF3_notTOC1_notCCA1_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')


# tidy up dataset
# select first two columns
# Column 16 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col16 <- TF_common_LUX_LHY_ELF3_notTOC1_notCCA1_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col16$clock <- "LUX and LHY and ELF3"

UpSet_trimmed_col16 <- UpSet_trimmed_col16 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col16 <- UpSet_trimmed_col16 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col17----
# LUX and LHY and CCA1 and PRR7 common targets - not targets of TOC1, ELF3, PRR5 and ELF4 
# firstly take TF_common_LUX_LHY_CCA1 and do an inner_join with PRR7
TF_common_LUX_LHY_CCA1_PRR7 <- inner_join(TF_common_LUX_LHY_CCA1[,1:2], TF_liu_merge_trimmed[,2], by="gene_ID")

# take this and rule out TOC1 targets
TF_common_LUX_LHY_CCA1_PRR7_notTOC1 <- anti_join(TF_common_LUX_LHY_CCA1_PRR7, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LUX_LHY_CCA1_PRR7_notTOC1_notELF3 <- anti_join(TF_common_LUX_LHY_CCA1_PRR7_notTOC1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_LHY_CCA1_PRR7_notTOC1_notELF3_notPRR5 <- anti_join(TF_common_LUX_LHY_CCA1_PRR7_notTOC1_notELF3, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_LHY_CCA1_PRR7_notTOC1_notELF3_notPRR5_notELF4 <- anti_join(TF_common_LUX_LHY_CCA1_PRR7_notTOC1_notELF3_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')


# tidy up dataset
# select first two columns
# Column 17 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col17 <- TF_common_LUX_LHY_CCA1_PRR7_notTOC1_notELF3_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col17$clock <- "LUX and LHY and CCA1 and PRR7"

UpSet_trimmed_col17 <- UpSet_trimmed_col17 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col17 <- UpSet_trimmed_col17 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)


# col18----
# ELF3 targets alone - not targets of LUX, LHY, TOC1, CCA1, PRR7, PRR5 and ELF4 
# firstly an anti_join ELF3 with LUX
TF_ELF3_LUX <- anti_join(TF_ezer_ELF3_merge_trimmed, TF_ezer_LUX_merge_trimmed, by="gene_ID")

# take this and rule out LHY targets
TF_ELF3_LUX_notLHY <- anti_join(TF_ELF3_LUX, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_ELF3_LUX_notLHY_notTOC1 <- anti_join(TF_ELF3_LUX_notLHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_ELF3_LUX_notLHY_notTOC1_notCCA1 <- anti_join(TF_ELF3_LUX_notLHY_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_ELF3_LUX_notLHY_notTOC1_notCCA1_notPRR7 <- anti_join(TF_ELF3_LUX_notLHY_notTOC1_notCCA1, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_ELF3_LUX_notLHY_notTOC1_notCCA1_notPRR7_notPRR5 <- anti_join(TF_ELF3_LUX_notLHY_notTOC1_notCCA1_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_ELF3_LUX_notLHY_notTOC1_notCCA1_notPRR7_notPRR5_notELF4 <- anti_join(TF_ELF3_LUX_notLHY_notTOC1_notCCA1_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')


# tidy up dataset
# select first two columns
# Column 18 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col18 <- TF_ELF3_LUX_notLHY_notTOC1_notCCA1_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col18$clock <- "ELF3"

UpSet_trimmed_col18 <- UpSet_trimmed_col18 %>% 
  left_join(ELF3_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col18 <- UpSet_trimmed_col18 %>% 
  left_join(ELF3_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)


# col19----
# LUX and ELF3 and PRR5 and ELF4 common targets - not targets of LHY, TOC1, CCA1, and PRR7
# firstly start with TF_common_LUX_ELF3_ELF4 and do an inner_join with PRR5 

TF_common_LUX_ELF3_ELF4_PRR5 <- inner_join(TF_common_LUX_ELF3_ELF4[,1:2], TF_nakamichi_merge_trimmed[,2], by="gene_ID")

# take this and rule out LHY targets
TF_common_LUX_ELF3_ELF4_PRR5_notLHY <- anti_join(TF_common_LUX_ELF3_ELF4_PRR5, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_common_LUX_ELF3_ELF4_PRR5_notLHY_notTOC1 <- anti_join(TF_common_LUX_ELF3_ELF4_PRR5_notLHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_ELF3_ELF4_PRR5_notLHY_notTOC1_notCCA1 <- anti_join(TF_common_LUX_ELF3_ELF4_PRR5_notLHY_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_ELF3_ELF4_PRR5_notLHY_notTOC1_notCCA1_notPRR7 <- anti_join(TF_common_LUX_ELF3_ELF4_PRR5_notLHY_notTOC1_notCCA1, TF_liu_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 19 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col19 <- TF_common_LUX_ELF3_ELF4_PRR5_notLHY_notTOC1_notCCA1_notPRR7[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col19$clock <- "LUX and ELF3 and PRR5 and ELF4"

UpSet_trimmed_col19 <- UpSet_trimmed_col19 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col19 <- UpSet_trimmed_col19 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col20----
# LUX and TOC1 and ELF3 and PRR7 and PRR5 and ELF4 common targets - not targets of LHY and CCA1
# firstly start with TF_common_LUX_ELF3_ELF4_PRR5 and do an inner_join with TOC1 

TF_common_LUX_ELF3_ELF4_PRR5_TOC1 <- inner_join(TF_common_LUX_ELF3_ELF4_PRR5[,1:2], TF_huang_merge_trimmed[,2], by="gene_ID")

# take this and include PRR7 targets
TF_common_LUX_ELF3_ELF4_PRR5_TOC1_PRR7 <- inner_join(TF_common_LUX_ELF3_ELF4_PRR5_TOC1[,1:2], TF_liu_merge_trimmed[,2], by='gene_ID')

# take previous and rule out LHY targets
TF_common_LUX_ELF3_ELF4_PRR5_TOC1_PRR7_notLHY <- anti_join(TF_common_LUX_ELF3_ELF4_PRR5_TOC1_PRR7, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_ELF3_ELF4_PRR5_TOC1_PRR7_notLHY_notCCA1 <- anti_join(TF_common_LUX_ELF3_ELF4_PRR5_TOC1_PRR7_notLHY, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 20 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col20 <- TF_common_LUX_ELF3_ELF4_PRR5_TOC1_PRR7_notLHY_notCCA1[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col20$clock <- "LUX and TOC1 and ELF3 and PRR7 and PRR5 and ELF4"

UpSet_trimmed_col20 <- UpSet_trimmed_col20 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col20 <- UpSet_trimmed_col20 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col21----
# LUX and LHY and ELF3 and ELF4 common targets - not targets of TOC1, CCA1, PRR7 and PRR5
# firstly start with TF_common_LUX_LHY_ELF3 and do an inner_join with ELF4 

TF_common_LUX_LHY_ELF3_ELF4 <- inner_join(TF_common_LUX_LHY_ELF3[,1:2], TF_ezer_ELF4_merge_trimmed[,2], by="gene_ID")

# take this and rule out TOC1 targets
TF_common_LUX_LHY_ELF3_ELF4_notTOC1 <- anti_join(TF_common_LUX_LHY_ELF3_ELF4, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_LHY_ELF3_ELF4_notTOC1_notCCA1 <- anti_join(TF_common_LUX_LHY_ELF3_ELF4_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_LHY_ELF3_ELF4_notTOC1_notCCA1_notPRR7 <- anti_join(TF_common_LUX_LHY_ELF3_ELF4_notTOC1_notCCA1, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_LHY_ELF3_ELF4_notTOC1_notCCA1_notPRR7_notPRR5 <- anti_join(TF_common_LUX_LHY_ELF3_ELF4_notTOC1_notCCA1_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 21 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col21 <- TF_common_LUX_LHY_ELF3_ELF4_notTOC1_notCCA1_notPRR7_notPRR5[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col21$clock <- "LUX and LHY and ELF3 and ELF4"

UpSet_trimmed_col21 <- UpSet_trimmed_col21 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col21 <- UpSet_trimmed_col21 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col22----
# LUX and LHY and CCA1 and ELF3 and PRR7 and PRR5 targets - not targets of TOC1 and ELF4
# firstly start with TF_common_LUX_LHY_CCA1_PRR7 and do an inner_join with ELF3 and PRR5 

TF_common_LUX_LHY_CCA1_PRR7_ELF3 <- inner_join(TF_common_LUX_LHY_CCA1_PRR7[,1:2], TF_ezer_ELF3_merge_trimmed[,2], by="gene_ID")

# take this and include PRR5 targets
TF_common_LUX_LHY_CCA1_PRR7_ELF3_PRR5 <- inner_join(TF_common_LUX_LHY_CCA1_PRR7_ELF3[,1:2], TF_nakamichi_merge_trimmed[,2], by='gene_ID')

# take previous and rule out TOC1 targets
TF_common_LUX_LHY_CCA1_PRR7_ELF3_PRR5_notTOC1 <- anti_join(TF_common_LUX_LHY_CCA1_PRR7_ELF3_PRR5, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_LHY_CCA1_PRR7_ELF3_PRR5_notTOC1_notELF4 <- anti_join(TF_common_LUX_LHY_CCA1_PRR7_ELF3_PRR5_notTOC1, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 22 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col22 <- TF_common_LUX_LHY_CCA1_PRR7_ELF3_PRR5_notTOC1_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col22$clock <- "LUX and LHY and CCA1 and ELF3 and PRR7 and PRR5"

UpSet_trimmed_col22 <- UpSet_trimmed_col22 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col22 <- UpSet_trimmed_col22 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col23----
# LUX and LHY and TOC1 common targets - not targets of CCA1, ELF3, PRR7, PRR5 and ELF4
# firstly start with TF_common_LUX_LHY and do an inner_join with TOC1 

TF_common_LUX_LHY_TOC1 <- inner_join(TF_common_LUX_LHY[,1:2], TF_huang_merge_trimmed[,2], by="gene_ID")

# take this and rule out CCA1 targets
TF_common_LUX_LHY_TOC1_notCCA1 <- anti_join(TF_common_LUX_LHY_TOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LUX_LHY_TOC1_notCCA1_notELF3 <- anti_join(TF_common_LUX_LHY_TOC1_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_LHY_TOC1_notCCA1_notELF3_notPRR7 <- anti_join(TF_common_LUX_LHY_TOC1_notCCA1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_LHY_TOC1_notCCA1_notELF3_notPRR7_notPRR5 <- anti_join(TF_common_LUX_LHY_TOC1_notCCA1_notELF3_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_LHY_TOC1_notCCA1_notELF3_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LUX_LHY_TOC1_notCCA1_notELF3_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 23 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col23 <- TF_common_LUX_LHY_TOC1_notCCA1_notELF3_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col23$clock <- "LUX and LHY and TOC1"

UpSet_trimmed_col23 <- UpSet_trimmed_col23 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col23 <- UpSet_trimmed_col23 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col24----
# PRR7 and PRR5 common targets - not targets of LUX, LHY, TOC1, CCA1, ELF3 and ELF4
# # firstly an inner_join of PRR5 with PRR7 

TF_common_PRR7_PRR5 <- inner_join(TF_liu_merge_trimmed[,1:2], TF_nakamichi_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_common_PRR7_PRR5_notLUX <- anti_join(TF_common_PRR7_PRR5, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out LHY targets
TF_common_PRR7_PRR5_notLUX_notLHY <- anti_join(TF_common_PRR7_PRR5_notLUX, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_common_PRR7_PRR5_notLUX_notLHY_notTOC1 <- anti_join(TF_common_PRR7_PRR5_notLUX_notLHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_PRR7_PRR5_notLUX_notLHY_notTOC1_notCCA1 <- anti_join(TF_common_PRR7_PRR5_notLUX_notLHY_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_PRR7_PRR5_notLUX_notLHY_notTOC1_notCCA1_notELF3 <- anti_join(TF_common_PRR7_PRR5_notLUX_notLHY_notTOC1_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_PRR7_PRR5_notLUX_notLHY_notTOC1_notCCA1_notELF3_notELF4 <- anti_join(TF_common_PRR7_PRR5_notLUX_notLHY_notTOC1_notCCA1_notELF3, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 24 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col24 <- TF_common_PRR7_PRR5_notLUX_notLHY_notTOC1_notCCA1_notELF3_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col24$clock <- "PRR7 and PRR5"

UpSet_trimmed_col24 <- UpSet_trimmed_col24 %>% 
  left_join(PRR7_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col24 <- UpSet_trimmed_col24 %>% 
  left_join(PRR7_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col25----
# TOC1 and PRR7 common targets - not targets of LUX, LHY, CCA1, ELF3, PRR5 and ELF4
# # firstly an inner_join of PRR7 with TOC1 

TF_common_TOC1_PRR7 <- inner_join(TF_huang_merge_trimmed[,1:2], TF_liu_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_common_TOC1_PRR7_notLUX <- anti_join(TF_common_TOC1_PRR7, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out LHY targets
TF_common_TOC1_PRR7_notLUX_notLHY <- anti_join(TF_common_TOC1_PRR7_notLUX, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_TOC1_PRR7_notLUX_notLHY_notCCA1 <- anti_join(TF_common_TOC1_PRR7_notLUX_notLHY, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_TOC1_PRR7_notLUX_notLHY_notCCA1_notELF3 <- anti_join(TF_common_TOC1_PRR7_notLUX_notLHY_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_TOC1_PRR7_notLUX_notLHY_notCCA1_notELF3_notPRR5 <- anti_join(TF_common_TOC1_PRR7_notLUX_notLHY_notCCA1_notELF3, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_TOC1_PRR7_notLUX_notLHY_notCCA1_notELF3_notPRR5_notELF4 <- anti_join(TF_common_TOC1_PRR7_notLUX_notLHY_notCCA1_notELF3_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 25 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col25 <- TF_common_TOC1_PRR7_notLUX_notLHY_notCCA1_notELF3_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col25$clock <- "TOC1 and PRR7"

UpSet_trimmed_col25 <- UpSet_trimmed_col25 %>% 
  left_join(TOC1_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col25 <- UpSet_trimmed_col25 %>% 
  left_join(TOC1_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col26----
# LHY and TOC1 and PRR7 common targets - not targets of LUX, CCA1, ELF3, PRR5 and ELF4
# # firstly take TF_common_LHY_TOC1 and do an inner_join with PRR7 

TF_common_LHY_TOC1_PRR7 <- inner_join(TF_common_LHY_TOC1[,1:2], TF_liu_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_common_LHY_TOC1_PRR7_notLUX <- anti_join(TF_common_LHY_TOC1_PRR7, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LHY_TOC1_PRR7_notLUX_notCCA1 <- anti_join(TF_common_LHY_TOC1_PRR7_notLUX, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LHY_TOC1_PRR7_notLUX_notCCA1_notELF3 <- anti_join(TF_common_LHY_TOC1_PRR7_notLUX_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LHY_TOC1_PRR7_notLUX_notCCA1_notELF3_notPRR5 <- anti_join(TF_common_LHY_TOC1_PRR7_notLUX_notCCA1_notELF3, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LHY_TOC1_PRR7_notLUX_notCCA1_notELF3_notPRR5_notELF4 <- anti_join(TF_common_LHY_TOC1_PRR7_notLUX_notCCA1_notELF3_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 26 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col26 <- TF_common_LHY_TOC1_PRR7_notLUX_notCCA1_notELF3_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col26$clock <- "LHY and TOC1 and PRR7"

UpSet_trimmed_col26 <- UpSet_trimmed_col26 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col26 <- UpSet_trimmed_col26 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col27----
# LUX and ELF3 and PRR5 common targets - not targets of LHY, TOC1, CCA1, PRR7 and ELF4
# # firstly take TF_common_LUX_ELF3 and do an inner_join with PRR5 

TF_common_LUX_ELF3_PRR5 <- inner_join(TF_common_LUX_ELF3[,1:2], TF_nakamichi_merge_trimmed[,2], by="gene_ID")

# take this and rule out LHY targets
TF_common_LUX_ELF3_PRR5_notLHY <- anti_join(TF_common_LUX_ELF3_PRR5, TF_adams_merge_trimmed, by='gene_ID')

# take this and rule out TOC1 targets
TF_common_LUX_ELF3_PRR5_notLHY_notTOC1 <- anti_join(TF_common_LUX_ELF3_PRR5_notLHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_ELF3_PRR5_notLHY_notTOC1_notCCA1 <- anti_join(TF_common_LUX_ELF3_PRR5_notLHY_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_ELF3_PRR5_notLHY_notTOC1_notCCA1_notPRR7 <- anti_join(TF_common_LUX_ELF3_PRR5_notLHY_notTOC1_notCCA1, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_ELF3_PRR5_notLHY_notTOC1_notCCA1_notPRR7_notELF4 <- anti_join(TF_common_LUX_ELF3_PRR5_notLHY_notTOC1_notCCA1_notPRR7, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 27 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col27 <- TF_common_LUX_ELF3_PRR5_notLHY_notTOC1_notCCA1_notPRR7_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col27$clock <- "LUX and ELF3 and PRR5"

UpSet_trimmed_col27 <- UpSet_trimmed_col27 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col27 <- UpSet_trimmed_col27 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col28----
# LUX and CCA1 and ELF3 - not targets of LHY, TOC1, PRR7, PRR5 and ELF4
# # firstly take TF_common_LUX_ELF3 and do an inner_join with CCA1 

TF_common_LUX_ELF3_CCA1 <- inner_join(TF_common_LUX_ELF3[,1:2], TF_kamioka_nagel_merge_trimmed[,2], by="gene_ID")

# take this and rule out LHY targets
TF_common_LUX_ELF3_CCA1_notLHY <- anti_join(TF_common_LUX_ELF3_CCA1, TF_adams_merge_trimmed, by='gene_ID')

# take this and rule out TOC1 targets
TF_common_LUX_ELF3_CCA1_notLHY_notTOC1 <- anti_join(TF_common_LUX_ELF3_CCA1_notLHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_ELF3_CCA1_notLHY_notTOC1_notPRR7 <- anti_join(TF_common_LUX_ELF3_CCA1_notLHY_notTOC1, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_ELF3_CCA1_notLHY_notTOC1_notPRR7_notPRR5 <- anti_join(TF_common_LUX_ELF3_CCA1_notLHY_notTOC1_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_ELF3_CCA1_notLHY_notTOC1_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LUX_ELF3_CCA1_notLHY_notTOC1_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 28 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col28 <- TF_common_LUX_ELF3_CCA1_notLHY_notTOC1_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col28$clock <- "LUX and CCA1 and ELF3"

UpSet_trimmed_col28 <- UpSet_trimmed_col28 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col28 <- UpSet_trimmed_col28 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col29----
# LUX and TOC1 and ELF3 - not targets of LHY, CCA1, PRR7, PRR5 and ELF4
# # firstly take TF_common_LUX_TOC1 and do an inner_join with ELF3 

TF_common_LUX_TOC1_ELF3 <- inner_join(TF_common_LUX_TOC1[,1:2], TF_ezer_ELF3_merge_trimmed[,2], by="gene_ID")

# take this and rule out LHY targets
TF_common_LUX_TOC1_ELF3_notLHY <- anti_join(TF_common_LUX_TOC1_ELF3, TF_adams_merge_trimmed, by='gene_ID')

# take this and rule out CCA1 targets
TF_common_LUX_TOC1_ELF3_notLHY_notCCA1 <- anti_join(TF_common_LUX_TOC1_ELF3_notLHY, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_TOC1_ELF3_notLHY_notCCA1_notPRR7 <- anti_join(TF_common_LUX_TOC1_ELF3_notLHY_notCCA1, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_TOC1_ELF3_notLHY_notCCA1_notPRR7_notPRR5 <- anti_join(TF_common_LUX_TOC1_ELF3_notLHY_notCCA1_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_TOC1_ELF3_notLHY_notCCA1_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LUX_TOC1_ELF3_notLHY_notCCA1_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 29 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col29 <- TF_common_LUX_TOC1_ELF3_notLHY_notCCA1_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col29$clock <- "LUX and TOC1 and ELF3"

UpSet_trimmed_col29 <- UpSet_trimmed_col29 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col29 <- UpSet_trimmed_col29 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col30----
# LUX and TOC1 and CCA1 - not targets of LHY, ELF3, PRR7, PRR5 and ELF4
# # firstly take TF_common_LUX_TOC1 and do an inner_join with CCA1 

TF_common_LUX_TOC1_CCA1 <- inner_join(TF_common_LUX_TOC1[,1:2], TF_kamioka_nagel_merge_trimmed[,2], by="gene_ID")

# take this and rule out LHY targets
TF_common_LUX_TOC1_CCA1_notLHY <- anti_join(TF_common_LUX_TOC1_CCA1, TF_adams_merge_trimmed, by='gene_ID')

# take this and rule out ELF3 targets
TF_common_LUX_TOC1_CCA1_notLHY_notELF3 <- anti_join(TF_common_LUX_TOC1_CCA1_notLHY, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_TOC1_CCA1_notLHY_notELF3_notPRR7 <- anti_join(TF_common_LUX_TOC1_CCA1_notLHY_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_TOC1_CCA1_notLHY_notELF3_notPRR7_notPRR5 <- anti_join(TF_common_LUX_TOC1_CCA1_notLHY_notELF3_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_TOC1_CCA1_notLHY_notELF3_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LUX_TOC1_CCA1_notLHY_notELF3_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 30 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col30 <- TF_common_LUX_TOC1_CCA1_notLHY_notELF3_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col30$clock <- "LUX and TOC1 and CCA1"

UpSet_trimmed_col30 <- UpSet_trimmed_col30 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col30 <- UpSet_trimmed_col30 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col31----
# LUX and LHY and ELF3 and PRR7 and PRR5 and ELF4 - not targets of TOC1 and CCA1
# # firstly take TF_common_LUX_LHY_ELF3 and do an inner_join with TF_common_PRR7_PRR5 

TF_common_LUX_LHY_ELF3_PRR7_PRR5 <- inner_join(TF_common_LUX_LHY_ELF3[,1:2], TF_common_PRR7_PRR5[,2], by="gene_ID")

# take this and include ELF4 targets
TF_common_LUX_LHY_ELF3_PRR7_PRR5_ELF4 <- inner_join(TF_common_LUX_LHY_ELF3_PRR7_PRR5[,1:2], TF_ezer_ELF4_merge_trimmed[,2], by='gene_ID')

# take this and rule out TOC1 targets
TF_common_LUX_LHY_ELF3_PRR7_PRR5_ELF4_notTOC1 <- anti_join(TF_common_LUX_LHY_ELF3_PRR7_PRR5_ELF4, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_LHY_ELF3_PRR7_PRR5_ELF4_notTOC1_notCCA1 <- anti_join(TF_common_LUX_LHY_ELF3_PRR7_PRR5_ELF4_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 31 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col31 <- TF_common_LUX_LHY_ELF3_PRR7_PRR5_ELF4_notTOC1_notCCA1[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col31$clock <- "LUX and LHY and ELF3 and PRR7 and PRR5 and ELF4"

UpSet_trimmed_col31 <- UpSet_trimmed_col31 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col31 <- UpSet_trimmed_col31 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col32----
# LUX and LHY and CCA1 and ELF3 common targets - not targets of TOC1, PRR7, PRR5 and ELF4
# # firstly take TF_common_LUX_LHY_CCA1 and do an inner_join with ELF3 

TF_common_LUX_LHY_CCA1_ELF3 <- inner_join(TF_common_LUX_LHY_CCA1[,1:2], TF_ezer_ELF3_merge_trimmed[,2], by="gene_ID")

# take this and rule out TOC1 targets
TF_common_LUX_LHY_CCA1_ELF3_notTOC1 <- anti_join(TF_common_LUX_LHY_CCA1_ELF3, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_LHY_CCA1_ELF3_notTOC1_notPRR7 <- anti_join(TF_common_LUX_LHY_CCA1_ELF3_notTOC1, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_LHY_CCA1_ELF3_notTOC1_notPRR7_notPRR5 <- anti_join(TF_common_LUX_LHY_CCA1_ELF3_notTOC1_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_LHY_CCA1_ELF3_notTOC1_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LUX_LHY_CCA1_ELF3_notTOC1_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 32 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col32 <- TF_common_LUX_LHY_CCA1_ELF3_notTOC1_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col32$clock <- "LUX and LHY and CCA1 and ELF3"

UpSet_trimmed_col32 <- UpSet_trimmed_col32 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col32 <- UpSet_trimmed_col32 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col33----
# LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and PRR5 common targets - not targets of ELF4
# # firstly take TF_common_LUX_LHY_TOC1 and do an inner_join with CCA1, ELF3, PRR7 and PRR5 

TF_common_LUX_LHY_TOC1_CCA1 <- inner_join(TF_common_LUX_LHY_TOC1[,1:2], TF_kamioka_nagel_merge_trimmed[,2], by="gene_ID")

# take this and include ELF3 targets
TF_common_LUX_LHY_TOC1_CCA1_ELF3 <- inner_join(TF_common_LUX_LHY_TOC1_CCA1[,1:2], TF_ezer_ELF3_merge_trimmed[,2], by='gene_ID')

# take previous and include PRR7 targets
TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7 <- inner_join(TF_common_LUX_LHY_TOC1_CCA1_ELF3[,1:2], TF_liu_merge_trimmed[,2], by='gene_ID')

# take previous and include PRR5 targets
TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_PRR5 <- inner_join(TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7[,1:2], TF_nakamichi_merge_trimmed[,2], by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_PRR5_notELF4 <- anti_join(TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_PRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 33 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col33 <- TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_PRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col33$clock <- "LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and PRR5"

UpSet_trimmed_col33 <- UpSet_trimmed_col33 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col33 <- UpSet_trimmed_col33 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col34----
# LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and PRR5 and ELF4 common targets
# firstly take TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_PRR5 and do an inner_join with ELF4 

TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_PRR5_ELF4 <- inner_join(TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_PRR5[,1:2], TF_ezer_ELF4_merge_trimmed[,2], by="gene_ID")

# tidy up dataset
# select first two columns
# Column 34 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col34 <- TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_PRR5_ELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col34$clock <- "LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and PRR5 and ELF4"

UpSet_trimmed_col34 <- UpSet_trimmed_col34 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col34 <- UpSet_trimmed_col34 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col35----
# ELF3 and ELF4 common targets - not targets of LUX, LHY, TOC1, CCA1, PRR7 and PRR5
# firstly an inner_join of ELF3 with ELF4 

TF_ELF3_ELF4 <- inner_join(TF_ezer_ELF3_merge_trimmed[,1:2], TF_ezer_ELF4_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_ELF3_ELF4_notLUX <- anti_join(TF_ELF3_ELF4, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out LHY targets
TF_ELF3_ELF4_notLUX_notLHY <- anti_join(TF_ELF3_ELF4_notLUX, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_ELF3_ELF4_notLUX_notLHY_notTOC1 <- anti_join(TF_ELF3_ELF4_notLUX_notLHY, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_ELF3_ELF4_notLUX_notLHY_notTOC1_notCCA1 <- anti_join(TF_ELF3_ELF4_notLUX_notLHY_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_ELF3_ELF4_notLUX_notLHY_notTOC1_notCCA1_notPRR7 <- anti_join(TF_ELF3_ELF4_notLUX_notLHY_notTOC1_notCCA1, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_ELF3_ELF4_notLUX_notLHY_notTOC1_notCCA1_notPRR7_notPRR5 <- anti_join(TF_ELF3_ELF4_notLUX_notLHY_notTOC1_notCCA1_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 35 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col35 <- TF_ELF3_ELF4_notLUX_notLHY_notTOC1_notCCA1_notPRR7_notPRR5[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col35$clock <- "ELF3 and ELF4"

UpSet_trimmed_col35 <- UpSet_trimmed_col35 %>% 
  left_join(ELF3_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col35 <- UpSet_trimmed_col35 %>% 
  left_join(ELF3_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col36----
# TOC1 and PRR5 common targets - not targets of LUX, LHY, CCA1, ELF3, PRR7 and ELF4
# firstly an inner_join of TOC1 with PRR5 

TF_TOC1_PRR5 <- inner_join(TF_huang_merge_trimmed[,1:2], TF_nakamichi_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_TOC1_PRR5_notLUX <- anti_join(TF_TOC1_PRR5, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out LHY targets
TF_TOC1_PRR5_notLUX_notLHY <- anti_join(TF_TOC1_PRR5_notLUX, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_TOC1_PRR5_notLUX_notLHY_notCCA1 <- anti_join(TF_TOC1_PRR5_notLUX_notLHY, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_TOC1_PRR5_notLUX_notLHY_notCCA1_notELF3 <- anti_join(TF_TOC1_PRR5_notLUX_notLHY_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_TOC1_PRR5_notLUX_notLHY_notCCA1_notELF3_notPRR7 <- anti_join(TF_TOC1_PRR5_notLUX_notLHY_notCCA1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_TOC1_PRR5_notLUX_notLHY_notCCA1_notELF3_notPRR7_notELF4 <- anti_join(TF_TOC1_PRR5_notLUX_notLHY_notCCA1_notELF3_notPRR7, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 36 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col36 <- TF_TOC1_PRR5_notLUX_notLHY_notCCA1_notELF3_notPRR7_notELF4[,1:2]


# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col36$clock <- "TOC1 and PRR5"

UpSet_trimmed_col36 <- UpSet_trimmed_col36 %>% 
  left_join(TOC1_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col36 <- UpSet_trimmed_col36 %>% 
  left_join(TOC1_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col37----
# TOC1 and PRR7 and PRR5 common targets - not targets of LUX, LHY, CCA1, ELF3 and ELF4
# firstly take TF_common_TOC1_PRR7 and do an inner_join with PRR5 

TF_common_TOC1_PRR7_PRR5 <- inner_join(TF_common_TOC1_PRR7[,1:2], TF_nakamichi_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_common_TOC1_PRR7_PRR5_notLUX <- anti_join(TF_common_TOC1_PRR7_PRR5, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out LHY targets
TF_common_TOC1_PRR7_PRR5_notLUX_notLHY <- anti_join(TF_common_TOC1_PRR7_PRR5_notLUX, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_TOC1_PRR7_PRR5_notLUX_notLHY_notCCA1 <- anti_join(TF_common_TOC1_PRR7_PRR5_notLUX_notLHY, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_TOC1_PRR7_PRR5_notLUX_notLHY_notCCA1_notELF3 <- anti_join(TF_common_TOC1_PRR7_PRR5_notLUX_notLHY_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_TOC1_PRR7_PRR5_notLUX_notLHY_notCCA1_notELF3_notELF4 <- anti_join(TF_common_TOC1_PRR7_PRR5_notLUX_notLHY_notCCA1_notELF3, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 37 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col37 <- TF_common_TOC1_PRR7_PRR5_notLUX_notLHY_notCCA1_notELF3_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col37$clock <- "TOC1 and PRR7 and PRR5"

UpSet_trimmed_col37 <- UpSet_trimmed_col37 %>% 
  left_join(TOC1_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col37 <- UpSet_trimmed_col37 %>% 
  left_join(TOC1_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col38----
# TOC1 and CCA1 common targets - not targets of LUX, LHY, ELF3, PRR7, PRR5 and ELF4
# firstly an inner_join of TOC1 with CCA1

TF_common_TOC1_CCA1 <- inner_join(TF_huang_merge_trimmed[,1:2], TF_kamioka_nagel_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_common_TOC1_CCA1_notLUX <- anti_join(TF_common_TOC1_CCA1, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out LHY targets
TF_common_TOC1_CCA1_notLUX_notLHY <- anti_join(TF_common_TOC1_CCA1_notLUX, TF_adams_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_TOC1_CCA1_notLUX_notLHY_notELF3 <- anti_join(TF_common_TOC1_CCA1_notLUX_notLHY, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_TOC1_CCA1_notLUX_notLHY_notELF3_notPRR7 <- anti_join(TF_common_TOC1_CCA1_notLUX_notLHY_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_TOC1_CCA1_notLUX_notLHY_notELF3_notPRR7_notPRR5 <- anti_join(TF_common_TOC1_CCA1_notLUX_notLHY_notELF3_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_TOC1_CCA1_notLUX_notLHY_notELF3_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_TOC1_CCA1_notLUX_notLHY_notELF3_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 38 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col38 <- TF_common_TOC1_CCA1_notLUX_notLHY_notELF3_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col38$clock <- "TOC1 and CCA1"

UpSet_trimmed_col38 <- UpSet_trimmed_col38 %>% 
  left_join(TOC1_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col38 <- UpSet_trimmed_col38 %>% 
  left_join(TOC1_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col39----
# LHY and PRR7 common targets - not targets of LUX, TOC1, CCA1, ELF3, PRR5 and ELF4
# firstly an inner_join of LHY with PRR7

TF_common_LHY_PRR7 <- inner_join(TF_adams_merge_trimmed[,1:2], TF_liu_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_common_LHY_PRR7_notLUX <- anti_join(TF_common_LHY_PRR7, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_common_LHY_PRR7_notLUX_notTOC1 <- anti_join(TF_common_LHY_PRR7_notLUX, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LHY_PRR7_notLUX_notTOC1_notCCA1 <- anti_join(TF_common_LHY_PRR7_notLUX_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LHY_PRR7_notLUX_notTOC1_notCCA1_notELF3 <- anti_join(TF_common_LHY_PRR7_notLUX_notTOC1_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LHY_PRR7_notLUX_notTOC1_notCCA1_notELF3_notPRR5 <- anti_join(TF_common_LHY_PRR7_notLUX_notTOC1_notCCA1_notELF3, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LHY_PRR7_notLUX_notTOC1_notCCA1_notELF3_notPRR5_notELF4 <- anti_join(TF_common_LHY_PRR7_notLUX_notTOC1_notCCA1_notELF3_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 39 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col39 <- TF_common_LHY_PRR7_notLUX_notTOC1_notCCA1_notELF3_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col39$clock <- "LHY and PRR7"

UpSet_trimmed_col39 <- UpSet_trimmed_col39 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col39 <- UpSet_trimmed_col39 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col40----
# LHY and PRR7 and PRR5 common targets - not targets of LUX, TOC1, CCA1, ELF3 and ELF4
# firstly take TF_common_LHY_PRR7 an then do an inner_join with PRR5

TF_common_LHY_PRR7_PRR5 <- inner_join(TF_common_LHY_PRR7[,1:2], TF_nakamichi_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_common_LHY_PRR7_PRR5_notLUX <- anti_join(TF_common_LHY_PRR7_PRR5, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_common_LHY_PRR7_PRR5_notLUX_notTOC1 <- anti_join(TF_common_LHY_PRR7_PRR5_notLUX, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LHY_PRR7_PRR5_notLUX_notTOC1_notCCA1 <- anti_join(TF_common_LHY_PRR7_PRR5_notLUX_notTOC1, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LHY_PRR7_PRR5_notLUX_notTOC1_notCCA1_notELF3 <- anti_join(TF_common_LHY_PRR7_PRR5_notLUX_notTOC1_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LHY_PRR7_PRR5_notLUX_notTOC1_notCCA1_notELF3_notELF4 <- anti_join(TF_common_LHY_PRR7_PRR5_notLUX_notTOC1_notCCA1_notELF3, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 40 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col40 <- TF_common_LHY_PRR7_PRR5_notLUX_notTOC1_notCCA1_notELF3_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col40$clock <- "LHY and PRR7 and PRR5"

UpSet_trimmed_col40 <- UpSet_trimmed_col40 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col40 <- UpSet_trimmed_col40 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col41----
# ELF4 targets alone i.e. not targets of LUX, LHY, CCA1, TOC1, CCA1, ELF3, PRR7, and PRR5 
# An anti_join of ELF4 and LUX gives 1 obs

TF_ELF4 <- anti_join(TF_ezer_ELF4_merge_trimmed, TF_ezer_LUX_merge_trimmed, by="gene_ID")

# tidy up dataset
# select first two columns
# Column 41 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col41 <- TF_ELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col41$clock <- "ELF4"

UpSet_trimmed_col41 <- UpSet_trimmed_col41 %>% 
  left_join(ELF4_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col41 <- UpSet_trimmed_col41 %>% 
  left_join(ELF4_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col42----
# LHY and CCA1 and PRR7 - not targets of LUX, TOC1, ELF3, PRR5 and ELF4
# firstly take TF_common_LHY_CCA1 and do an inner_join with PRR7 

TF_common_LHY_CCA1_PRR7 <- inner_join(TF_common_LHY_CCA1[,1:2], TF_liu_merge_trimmed[,2], by="gene_ID")

# take this and rule out LUX targets
TF_common_LHY_CCA1_PRR7_notLUX <- anti_join(TF_common_LHY_CCA1_PRR7, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_common_LHY_CCA1_PRR7_notLUX_notTOC1 <- anti_join(TF_common_LHY_CCA1_PRR7_notLUX, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LHY_CCA1_PRR7_notLUX_notTOC1_not_ELF3 <- anti_join(TF_common_LHY_CCA1_PRR7_notLUX_notTOC1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LHY_CCA1_PRR7_notLUX_notTOC1_not_ELF3_notPRR5 <- anti_join(TF_common_LHY_CCA1_PRR7_notLUX_notTOC1_not_ELF3, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LHY_CCA1_PRR7_notLUX_notTOC1_not_ELF3_notPRR5_notELF4 <- anti_join(TF_common_LHY_CCA1_PRR7_notLUX_notTOC1_not_ELF3_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 42 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col42 <- TF_common_LHY_CCA1_PRR7_notLUX_notTOC1_not_ELF3_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col42$clock <- "LHY and CCA1 and PRR7"

UpSet_trimmed_col42 <- UpSet_trimmed_col42 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col42 <- UpSet_trimmed_col42 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col43----
# LHY and TOC1 and CCA1 and PRR7 and PRR5 - not targets of LUX, ELF3 and ELF4
# firstly take TF_common_LHY_TOC1_CCA1 and do an inner join with PRR7 

TF_common_LHY_TOC1_CCA1_PRR7 <- inner_join(TF_common_LHY_TOC1_CCA1[,1:2], TF_liu_merge_trimmed[,2], by="gene_ID")

# take this and include PRR5 targets
TF_common_LHY_TOC1_CCA1_PRR7_PRR5 <- inner_join(TF_common_LHY_TOC1_CCA1_PRR7[,1:2], TF_nakamichi_merge_trimmed[,2], by='gene_ID')

# take previous and rule out LUX targets
TF_common_LHY_TOC1_CCA1_PRR7_PRR5_notLUX <- anti_join(TF_common_LHY_TOC1_CCA1_PRR7_PRR5, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LHY_TOC1_CCA1_PRR7_PRR5_notLUX_not_ELF3 <- anti_join(TF_common_LHY_TOC1_CCA1_PRR7_PRR5_notLUX, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LHY_TOC1_CCA1_PRR7_PRR5_notLUX_not_ELF3_notELF4 <- anti_join(TF_common_LHY_TOC1_CCA1_PRR7_PRR5_notLUX_not_ELF3, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 43 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col43 <- TF_common_LHY_TOC1_CCA1_PRR7_PRR5_notLUX_not_ELF3_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col43$clock <- "LHY and TOC1 and CCA1 and PRR7 and PRR5"

UpSet_trimmed_col43 <- UpSet_trimmed_col43 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col43 <- UpSet_trimmed_col43 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col44----
# LHY and CCA1 and ELF3 - not targets of LUX, TOC1, PRR7, PRR5 and ELF4
# firstly take TF_common_LHY_CCA1 and do an inner join with ELF3 

TF_common_LHY_CCA1_ELF3 <- inner_join(TF_common_LHY_CCA1[,1:2], TF_ezer_ELF3_merge_trimmed[,2], by="gene_ID")

# take previous and rule out LUX targets
TF_common_LHY_CCA1_ELF3_notLUX <- anti_join(TF_common_LHY_CCA1_ELF3, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_common_LHY_CCA1_ELF3_notLUX_notTOC1 <- anti_join(TF_common_LHY_CCA1_ELF3_notLUX, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LHY_CCA1_ELF3_notLUX_notTOC1_notPRR7 <- anti_join(TF_common_LHY_CCA1_ELF3_notLUX_notTOC1, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LHY_CCA1_ELF3_notLUX_notTOC1_notPRR7_notPRR5 <- anti_join(TF_common_LHY_CCA1_ELF3_notLUX_notTOC1_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LHY_CCA1_ELF3_notLUX_notTOC1_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LHY_CCA1_ELF3_notLUX_notTOC1_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 44 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col44 <- TF_common_LHY_CCA1_ELF3_notLUX_notTOC1_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col44$clock <- "LHY and CCA1 and ELF3"

UpSet_trimmed_col44 <- UpSet_trimmed_col44 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col44 <- UpSet_trimmed_col44 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col45----
# LUX and LHY and TOC1 and CCA1 - not targets of ELF3, PRR7, PRR5 and ELF4
# firstly take TF_common_LUX_LHY_TOC1 and do an inner join with CCA1 

TF_common_LUX_LHY_TOC1_CCA1 <- inner_join(TF_common_LUX_LHY_TOC1[,1:2], TF_kamioka_nagel_merge_trimmed[,2], by="gene_ID")

# take previous and rule out ELF3 targets
TF_common_LUX_LHY_TOC1_CCA1_notELF3 <- anti_join(TF_common_LUX_LHY_TOC1_CCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_LHY_TOC1_CCA1_notELF3_notPRR7 <- anti_join(TF_common_LUX_LHY_TOC1_CCA1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_LHY_TOC1_CCA1_notELF3_notPRR7_notPRR5 <- anti_join(TF_common_LUX_LHY_TOC1_CCA1_notELF3_notPRR7, TF_nakamichi_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_LHY_TOC1_CCA1_notELF3_notPRR7_notPRR5_notELF4 <- anti_join(TF_common_LUX_LHY_TOC1_CCA1_notELF3_notPRR7_notPRR5, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 45 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col45 <- TF_common_LUX_LHY_TOC1_CCA1_notELF3_notPRR7_notPRR5_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col45$clock <- "LUX and LHY and TOC1 and CCA1"

UpSet_trimmed_col45 <- UpSet_trimmed_col45 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col45 <- UpSet_trimmed_col45 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col46----
# LHY and TOC1 and PRR7 and PRR5 - not targets of LUX, CCA1, ELF3 and ELF4
# firstly take TF_common_LHY_TOC1_PRR7 and do an inner join with PRR5 

TF_common_LHY_TOC1_PRR7_PRR5 <- inner_join(TF_common_LHY_TOC1_PRR7[,1:2], TF_nakamichi_merge_trimmed[,2], by="gene_ID")

# take previous and rule out LUX targets
TF_common_LHY_TOC1_PRR7_PRR5_notLUX <- anti_join(TF_common_LHY_TOC1_PRR7_PRR5, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LHY_TOC1_PRR7_PRR5_notLUX_notCCA1 <- anti_join(TF_common_LHY_TOC1_PRR7_PRR5_notLUX, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LHY_TOC1_PRR7_PRR5_notLUX_notCCA1_notELF3 <- anti_join(TF_common_LHY_TOC1_PRR7_PRR5_notLUX_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LHY_TOC1_PRR7_PRR5_notLUX_notCCA1_notELF3_notELF4 <- anti_join(TF_common_LHY_TOC1_PRR7_PRR5_notLUX_notCCA1_notELF3, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 46 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col46 <- TF_common_LHY_TOC1_PRR7_PRR5_notLUX_notCCA1_notELF3_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col46$clock <- "LHY and TOC1 and PRR7 and PRR5"

UpSet_trimmed_col46 <- UpSet_trimmed_col46 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col46 <- UpSet_trimmed_col46 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col47----
# LHY and CCA1 and PRR5 - not targets of LUX, TOC1, ELF3, PRR7 and ELF4
# firstly take TF_common_LHY_CCA1 and do an inner join with PRR5 

TF_common_LHY_CCA1_PRR5 <- inner_join(TF_common_LHY_CCA1[,1:2], TF_nakamichi_merge_trimmed[,2], by="gene_ID")

# take previous and rule out LUX targets
TF_common_LHY_CCA1_PRR5_notLUX <- anti_join(TF_common_LHY_CCA1_PRR5, TF_ezer_LUX_merge_trimmed, by='gene_ID')

# take previous and rule out TOC1 targets
TF_common_LHY_CCA1_PRR5_notLUX_notTOC1 <- anti_join(TF_common_LHY_CCA1_PRR5_notLUX, TF_huang_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LHY_CCA1_PRR5_notLUX_notTOC1_notELF3 <- anti_join(TF_common_LHY_CCA1_PRR5_notLUX_notTOC1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LHY_CCA1_PRR5_notLUX_notTOC1_notELF3_notPRR7 <- anti_join(TF_common_LHY_CCA1_PRR5_notLUX_notTOC1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LHY_CCA1_PRR5_notLUX_notTOC1_notELF3_notPRR7_notELF4 <- anti_join(TF_common_LHY_CCA1_PRR5_notLUX_notTOC1_notELF3_notPRR7, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 47 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col47 <- TF_common_LHY_CCA1_PRR5_notLUX_notTOC1_notELF3_notPRR7_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col47$clock <- "LHY and CCA1 and PRR5"

UpSet_trimmed_col47 <- UpSet_trimmed_col47 %>% 
  left_join(LHY_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col47 <- UpSet_trimmed_col47 %>% 
  left_join(LHY_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col48----
# LUX and LHY and TOC1 and PRR5 - not targets of CCA1, ELF3, PRR7 and ELF4
# firstly take TF_common_LUX_LHY_TOC1 and do an inner join with PRR5 

TF_common_LUX_LHY_TOC1_PRR5 <- inner_join(TF_common_LUX_LHY_TOC1[,1:2], TF_nakamichi_merge_trimmed[,2], by="gene_ID")

# take previous and rule out CCA1 targets
TF_common_LUX_LHY_TOC1_PRR5_notCCA1 <- anti_join(TF_common_LUX_LHY_TOC1_PRR5, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# take previous and rule out ELF3 targets
TF_common_LUX_LHY_TOC1_PRR5_notCCA1_notELF3 <- anti_join(TF_common_LUX_LHY_TOC1_PRR5_notCCA1, TF_ezer_ELF3_merge_trimmed, by='gene_ID')

# take previous and rule out PRR7 targets
TF_common_LUX_LHY_TOC1_PRR5_notCCA1_notELF3_notPRR7 <- anti_join(TF_common_LUX_LHY_TOC1_PRR5_notCCA1_notELF3, TF_liu_merge_trimmed, by='gene_ID')

# take previous and rule out ELF4 targets
TF_common_LUX_LHY_TOC1_PRR5_notCCA1_notELF3_notPRR7_notELF4 <- anti_join(TF_common_LUX_LHY_TOC1_PRR5_notCCA1_notELF3_notPRR7, TF_ezer_ELF4_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 48 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col48 <- TF_common_LUX_LHY_TOC1_PRR5_notCCA1_notELF3_notPRR7_notELF4[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col48$clock <- "LUX and LHY and TOC1 and PRR5"

UpSet_trimmed_col48 <- UpSet_trimmed_col48 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col48 <- UpSet_trimmed_col48 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col49----
# LUX and LHY and TOC1 and ELF3 and PRR7 and PRR5 and ELF4 - not targets of CCA1
# firstly take TF_common_LUX_LHY_TOC1 and do an inner join with ELF3 

TF_common_LUX_LHY_TOC1_ELF3 <- inner_join(TF_common_LUX_LHY_TOC1[,1:2], TF_ezer_ELF3_merge_trimmed[,2], by="gene_ID")

# take previous and include PRR7 targets
TF_common_LUX_LHY_TOC1_ELF3_PRR7 <- inner_join(TF_common_LUX_LHY_TOC1_ELF3[,1:2], TF_liu_merge_trimmed[,2], by='gene_ID')

# take previous and include PRR5 targets
TF_common_LUX_LHY_TOC1_ELF3_PRR7_PRR5 <- inner_join(TF_common_LUX_LHY_TOC1_ELF3_PRR7[,1:2], TF_nakamichi_merge_trimmed[,2], by='gene_ID')

# take previous and include ELF4 targets
TF_common_LUX_LHY_TOC1_ELF3_PRR7_PRR5_ELF4 <- inner_join(TF_common_LUX_LHY_TOC1_ELF3_PRR7_PRR5[,1:2], TF_ezer_ELF4_merge_trimmed[,2], by='gene_ID')

# take previous and rule out CCA1 targets
TF_common_LUX_LHY_TOC1_ELF3_PRR7_PRR5_ELF4_notCCA1 <- anti_join(TF_common_LUX_LHY_TOC1_ELF3_PRR7_PRR5_ELF4, TF_kamioka_nagel_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 49 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col49 <- TF_common_LUX_LHY_TOC1_ELF3_PRR7_PRR5_ELF4_notCCA1[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col49$clock <- "LUX and LHY and TOC1 and ELF3 and PRR7 and PRR5 and ELF4"

UpSet_trimmed_col49 <- UpSet_trimmed_col49 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col49 <- UpSet_trimmed_col49 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# col50----
# LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and ELF4 - not targets of PRR5
# firstly take TF_common_LUX_LHY_TOC1 and do an inner join with CCA1 

TF_common_LUX_LHY_TOC1_CCA1 <- inner_join(TF_common_LUX_LHY_TOC1[,1:2], TF_kamioka_nagel_merge_trimmed[,2], by="gene_ID")

# take previous and include ELF3 targets
TF_common_LUX_LHY_TOC1_CCA1_ELF3 <- inner_join(TF_common_LUX_LHY_TOC1_CCA1[,1:2], TF_ezer_ELF3_merge_trimmed[,2], by='gene_ID')

# take previous and include PRR7 targets
TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7 <- inner_join(TF_common_LUX_LHY_TOC1_CCA1_ELF3[,1:2], TF_liu_merge_trimmed[,2], by='gene_ID')

# take previous and include ELF4 targets
TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_ELF4 <- inner_join(TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7[,1:2], TF_ezer_ELF4_merge_trimmed[,2], by='gene_ID')

# take previous and rule out PRR5 targets
TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_ELF4_notPRR5 <- anti_join(TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_ELF4, TF_nakamichi_merge_trimmed, by='gene_ID')

# tidy up dataset
# select first two columns
# Column 50 of UpSetR_trimmed plot for 8 clock components
UpSet_trimmed_col50 <- TF_common_LUX_LHY_TOC1_CCA1_ELF3_PRR7_ELF4_notPRR5[,1:2]

# add a column describing clock TFs targeting gene_IDs
UpSet_trimmed_col50$clock <- "LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and ELF4"

UpSet_trimmed_col50 <- UpSet_trimmed_col50 %>% 
  left_join(LUX_bind_d1d2, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-4) %>% 
  dplyr::rename(type_d1d2 = type)

UpSet_trimmed_col50 <- UpSet_trimmed_col50 %>% 
  left_join(LUX_bind_d1d5, by = c('gene_ID', 'cluster')) %>% 
  dplyr::select(-5) %>% 
  dplyr::rename(type_d1d5 = type)

# UpSet plot overlap summaries
get_summary <- function(df, group_col, clock_id, col_number){
  
  s <- df %>% 
    group_by({{ group_col }}) %>%
    dplyr::summarise(count = n()) %>% mutate(percent = (count/sum(count)) * 100) %>% 
    mutate(clock = {{ clock_id }}) %>% 
    mutate(column = {{ col_number }})
}

# d1d2----
col1_summary_d1d2 <- get_summary(UpSet_trimmed_col1, type_d1d2, 'LHY', 1)
col2_summary_d1d2 <- get_summary(UpSet_trimmed_col2, type_d1d2, 'TOC1', 2)
col3_summary_d1d2 <- get_summary(UpSet_trimmed_col3, type_d1d2, 'LUX', 3)
col4_summary_d1d2 <- get_summary(UpSet_trimmed_col4, type_d1d2, 'LHY and CCA1', 4)
col5_summary_d1d2 <- get_summary(UpSet_trimmed_col5, type_d1d2, 'LUX and ELF3', 5)
col6_summary_d1d2 <- get_summary(UpSet_trimmed_col6, type_d1d2, 'PRR5', 6)
col7_summary_d1d2 <- get_summary(UpSet_trimmed_col7, type_d1d2, 'PRR7', 7)
col8_summary_d1d2 <- get_summary(UpSet_trimmed_col8, type_d1d2, 'LUX and LHY', 8)
col9_summary_d1d2 <- get_summary(UpSet_trimmed_col9, type_d1d2, 'LUX and LHY and CCA1', 9)
col10_summary_d1d2 <- get_summary(UpSet_trimmed_col10, type_d1d2, 'LUX and TOC1', 10)

col11_summary_d1d2 <- get_summary(UpSet_trimmed_col11, type_d1d2, 'CCA1', 11)
col12_summary_d1d2 <- get_summary(UpSet_trimmed_col12, type_d1d2, 'LUX and ELF3 and ELF4', 12)
col13_summary_d1d2 <- get_summary(UpSet_trimmed_col13, type_d1d2, 'LHY and TOC1', 13)
col14_summary_d1d2 <- get_summary(UpSet_trimmed_col14, type_d1d2, 'LHY and TOC1 and CCA1', 14)
col15_summary_d1d2 <- get_summary(UpSet_trimmed_col15, type_d1d2, 'LUX and TOC1 and PRR5', 15)
col16_summary_d1d2 <- get_summary(UpSet_trimmed_col16, type_d1d2, 'LUX and LHY and ELF3', 16)
col17_summary_d1d2 <- get_summary(UpSet_trimmed_col17, type_d1d2, 'LUX and LHY and CCA1 and PRR7', 17)
col18_summary_d1d2 <- get_summary(UpSet_trimmed_col18, type_d1d2, 'ELF3', 18)
col19_summary_d1d2 <- get_summary(UpSet_trimmed_col19, type_d1d2, 'LUX and ELF3 and PRR5 and ELF4', 19)
col20_summary_d1d2 <- get_summary(UpSet_trimmed_col20, type_d1d2, 'LUX and TOC1 and ELF3 and PRR7 and PRR5 and ELF4', 20)

col21_summary_d1d2 <- get_summary(UpSet_trimmed_col21, type_d1d2, 'LUX and LHY and ELF3 and ELF4', 21)
col22_summary_d1d2 <- get_summary(UpSet_trimmed_col22, type_d1d2, 'LUX and LHY and CCA1 and ELF3 and PRR7 and PRR5', 22)
col23_summary_d1d2 <- get_summary(UpSet_trimmed_col23, type_d1d2, 'LUX and LHY and TOC1', 23)
col24_summary_d1d2 <- get_summary(UpSet_trimmed_col24, type_d1d2, 'PRR7 and PRR5', 24)
col25_summary_d1d2 <- get_summary(UpSet_trimmed_col25, type_d1d2, 'TOC1 and PRR7', 25)
col26_summary_d1d2 <- get_summary(UpSet_trimmed_col26, type_d1d2, 'LHY and TOC1 and PRR7', 26)
col27_summary_d1d2 <- get_summary(UpSet_trimmed_col27, type_d1d2, 'LUX and ELF3 and PRR5', 27)
col28_summary_d1d2 <- get_summary(UpSet_trimmed_col28, type_d1d2, 'LUX and CCA1 and ELF3', 28)
col29_summary_d1d2 <- get_summary(UpSet_trimmed_col29, type_d1d2, 'LUX and TOC1 and ELF3', 29)
col30_summary_d1d2 <- get_summary(UpSet_trimmed_col30, type_d1d2, 'LUX and TOC1 and CCA1', 30)

col31_summary_d1d2 <- get_summary(UpSet_trimmed_col31, type_d1d2, 'LUX and LHY and ELF3 and PRR7 and PRR5 and ELF4', 31)
col32_summary_d1d2 <- get_summary(UpSet_trimmed_col32, type_d1d2, 'LUX and LHY and CCA1 and ELF3', 32)
col33_summary_d1d2 <- get_summary(UpSet_trimmed_col33, type_d1d2, 'LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and PRR5', 33)
col34_summary_d1d2 <- get_summary(UpSet_trimmed_col34, type_d1d2, 'LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and PRR5 and ELF4', 34)
col35_summary_d1d2 <- get_summary(UpSet_trimmed_col35, type_d1d2, 'ELF3 and ELF4', 35)
col36_summary_d1d2 <- get_summary(UpSet_trimmed_col36, type_d1d2, 'TOC1 and PRR5', 36)
col37_summary_d1d2 <- get_summary(UpSet_trimmed_col37, type_d1d2, 'TOC1 and PRR7 and PRR5', 37)
col38_summary_d1d2 <- get_summary(UpSet_trimmed_col38, type_d1d2, 'TOC1 and CCA1', 38)
col39_summary_d1d2 <- get_summary(UpSet_trimmed_col39, type_d1d2, 'LHY and PRR7', 39)
col40_summary_d1d2 <- get_summary(UpSet_trimmed_col40, type_d1d2, 'LHY and PRR7 and PRR5', 40)

# d1d5----
col1_summary_d1d5 <- get_summary(UpSet_trimmed_col1, type_d1d5, 'LHY', 1)
col2_summary_d1d5 <- get_summary(UpSet_trimmed_col2, type_d1d5, 'TOC1', 2)
col3_summary_d1d5 <- get_summary(UpSet_trimmed_col3, type_d1d5, 'LUX', 3)
col4_summary_d1d5 <- get_summary(UpSet_trimmed_col4, type_d1d5, 'LHY and CCA1', 4)
col5_summary_d1d5 <- get_summary(UpSet_trimmed_col5, type_d1d5, 'LUX and ELF3', 5)
col6_summary_d1d5 <- get_summary(UpSet_trimmed_col6, type_d1d5, 'PRR5', 6)
col7_summary_d1d5 <- get_summary(UpSet_trimmed_col7, type_d1d5, 'PRR7', 7)
col8_summary_d1d5 <- get_summary(UpSet_trimmed_col8, type_d1d5, 'LUX and LHY', 8)
col9_summary_d1d5 <- get_summary(UpSet_trimmed_col9, type_d1d5, 'LUX and LHY and CCA1', 9)
col10_summary_d1d5 <- get_summary(UpSet_trimmed_col10, type_d1d5, 'LUX and TOC1', 10)

col11_summary_d1d5 <- get_summary(UpSet_trimmed_col11, type_d1d5, 'CCA1', 11)
col12_summary_d1d5 <- get_summary(UpSet_trimmed_col12, type_d1d5, 'LUX and ELF3 and ELF4', 12)
col13_summary_d1d5 <- get_summary(UpSet_trimmed_col13, type_d1d5, 'LHY and TOC1', 13)
col14_summary_d1d5 <- get_summary(UpSet_trimmed_col14, type_d1d5, 'LHY and TOC1 and CCA1', 14)
col15_summary_d1d5 <- get_summary(UpSet_trimmed_col15, type_d1d5, 'LUX and TOC1 and PRR5', 15)
col16_summary_d1d5 <- get_summary(UpSet_trimmed_col16, type_d1d5, 'LUX and LHY and ELF3', 16)
col17_summary_d1d5 <- get_summary(UpSet_trimmed_col17, type_d1d5, 'LUX and LHY and CCA1 and PRR7', 17)
col18_summary_d1d5 <- get_summary(UpSet_trimmed_col18, type_d1d5, 'ELF3', 18)
col19_summary_d1d5 <- get_summary(UpSet_trimmed_col19, type_d1d5, 'LUX and ELF3 and PRR5 and ELF4', 19)
col20_summary_d1d5 <- get_summary(UpSet_trimmed_col20, type_d1d5, 'LUX and TOC1 and ELF3 and PRR7 and PRR5 and ELF4', 20)

col21_summary_d1d5 <- get_summary(UpSet_trimmed_col21, type_d1d5, 'LUX and LHY and ELF3 and ELF4', 21)
col22_summary_d1d5 <- get_summary(UpSet_trimmed_col22, type_d1d5, 'LUX and LHY and CCA1 and ELF3 and PRR7 and PRR5', 22)
col23_summary_d1d5 <- get_summary(UpSet_trimmed_col23, type_d1d5, 'LUX and LHY and TOC1', 23)
col24_summary_d1d5 <- get_summary(UpSet_trimmed_col24, type_d1d5, 'PRR7 and PRR5', 24)
col25_summary_d1d5 <- get_summary(UpSet_trimmed_col25, type_d1d5, 'TOC1 and PRR7', 25)
col26_summary_d1d5 <- get_summary(UpSet_trimmed_col26, type_d1d5, 'LHY and TOC1 and PRR7', 26)
col27_summary_d1d5 <- get_summary(UpSet_trimmed_col27, type_d1d5, 'LUX and ELF3 and PRR5', 27)
col28_summary_d1d5 <- get_summary(UpSet_trimmed_col28, type_d1d5, 'LUX and CCA1 and ELF3', 28)
col29_summary_d1d5 <- get_summary(UpSet_trimmed_col29, type_d1d5, 'LUX and TOC1 and ELF3', 29)
col30_summary_d1d5 <- get_summary(UpSet_trimmed_col30, type_d1d5, 'LUX and TOC1 and CCA1', 30)

col31_summary_d1d5 <- get_summary(UpSet_trimmed_col31, type_d1d5, 'LUX and LHY and ELF3 and PRR7 and PRR5 and ELF4', 31)
col32_summary_d1d5 <- get_summary(UpSet_trimmed_col32, type_d1d5, 'LUX and LHY and CCA1 and ELF3', 32)
col33_summary_d1d5 <- get_summary(UpSet_trimmed_col33, type_d1d5, 'LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and PRR5', 33)
col34_summary_d1d5 <- get_summary(UpSet_trimmed_col34, type_d1d5, 'LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and PRR5 and ELF4', 34)
col35_summary_d1d5 <- get_summary(UpSet_trimmed_col35, type_d1d5, 'ELF3 and ELF4', 35)
col36_summary_d1d5 <- get_summary(UpSet_trimmed_col36, type_d1d5, 'TOC1 and PRR5', 36)
col37_summary_d1d5 <- get_summary(UpSet_trimmed_col37, type_d1d5, 'TOC1 and PRR7 and PRR5', 37)
col38_summary_d1d5 <- get_summary(UpSet_trimmed_col38, type_d1d5, 'TOC1 and CCA1', 38)
col39_summary_d1d5 <- get_summary(UpSet_trimmed_col39, type_d1d5, 'LHY and PRR7', 39)
col40_summary_d1d5 <- get_summary(UpSet_trimmed_col40, type_d1d5, 'LHY and PRR7 and PRR5', 40)

upset_columns_d1_d2 <- bind_rows(col1_summary_d1d2,
                                 col2_summary_d1d2,
                                 col3_summary_d1d2,
                                 col4_summary_d1d2,
                                 col5_summary_d1d2,
                                 col6_summary_d1d2,
                                 col7_summary_d1d2,
                                 col8_summary_d1d2,
                                 col9_summary_d1d2,
                                 col10_summary_d1d2,
                                 col11_summary_d1d2,
                                 col12_summary_d1d2,
                                 col13_summary_d1d2,
                                 col14_summary_d1d2,
                                 col15_summary_d1d2,
                                 col16_summary_d1d2,
                                 col17_summary_d1d2,
                                 col18_summary_d1d2,
                                 col19_summary_d1d2,
                                 col20_summary_d1d2,
                                 col21_summary_d1d2,
                                 col22_summary_d1d2,
                                 col23_summary_d1d2,
                                 col24_summary_d1d2,
                                 col25_summary_d1d2,
                                 col26_summary_d1d2,
                                 col27_summary_d1d2,
                                 col28_summary_d1d2,
                                 col29_summary_d1d2,
                                 col30_summary_d1d2,
                                 col31_summary_d1d2,
                                 col32_summary_d1d2,
                                 col33_summary_d1d2,
                                 col34_summary_d1d2,
                                 col35_summary_d1d2,
                                 col36_summary_d1d2,
                                 col37_summary_d1d2,
                                 col38_summary_d1d2,
                                 col39_summary_d1d2,
                                 col40_summary_d1d2)

upset_columns_d1_d5 <- bind_rows(col1_summary_d1d5,
                                 col2_summary_d1d5,
                                 col3_summary_d1d5,
                                 col4_summary_d1d5,
                                 col5_summary_d1d5,
                                 col6_summary_d1d5,
                                 col7_summary_d1d5,
                                 col8_summary_d1d5,
                                 col9_summary_d1d5,
                                 col10_summary_d1d5,
                                 col11_summary_d1d5,
                                 col12_summary_d1d5,
                                 col13_summary_d1d5,
                                 col14_summary_d1d5,
                                 col15_summary_d1d5,
                                 col16_summary_d1d5,
                                 col17_summary_d1d5,
                                 col18_summary_d1d5,
                                 col19_summary_d1d5,
                                 col20_summary_d1d5,
                                 col21_summary_d1d5,
                                 col22_summary_d1d5,
                                 col23_summary_d1d5,
                                 col24_summary_d1d5,
                                 col25_summary_d1d5,
                                 col26_summary_d1d5,
                                 col27_summary_d1d5,
                                 col28_summary_d1d5,
                                 col29_summary_d1d5,
                                 col30_summary_d1d5,
                                 col31_summary_d1d5,
                                 col32_summary_d1d5,
                                 col33_summary_d1d5,
                                 col34_summary_d1d5,
                                 col35_summary_d1d5,
                                 col36_summary_d1d5,
                                 col37_summary_d1d5,
                                 col38_summary_d1d5,
                                 col39_summary_d1d5,
                                 col40_summary_d1d5)


# plot d1d2----
plot_upset_columns_d1_d2 <- upset_columns_d1_d2 %>% 
  mutate(type = factor(type_d1d2, levels = c('gain_high_d1_d2', 'gain_medium_d1_d2', 'other_d1_d2', 'lose_medium_d1_d2', 'lose_high_d1_d2')),
         clock = fct_reorder(clock, column)) %>% 
  ggplot(aes(x = clock, y = count, fill = type)) +
  scale_y_continuous(limits = c(0, 80), breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_fill_manual(values = c('#31a354', '#74c476', '#cccccc', '#fb6a4a', '#de2d26'), labels = c('gain - high', 'gain - medium', 'other', 'lose - medium', 'lose - high')) +
  geom_bar(position = 'stack', stat = 'identity') +
  ggpubr::theme_pubr() +
  theme(legend.position = "right") +
  theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1)) +
  labs(fill = 'Amplitude', y = 'number', x = '') +
  ggtitle("Amplitude patterns for clock targets",
          subtitle = "Compare day 1 (20C steady state) with day 2 (to 4C transient-cooling)") 

plot_upset_columns_d1_d2

ggsave('./03_plots/plot_clock_d1d2.png', dpi = 300, height = 6, width = 6, units = 'in')

# plot d1d5----
plot_upset_columns_d1_d5 <- upset_columns_d1_d5 %>% 
  mutate(type = factor(type_d1d5, levels = c('gain_high_d1_d5', 'gain_medium_d1_d5', 'other_d1_d5', 'lose_medium_d1_d5', 'lose_high_d1_d5')),
         clock = fct_reorder(clock, column)) %>% 
  ggplot(aes(x = clock, y = count, fill = type)) +
  scale_y_continuous(limits = c(0, 80), breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_fill_manual(values = c('#31a354', '#74c476', '#cccccc', '#fb6a4a', '#de2d26'), labels = c('gain - high', 'gain - medium', 'other', 'lose - medium', 'lose - high')) +
  geom_bar(position = 'stack', stat = 'identity') +
  ggpubr::theme_pubr() +
  theme(legend.position = "right") +
  theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1)) +
  labs(fill = 'Amplitude', y = 'number', x = '') +
  ggtitle("Amplitude patterns for clock targets",
          subtitle = "Compare day 1 (20C steady state) with day 2 (to 4C transient-cooling)") 

plot_upset_columns_d1_d5

ggsave('./03_plots/plot_clock_d1d2.png', dpi = 300, height = 6, width = 6, units = 'in')

sets_trimmed_clock <- sets_trimmed %>% 
  mutate(clock = case_when(LHY == 1 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'LHY',
                           LHY == 0 & CCA1 == 1 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'CCA1',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 0 & PRR7 == 0 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'TOC1',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 1 & PRR7 == 0 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'PRR5',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 1 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'PRR7',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'LUX',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 0 & ELF3 == 1 & ELF4 == 0 ~ 'ELF3',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 0 & ELF3 == 0 & ELF4 == 1 ~ 'ELF4',
                           LHY == 1 & CCA1 == 1 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'LHY and CCA1',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'LUX and ELF3',
                           LHY == 1 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'LUX and LHY',
                           LHY == 1 & CCA1 == 1 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'LUX and LHY and CCA1',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'LUX and TOC1',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'LUX and ELF3 and ELF4',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 1 & PRR7 == 0 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'LUX and TOC1 and PRR5',
                           LHY == 1 & CCA1 == 0 & TOC1 == 1 & PRR5 == 0 & PRR7 == 0 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'LHY and TOC1',
                           LHY == 1 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'LUX and LHY and ELF3',
                           LHY == 1 & CCA1 == 1 & TOC1 == 1 & PRR5 == 0 & PRR7 == 0 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'LHY and TOC1 and CCA1',
                           LHY == 1 & CCA1 == 1 & TOC1 == 0 & PRR5 == 0 & PRR7 == 1 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'LUX and LHY and CCA1 and PRR7',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 1 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'LUX and TOC1 and ELF3 and PRR7 and PRR5 and ELF4',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 1 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'LUX and ELF3 and PRR5 and ELF4',
                           LHY == 1 & CCA1 == 0 & TOC1 == 1 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'LUX and LHY and TOC1',
                           LHY == 1 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'LUX and LHY and ELF3 and ELF4',
                           LHY == 1 & CCA1 == 1 & TOC1 == 0 & PRR5 == 1 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'LUX and LHY and CCA1 and ELF3 and PRR7 and PRR5',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 0 & PRR7 == 1 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'TOC1 and PRR7',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'LUX and TOC1 and ELF3',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 1 & PRR7 == 1 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'PRR7 and PRR5',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 1 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'LUX and ELF3 and PRR5',
                           LHY == 1 & CCA1 == 0 & TOC1 == 1 & PRR5 == 0 & PRR7 == 1 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'LHY and TOC1 and PRR7',
                           LHY == 1 & CCA1 == 0 & TOC1 == 0 & PRR5 == 1 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'LUX and LHY and ELF3 and PRR7 and PRR5 and ELF4',
                           LHY == 1 & CCA1 == 1 & TOC1 == 1 & PRR5 == 1 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and PRR5 and ELF4',
                           LHY == 1 & CCA1 == 1 & TOC1 == 1 & PRR5 == 1 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and PRR5',
                           LHY == 1 & CCA1 == 1 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'LUX and LHY and CCA1 and ELF3',
                           LHY == 0 & CCA1 == 1 & TOC1 == 1 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'LUX and TOC1 and CCA1',
                           LHY == 0 & CCA1 == 1 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'LUX and CCA1 and ELF3',
                           LHY == 1 & CCA1 == 1 & TOC1 == 0 & PRR5 == 0 & PRR7 == 1 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'LHY and CCA1 and PRR7',
                           LHY == 1 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 1 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'LHY and PRR7',
                           LHY == 1 & CCA1 == 1 & TOC1 == 1 & PRR5 == 1 & PRR7 == 1 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'LHY and TOC1 and CCA1 and PRR7 and PRR5',
                           LHY == 1 & CCA1 == 0 & TOC1 == 0 & PRR5 == 1 & PRR7 == 1 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'LHY and PRR7 and PRR5',
                           LHY == 1 & CCA1 == 1 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 0 & ELF3 == 1 & ELF4 == 0 ~ 'LHY and CCA1 and ELF3',
                           LHY == 1 & CCA1 == 1 & TOC1 == 1 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'LUX and LHY and TOC1 and CCA1',
                           LHY == 1 & CCA1 == 0 & TOC1 == 1 & PRR5 == 1 & PRR7 == 1 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'LHY and TOC1 and PRR7 and PRR5',
                           LHY == 1 & CCA1 == 1 & TOC1 == 0 & PRR5 == 1 & PRR7 == 0 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'LHY and CCA1 and PRR5',
                           LHY == 1 & CCA1 == 0 & TOC1 == 1 & PRR5 == 1 & PRR7 == 0 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'LUX and LHY and TOC1 and PRR5',
                           LHY == 1 & CCA1 == 0 & TOC1 == 1 & PRR5 == 1 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'LUX and LHY and TOC1 and ELF3 and PRR7 and PRR5 and ELF4',
                           LHY == 1 & CCA1 == 1 & TOC1 == 1 & PRR5 == 0 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'LUX and LHY and TOC1 and CCA1 and ELF3 and PRR7 and ELF4',
                           LHY == 0 & CCA1 == 1 & TOC1 == 1 & PRR5 == 1 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'CCA1 and TOC1 and PRR5 and PRR7 and LUX and ELF3',
                           LHY == 0 & CCA1 == 1 & TOC1 == 1 & PRR5 == 0 & PRR7 == 0 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'TOC1 and CCA1',
                           LHY == 0 & CCA1 == 1 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'CCA1 and LUX and ELF3 and ELF4',
                           LHY == 0 & CCA1 == 1 & TOC1 == 1 & PRR5 == 1 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'CCA1 and TOC1 and PRR5 and LUX and ELF3 and ELF4',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 1 & PRR7 == 0 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'TOC1 and PRR5',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 0 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'TOC1 and PRR7 and LUX and ELF3',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 1 & PRR7 == 0 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'TOC1 and PRR5 and LUX and ELF3',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 1 & PRR7 == 1 & LUX == 0 & ELF3 == 0 & ELF4 == 0 ~ 'TOC1 and PRR7 and PRR5',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 0 & PRR7 == 1 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'TOC1 and PRR7 and LUX',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 1 & PRR7 == 1 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'TOC1 and PRR5 and PRR7 and LUX',
                           LHY == 0 & CCA1 == 0 & TOC1 == 1 & PRR5 == 0 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'TOC1 and PRR7 and LUX and ELF3 and ELF4',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 1 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'PRR5 and PRR7 and LUX and ELF3',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 1 & PRR7 == 1 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'PRR5 and PRR7 and LUX',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 0 ~ 'PRR7 and LUX and ELF3',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 1 & LUX == 1 & ELF3 == 1 & ELF4 == 1 ~ 'PRR7 and LUX and ELF3 and ELF4',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 1 & LUX == 1 & ELF3 == 0 & ELF4 == 0 ~ 'PRR7 and LUX',
                           LHY == 0 & CCA1 == 0 & TOC1 == 0 & PRR5 == 0 & PRR7 == 0 & LUX == 0 & ELF3 == 1 & ELF4 == 1 ~ 'ELF3 and ELF4',
                           TRUE ~ 'NA'
                           ))

class(sets_trimmed)
table(sets_trimmed_clock$clock)
glimpse(UpSet_trimmed_col1_edit)

UpSet_trimmed_col1_edit <- UpSet_trimmed_col1 %>% select(-1)
UpSet_trimmed_col2_edit <- UpSet_trimmed_col2 %>% select(-1)

UpSet_trimmed_col1_col2_edit <- bind_rows(UpSet_trimmed_col1_edit, UpSet_trimmed_col2_edit)



sets_trimmed_clock_types <- sets_trimmed_clock %>% left_join(UpSet_trimmed_col1_edit, join_by(clock))
sets_trimmed_clock_types <- sets_trimmed_clock %>% merge(UpSet_trimmed_col1_edit, by = 'clock')

test_join <- left_join(sets_trimmed_clock %>% group_by(clock) %>% mutate(id = row_number()),
                       UpSet_trimmed_col1_edit %>% group_by(clock) %>% mutate(id = row_number()),
                       UpSet_trimmed_col2_edit %>% group_by(clock) %>% mutate(id = row_number()),
                       by = c("clock", "id"))

test_join3 <- sets_trimmed_clock %>% group_by(clock) %>% mutate(id = row_number()) %>% 
  left_join(UpSet_trimmed_col1_col2_edit %>% group_by(clock) %>% mutate(id = row_number()))

test_join3 %<>% left_join(UpSet_trimmed_col2_edit %>% group_by(clock) %>% mutate(id = row_number()))


#purrr::reduce(list(UpSet_trimmed_col1_edit, UpSet_trimmed_col2_edit), dplyr::left_join, by = 'Flag')

# test_join2 <- left_join(sets_trimmed_clock %>% group_by(clock) %>% mutate(id = row_number()) %>% 
#                           ,
#                        UpSet_trimmed_col2_edit %>% group_by(clock) %>% mutate(id = row_number()), 
#                        by = c("clock", "id"))

?left_join
clock_components = colnames(sets_trimmed)[1:8]
clock_components

sets_trimmed[clock_components] = sets_trimmed[clock_components] == 1
upset(sets_trimmed, 
      clock_components, 
      name='genre',
      width_ratio=0.1, min_size = 1)

data <- data.frame(
  entry = paste0("Entry.", 1:10),
  "A" = c(0,0,0,0,1,0,1,1,0,0),
  "B" = c(1,0,0,0,1,1,1,1,1,0),
  "C" = c(1,1,1,1,0,0,1,0,1,1)
)

data
data_with_intersection <- data %>%
  unite(col = "intersection", -c("entry"), sep = "")

data_with_intersection %>%
  dplyr::group_by(intersection) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::arrange(desc(n))

data_with_intersection %>%
  dplyr::group_by(intersection) %>%
  dplyr::summarise(list = list(entry)) %>%
  dplyr::mutate(list = setNames(list, intersection)) %>%
  pull(list)
