# a meta analysis of Arabidopsis clock ChIP data and overlap with cold TF network

library(tidyverse)
library(janitor)
library(MetaCycle)
devtools::install_github("vqv/ggbiplot")
library(ggbiplot)
library(ggrepel)

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
day1_output <- read_csv('./cluster_days_output_day1/meta2d_clusters_aggregated_day1.csv') %>% 
  select(CycID, meta2d_pvalue, meta2d_period, meta2d_phase, meta2d_AMP) %>% 
  rename(cluster = CycID, d1_meta2d_pvalue = meta2d_pvalue, d1_meta2d_period = meta2d_period, d1_meta2d_phase = meta2d_phase, d1_meta2d_AMP = meta2d_AMP)

day2_output <- read_csv('./cluster_days_output_day2/meta2d_clusters_aggregated_day2.csv') %>% 
  select(CycID, meta2d_pvalue, meta2d_period, meta2d_phase, meta2d_AMP) %>% 
  rename(cluster = CycID, d2_meta2d_pvalue = meta2d_pvalue, d2_meta2d_period = meta2d_period, d2_meta2d_phase = meta2d_phase, d2_meta2d_AMP = meta2d_AMP)

day5_output <- read_csv('./cluster_days_output_day5/meta2d_clusters_aggregated_day5.csv') %>% 
  select(CycID, meta2d_pvalue, meta2d_period, meta2d_phase, meta2d_AMP) %>% 
  rename(cluster = CycID, d5_meta2d_pvalue = meta2d_pvalue, d5_meta2d_period = meta2d_period, d5_meta2d_phase = meta2d_phase, d5_meta2d_AMP = meta2d_AMP)

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

write_csv(day1_vs_day2_150, 'day1_vs_day2_150.csv')

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

write_csv(day1_vs_day5_150, 'day1_vs_day5_150.csv')

# Scatter plot of the MetaCycle outputs for Amplitude coloured by the amp_flag grouping
plot_day1_vs_day2_150_AMP <- day1_vs_day2_150 %>% 
  #filter(d1_meta2d_AMP > 0.5 & d2_meta2d_AMP > 0.5) %>% 
  filter(pval_flag == 'd1_nr_d2_r' | pval_flag == 'd1_r_d2_nr' | pval_flag == 'd1_r_d2_r') %>%
  ggplot(aes(x = d1_meta2d_AMP, y = d2_meta2d_AMP, colour = amp_flag, shape = pval_flag)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  scale_x_continuous(limits = c(0, 1.8), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  geom_point(size = 2.5) +
  ggpubr::theme_pubr() +
  #ggpubr::theme_pubclean() +
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

plot_day1_vs_day5_150_AMP <- day1_vs_day5_150 %>%
  #filter(d1_meta2d_AMP > 0.5 & d2_meta2d_AMP > 0.5) %>% 
  filter(pval_flag == 'd1_nr_d5_r' | pval_flag == 'd1_r_d5_nr' | pval_flag == 'd1_r_d5_r') %>%
  ggplot(aes(x = d1_meta2d_AMP, y = d5_meta2d_AMP, colour = amp_flag, shape = pval_flag)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  scale_x_continuous(limits = c(0, 1.8), breaks = c(0, 0.4, 0.8, 1.2, 1.6)) +
  geom_point(size = 2.5) +
  ggpubr::theme_pubr() +
  #ggpubr::theme_pubclean() +
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


## An analysis of established and published Arabidopsis clock ChIP targets in TF network

# read the TF network cluster
# The TF network is 7302 genes over 75 clusters (clusters 0-74)

#setwd('/Users/Allan/Documents/Post crash 2016/Home laptop June 2013/Research/Splicing II/Network/Clock ChIP datasets/Clock_ChIP_datasets/')

TF <- read_csv("TF Network Cluster Nov2018.csv") %>%
  pivot_longer(cols = starts_with ('cluster'),
               names_to = "cluster", 
               values_to = "gene_ID",
               values_drop_na = TRUE) %>% 
  filter(grepl('AT', gene_ID)) %>%
  mutate(cluster = str_sub(cluster, 9, -1)) %>% 
  write_csv("TF Network Cluster Nov2018 pivot longer.csv")

# LHY dataset
# read in the Adams LHY paper dataset and skip first 2 lines
# The LHY paper is Adams et al. (2018) New Phytologist 220(3); 897
# supplemental data set (Table S2) 
adams <- read_csv("nph15415-sup-0002-tables2.csv",skip=2) %>% 
  dplyr::select(gene_ID = 1) %>%
  filter(!is.na(gene_ID)) %>%
  distinct(gene_ID) 
#%>% 
write_csv("LHY_targets.csv")

# CCA1 Nagel et al. dataset
# read in the Nagel CCA1 paper dataset and skip first 2 lines
# The CCA1 paper is Nagel et al. (2015) PNAS 112(34); E4802
# supplemental data set (Table S1) 
nagel <- read_csv("pnas.1513609112.sd01.csv",skip=2) %>% 
  dplyr::select(gene_ID = 10) %>% 
  mutate(gene_ID = str_sub(gene_ID, end = 9)) %>% 
  distinct(gene_ID) 
#%>% 
write_csv("./Clock ChIP targets/CCA1_nagel_targets.csv")

# CCA1 Kamioka et al. dataset
# read in the Kamioka CCA1 paper dataset and skip first 2 lines
# The CCA1 paper is Kamioka et al. (2016) Plant Cell 28(3); 696
# supplemental data set (Table S1C)
kamioka <- read_csv("TPC2015-00737-RAR3_Supplemental_Data_Set_1C.csv",skip=3) %>% 
  dplyr::select(gene_ID = 10) %>% 
  distinct(gene_ID) 
#%>% 
write_csv("./Clock ChIP targets/CCA1_kamioka_targets.csv")

# merge the nagel and kamioka CCA1 datasets
# use inner_join from dplyr
# 249 obs.
kamioka_nagel_merge <- inner_join(nagel, kamioka, by = "gene_ID")

# TOC1 dataset
# read in the Huang TOC1 paper dataset
# The TOC1 paper is Huang et al. (2012) Science 336:75
# supplemental data set (Table S1)
huang <- read_csv("Huang TOC1 CHiP TableS1.csv") %>% 
  dplyr::select(gene_ID = 14) %>% 
  distinct(gene_ID) 
#%>% 
write_csv("./Clock ChIP targets/TOC1_huang_targets.csv")

# PRR5 dataset
# read in the Nakamichi PRR5 paper dataset
# The PRR5 paper is Nakamichi et al. (2012) PNAS 109:17123
# supplemental data set (Table S3) 
nakamichi <- read_csv("Dataset S3 Nakamichi et al PRR5 binding targets PNAS 2012.csv", skip=2) %>% 
  dplyr::select(gene_ID = 3) %>% 
  distinct(gene_ID) 
#%>% 
write_csv("./Clock ChIP targets/PRR5_nakamichi_targets.csv")

# PRR7 dataset
# read in the Liu PRR7 paper dataset
# The PRR7 paper is Liu et al. (2013) The Plant Journal 76:101
# supplemental data set (Table S1)
liu <- read_csv("Dataset S1 Liu et al PRR7 edit.csv") %>% 
  dplyr::select(gene_ID = 17) %>% 
  distinct(gene_ID) 
#%>% 
write_csv("./Clock ChIP targets/PRR7_liu_targets.csv")

# LUX dataset
# read in the Ezer EC paper for the LUX dataset (LUX_17 tab)
# The Evening Complex (EC) paper is Ezer et al. (2017) Nature Plants 3: article 17087
# supplemental data set (LUX_17 tab of Table S6) 
ezer_LUX <- read_csv("Ezer et al nplants Suppl Table S6.csv") %>% 
  dplyr::select(gene_ID = 1) %>% 
  distinct(gene_ID) 
#%>% 
write_csv("./Clock ChIP targets/LUX_ezer_targets.csv")

# ELF3 dataset
# read in the Ezer EC paper for the ELF3 dataset (ELF3_22 tab)
# The Evening Complex (EC) paper is Ezer et al. (2017) Nature Plants 3: article 17087
# supplemental data set (ELF3_22 tab of Table S6) 
ezer_ELF3 <- read_csv("ELF3_22 Ezer Table S6.csv") %>% 
  dplyr::select(gene_ID = 1) %>% 
  distinct(gene_ID) 
#%>% 
write_csv("./Clock ChIP targets/ELF3_ezer_targets.csv")

# ELF4 dataset
# read in the Ezer EC paper for the ELF4 dataset (ELF4_22 tab)
# The Evening Complex (EC) paper is Ezer et al. (2017) Nature Plants 3: article 17087
# supplemental data set (ELF4_22 tab of Table S6) 
ezer_ELF4 <- read_csv("ELF4_22 Ezer Table S6.csv") %>% 
  dplyr::select(gene_ID = 1) %>% 
  distinct(gene_ID) 
#%>% 
write_csv("./Clock ChIP targets/ELF4_ezer_targets.csv")


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

table(LHY_bind_d1d2$type)

# CCA1-Nagel amp_gain_high (69 obs) + CCA1-Nagel amp_gain_medium (151 obs) + CCA1-Nagel amp_lose_high (97 obs) + CCA1-Nagel amp_lose_medium (33 obs) + CCA1-Nagel amp_other (183 obs): total equals 533 obs
CCA1_nagel_bind_d1d2 <- bind_rows(CCA1_nagel_gain_high_d1d2, CCA1_nagel_gain_medium_d1d2, CCA1_nagel_lose_high_d1d2, CCA1_nagel_lose_medium_d1d2, CCA1_nagel_other_d1d2)

names(CCA1_nagel_bind_d1d2)[3] <- "CCA1 Nagel"

table(CCA1_nagel_bind_d1d2$type)

# CCA1-Kamioka amp_gain_high (37 obs) + CCA1-Kamioka amp_gain_medium (59 obs) + CCA1-Kamioka amp_lose_high (24 obs) + CCA1-Kamioka amp_lose_medium (14 obs) + CCA1-Kamioka amp_other (68 obs): total equals 202 obs
CCA1_kamioka_bind_d1d2 <- bind_rows(CCA1_kamioka_gain_high_d1d2, CCA1_kamioka_gain_medium_d1d2, CCA1_kamioka_lose_high_d1d2, CCA1_kamioka_lose_medium_d1d2, CCA1_kamioka_other_d1d2)

names(CCA1_kamioka_bind_d1d2)[3] <- "CCA1 Kamioka"

table(CCA1_kamioka_bind_d1d2$type)

# CCA1-Nagel-Kamioka amp_gain_high (22 obs) + CCA1-Nagel-Kamioka amp_gain_medium (44 obs) + CCA1-Nagel-Kamioka amp_lose_high (10 obs) + CCA1-Nagel-Kamioka amp_lose_medium (7 obs) + CCA1-Nagel-Kamioka amp_other (50 obs): total equals 133 obs
CCA1_nagel_kamioka_bind_d1d2 <- bind_rows(CCA1_nagel_kamioka_gain_high_d1d2, CCA1_nagel_kamioka_gain_medium_d1d2, CCA1_nagel_kamioka_lose_high_d1d2, CCA1_nagel_kamioka_lose_medium_d1d2, CCA1_nagel_kamioka_other_d1d2)

names(CCA1_nagel_kamioka_bind_d1d2)[3] <- "CCA1 Nagel-Kamioka"

table(CCA1_nagel_kamioka_bind_d1d2$type)

# TOC1 amp_gain_high (79 obs) + TOC1 amp_gain_medium (73 obs) + TOC1 amp_lose_high (45 obs) + TOC1 amp_lose_medium (30 obs) + TOC1 amp_other (65 obs): total equals 292 obs
TOC1_bind_d1d2 <- bind_rows(TOC1_gain_high_d1d2, TOC1_gain_medium_d1d2, TOC1_lose_high_d1d2, TOC1_lose_medium_d1d2, TOC1_other_d1d2)

names(TOC1_bind_d1d2)[3] <- "TOC1"

table(TOC1_bind_d1d2$type)

# PRR5 amp_gain_high (30 obs) + PRR5 amp_gain_medium (14 obs) + PRR5 amp_lose_high (2 obs) + PRR5 amp_lose_medium (4 obs) + PRR5 amp_other (6 obs): total equals 56 obs
PRR5_bind_d1d2 <- bind_rows(PRR5_gain_high_d1d2, PRR5_gain_medium_d1d2, PRR5_lose_high_d1d2, PRR5_lose_medium_d1d2, PRR5_other_d1d2)

names(PRR5_bind_d1d2)[3] <- "PRR5"

table(PRR5_bind_d1d2$type)

# PRR7 amp_gain_high (29 obs) + PRR7 amp_gain_medium (10 obs) + PRR7 amp_lose_high (11 obs) + PRR7 amp_lose_medium (2 obs) + PRR7 amp_other (14 obs): total equals 66 obs
PRR7_bind_d1d2 <- bind_rows(PRR7_gain_high_d1d2, PRR7_gain_medium_d1d2, PRR7_lose_high_d1d2, PRR7_lose_medium_d1d2, PRR7_other_d1d2)

names(PRR7_bind_d1d2)[3] <- "PRR7"

table(PRR7_bind_d1d2$type)

# LUX amp_gain_high (58 obs) + LUX amp_gain_medium (86 obs) + LUX amp_lose_high (103 obs) + LUX amp_lose_medium (42 obs) + LUX amp_other (103 obs): total equals 392 obs
LUX_bind_d1d2 <- bind_rows(LUX_gain_high_d1d2, LUX_gain_medium_d1d2, LUX_lose_high_d1d2, LUX_lose_medium_d1d2, LUX_other_d1d2)

names(LUX_bind_d1d2)[3] <- "LUX"

table(LUX_bind_d1d2$type)

# ELF3 amp_gain_high (19 obs) + ELF3 amp_gain_medium (19 obs) + ELF3 amp_lose_high (44 obs) + ELF3 amp_lose_medium (11 obs) + ELF3 amp_other (27 obs): total equals 120 obs
ELF3_bind_d1d2 <- bind_rows(ELF3_gain_high_d1d2, ELF3_gain_medium_d1d2, ELF3_lose_high_d1d2, ELF3_lose_medium_d1d2, ELF3_other_d1d2)

names(ELF3_bind_d1d2)[3] <- "ELF3"

table(ELF3_bind_d1d2$type)

# ELF4 amp_gain_high (8 obs) + ELF4 amp_gain_medium (3 obs) + ELF4 amp_lose_high (12 obs) + ELF4 amp_lose_medium (2 obs) + ELF4 amp_other (5 obs): total equals 30 obs
ELF4_bind_d1d2 <- bind_rows(ELF4_gain_high_d1d2, ELF4_gain_medium_d1d2, ELF4_lose_high_d1d2, ELF4_lose_medium_d1d2, ELF4_other_d1d2)

names(ELF4_bind_d1d2)[3] <- "ELF4"

table(ELF4_bind_d1d2$type)

# * d1d5----
# append all the LHY targets:
# LHY amp_gain_high (0 obs) + LHY amp_gain_medium (16 obs) + LHY amp_lose_high (87 obs) + LHY amp_lose_medium (53 obs) + LHY amp_other (144 obs): total equals 300 obs
LHY_bind_d1d5 <- bind_rows(LHY_gain_high_d1d5, LHY_gain_medium_d1d5, LHY_lose_high_d1d5, LHY_lose_medium_d1d5, LHY_other_d1d5)

names(LHY_bind_d1d5)[3] <- "LHY"

table(LHY_bind_d1d5$type)

# CCA1-Nagel amp_gain_high (0 obs) + CCA1-Nagel amp_gain_medium (26 obs) + CCA1-Nagel amp_lose_high (166 obs) + CCA1-Nagel amp_lose_medium (82 obs) + CCA1-Nagel amp_other (210 obs): total equals 484 obs
CCA1_nagel_bind_d1d5 <- bind_rows(CCA1_nagel_gain_high_d1d5, CCA1_nagel_gain_medium_d1d5, CCA1_nagel_lose_high_d1d5, CCA1_nagel_lose_medium_d1d5, CCA1_nagel_other_d1d5)

names(CCA1_nagel_bind_d1d5)[3] <- "CCA1 Nagel"

table(CCA1_nagel_bind_d1d5$type)

# CCA1-Kamioka amp_gain_high (1 obs) + CCA1-Kamioka amp_gain_medium (6 obs) + CCA1-Kamioka amp_lose_high (53 obs) + CCA1-Kamioka amp_lose_medium (31 obs) + CCA1-Kamioka amp_other (88 obs): total equals 179 obs
CCA1_kamioka_bind_d1d5 <- bind_rows(CCA1_kamioka_gain_high_d1d5, CCA1_kamioka_gain_medium_d1d5, CCA1_kamioka_lose_high_d1d5, CCA1_kamioka_lose_medium_d1d5, CCA1_kamioka_other_d1d5)

names(CCA1_kamioka_bind_d1d5)[3] <- "CCA1 Kamioka"

table(CCA1_kamioka_bind_d1d5$type)

# CCA1-Nagel-Kamioka amp_gain_high (0 obs) + CCA1-Nagel-Kamioka amp_gain_medium (4 obs) + CCA1-Nagel-Kamioka amp_lose_high (21 obs) + CCA1-Nagel-Kamioka amp_lose_medium (28 obs) + CCA1-Nagel-Kamioka amp_other (64 obs): total equals 117 obs
CCA1_nagel_kamioka_bind_d1d5 <- bind_rows(CCA1_nagel_kamioka_gain_high_d1d5, CCA1_nagel_kamioka_gain_medium_d1d5, CCA1_nagel_kamioka_lose_high_d1d5, CCA1_nagel_kamioka_lose_medium_d1d5, CCA1_nagel_kamioka_other_d1d5)

names(CCA1_nagel_kamioka_bind_d1d5)[3] <- "CCA1 Nagel-Kamioka"

table(CCA1_nagel_kamioka_bind_d1d5$type)

# TOC1 amp_gain_high (0 obs) + TOC1 amp_gain_medium (13 obs) + TOC1 amp_lose_high (95 obs) + TOC1 amp_lose_medium (30 obs) + TOC1 amp_other (110 obs): total equals 248 obs
TOC1_bind_d1d5 <- bind_rows(TOC1_gain_high_d1d5, TOC1_gain_medium_d1d5, TOC1_lose_high_d1d5, TOC1_lose_medium_d1d5, TOC1_other_d1d5)

names(TOC1_bind_d1d5)[3] <- "TOC1"

table(TOC1_bind_d1d5$type)

# PRR5 amp_gain_high (0 obs) + PRR5 amp_gain_medium (1 obs) + PRR5 amp_lose_high (6 obs) + PRR5 amp_lose_medium (5 obs) + PRR5 amp_other (31 obs): total equals 43 obs
PRR5_bind_d1d5 <- bind_rows(PRR5_gain_high_d1d5, PRR5_gain_medium_d1d5, PRR5_lose_high_d1d5, PRR5_lose_medium_d1d5, PRR5_other_d1d5)

names(PRR5_bind_d1d5)[3] <- "PRR5"

table(PRR5_bind_d1d5$type)

# PRR7 amp_gain_high (0 obs) + PRR7 amp_gain_medium (3 obs) + PRR7 amp_lose_high (18 obs) + PRR7 amp_lose_medium (6 obs) + PRR7 amp_other (33 obs): total equals 60 obs
PRR7_bind_d1d5 <- bind_rows(PRR7_gain_high_d1d5, PRR7_gain_medium_d1d5, PRR7_lose_high_d1d5, PRR7_lose_medium_d1d5, PRR7_other_d1d5)

names(PRR7_bind_d1d5)[3] <- "PRR7"

table(PRR7_bind_d1d5$type)

# LUX amp_gain_high (3 obs) + LUX amp_gain_medium (21 obs) + LUX amp_lose_high (181 obs) + LUX amp_lose_medium (43 obs) + LUX amp_other (116 obs): total equals 364 obs
LUX_bind_d1d5 <- bind_rows(LUX_gain_high_d1d5, LUX_gain_medium_d1d5, LUX_lose_high_d1d5, LUX_lose_medium_d1d5, LUX_other_d1d5)

names(LUX_bind_d1d5)[3] <- "LUX"

table(LUX_bind_d1d5$type)

# ELF3 amp_gain_high (0 obs) + ELF3 amp_gain_medium (7 obs) + ELF3 amp_lose_high (59 obs) + ELF3 amp_lose_medium (16 obs) + ELF3 amp_other (31 obs): total equals 113 obs
ELF3_bind_d1d5 <- bind_rows(ELF3_gain_high_d1d5, ELF3_gain_medium_d1d5, ELF3_lose_high_d1d5, ELF3_lose_medium_d1d5, ELF3_other_d1d5)

names(ELF3_bind_d1d5)[3] <- "ELF3"

table(ELF3_bind_d1d5$type)

# ELF4 amp_gain_high (0 obs) + ELF4 amp_gain_medium (0 obs) + ELF4 amp_lose_high (14 obs) + ELF4 amp_lose_medium (3 obs) + ELF4 amp_other (10 obs): total equals 27 obs
ELF4_bind_d1d5 <- bind_rows(ELF4_gain_high_d1d5, ELF4_gain_medium_d1d5, ELF4_lose_high_d1d5, ELF4_lose_medium_d1d5, ELF4_other_d1d5)

names(ELF4_bind_d1d5)[3] <- "ELF4"

table(ELF4_bind_d1d5$type)


summarise_targets <- function(df, col_str, clock_id){
  
  summary <- df %>% 
    group_by({{col_str}}) %>% 
    summarise(n=n()) %>% 
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
  labs(fill = 'Amplitude')

plot_clock_d1d2

plot_clock_d1d5 <- clock_d1_d5 %>% 
  mutate(type = factor(type, levels = c('gain_high_d1_d5', 'gain_medium_d1_d5', 'other_d1_d5', 'lose_medium_d1_d5', 'lose_high_d1_d5')),
         clock = factor(clock, levels = c('CCA1', 'LHY', 'TOC1', 'PRR5', 'PRR7', 'LUX', 'ELF3', 'ELF4'))) %>% 
  ggplot(aes(x = clock, y = freq, fill = type)) +
  geom_bar(position = 'stack', stat = 'identity')


plot_clock_d1d5

TF_clock_summary <-
  bind_rows(TF_LHY_amp_append_summary, TF_CCA1_amp_append_summary, TF_TOC1_amp_append_summary, TF_PRR5_amp_append_summary,
            TF_PRR7_amp_append_summary, TF_LUX_amp_append_summary, TF_ELF3_amp_append_summary, TF_ELF4_amp_append_summary) %>% 
  mutate(order_type = case_when(type == 'gain_high_d1_d2' ~ 1,
                                type == 'gain_medium_d1_d2' ~ 2,
                                type == 'other_d1_d2' ~ 3,
                                type == 'lose_medium_d1_d2' ~ 4,
                                TRUE ~ 5),
         order_clock = case_when(clock == 'CCA1' ~ 1,
                                 clock == 'LHY' ~ 2,
                                 clock == 'PRR5' ~ 3,
                                 clock == 'PRR7' ~ 4,
                                 clock == 'TOC1' ~ 5,
                                 clock == 'LUX' ~ 6,
                                 clock == 'ELF3' ~ 7,
                                 TRUE ~ 8),
         type = factor(type),
         clock = factor(clock),
         type = fct_reorder(type, order_type, min),
         clock = fct_reorder(clock, order_clock, min)) %>% 
  ggplot(aes(fill = type, y = freq, x = clock)) +
  geom_col() 
#+geom_bar(position="stack", stat="identity")

TF_clock_summary

glimpse(TF_clock_summary)

#%>%
summarise(mean = mean)
group_by(clock, type)
mutate(type = factor(type),
       type = fct_reorder(type))

TF_clock_summary %>% 
  group_by(type) %>% 
  summarise(max = max(freq)) %>% 
  arrange(desc(max))

clock_rel_abund <-TF_clock_summary %>% 
  group_by(clock) %>% 
  summarize(max = max(freq))

TF_clock_summary <- as.factor(TF_clock_summary$clock, levels = c('LHY', 'CCA1', 'TOC1', 'PRR5', 'PRR7', 'LUX', 'ELF3', 'ELF4'))

TF_clock_summary_bar_plot <-
  ggplot(TF_clock_summary, aes(fill=factor(type, levels=c("gain_high_d1_d2", "gain_medium_d1_d2", "lose_high_d1_d2", "lose_medium_d1_d2", "other_d1_d2")), y=freq,
                               x=factor(clock, levels = c('LHY', 'CCA1', 'TOC1', 'PRR5', 'PRR7', 'LUX', 'ELF3', 'ELF4')))) +
  geom_bar(position="stack", stat="identity")

TF_clock_summary_bar_plot2 <-
  ggplot(TF_clock_summary, aes(fill = type, y = freq, x = clock)) +
  geom_bar(position="stack", stat="identity")

TF_clock_summary_bar_plot2
