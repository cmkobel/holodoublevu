#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)
source("workflow/scripts/utils.R")

# --- Inputs

proteome_intensities_raw <- read_rds(snakemake@input["proteome_intensities"] %>% as.character())
# proteome_intensities_raw <- read_rds("/glittertind/home/carl/PhD/26_proteomics_analysis/proteomic_integration/results/00_uniform_data/proteomics_long_v2.rds")

# --- Clean

# Check that the sample names are uniform.
proteome_intensities_raw %>%
    slice_sample(n = 100) %>%
    view()

# Rename L samples
proteome_intensities <- proteome_intensities_raw %>%
    select(-source) %>%
    supacow_separate_sample(keep_debug_columns = F) %>%
    mutate(timepoint = case_when(
        source == "L" ~ "T6",
        TRUE ~ timepoint
    )) %>%
    mutate(solution = case_when(
        source == "L" ~ "R",
        TRUE ~ solution
    )) %>%
    mutate(collection = case_when(
        timepoint == "T6" ~ "slaughter",
        TRUE ~ "tube"
    ))


proteome_intensities %>%
    random_head()


# --- Outlier removal


threshold_n_D_slaughter <- 5000
threshold_n_D_tube <- 0
threshold_n_L_slaughter <- 0
threshold_n_W_slaughter = 1000

filter1 <- proteome_intensities %>%
    group_by(source, animal, collection) %>%
    summarize(n = n(), sum = sum(intensity)) %>%
    mutate(pass_filter = case_when(
        source == "D" & collection == "slaughter" & n >= 5000 ~ T,
        source == "D" & collection == "tube" & n >= 0 ~ T,
        source == "L" & collection == "slaughter" & n >= 0 ~ T,
        source == "W" & collection == "slaughter" & n >= 1000 ~ T,
        TRUE ~ FALSE
    ))
filter1 %>%
    ggplot(aes(n, sum / n, color = pass_filter, label = animal)) +
    facet_wrap(~ paste(source, collection), scales = "free") +
    geom_text(size = 2, hjust = 0, vjust = 0, color = "black") +
    geom_point(alpha = 0.3)

ggsave("results/filtered/plots/plot1.png")



# --- Perform filtering


proteome_intensities_filtered <- proteome_intensities %>%
    left_join(filter1, by = c("source", "animal", "collection")) %>%
    filter(pass_filter)


paste("before:", proteome_intensities %>% nrow())
paste("after:", proteome_intensities_filtered %>% nrow())

proteome_intensities_filtered %>%
    write_rds(snakemake@output["filtered"] %>% as.character())
