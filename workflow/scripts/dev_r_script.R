#!/usr/bin/env Rscript

library(tidyverse)
source("workflow/scripts/utils.R")

# --- Inputs

proteome_intensities_raw <- read_rds(snakemake@input["proteome_intensities"] %>% as.character())
# proteome_intensities_raw <- read_rds("/glittertind/home/carl/PhD/26_proteomics_analysis/proteomic_integration/results/00_uniform_data/proteomics_long_v2.rds")
