#!/usr/bin/env Rscript

library(tidyverse)
library(missRanger)
source("workflow/scripts/utils.R")

# --- Inputs
metadata <- read_tsv(snakemake@input[["metadata"]] %>% as.character())
# metadata <- read_tsv("resources/metadata_v1.6.tsv")

proteome_intensities_raw <- read_rds(snakemake@input[["filtered"]] %>% as.character())
# proteome_intensities_raw <- read_rds("results/filtered/proteome_intensities.rds")

samples <- snakemake@params[["samples"]]
# samples <- c("D06T6S", "W06T6R", "L06T6R", "D08T6S", "W08T6R", "L08T6R", "D06T6S", "W06T6R", "L06T6R", "D08T6S", "W08T6R", "L08T6R", "D16T6S", "W16T6R", "L16T6R", "D25T6S", "W25T6R", "L25T6R", "D28T6S", "L28T6R", "D28T1S", "D28T3S", "D28T4S", "D28T5S", "D28T2S", "D29T6S", "W29T6R", "L29T6R", "D35T6S", "W35T6R", "L35T6R", "W42T6R", "L42T6R", "D51T6S", "W51T6R", "L51T6R", "D51T1S", "D51T2S", "D51T3S", "D51T4S", "D51T5S", "D60T6S", "W60T6R", "L60T6R", "D61T6S", "W61T6R", "L61T6R", "D76T6S", "W76T6R", "L76T6R")

tibble(samples) %>% show_some()

# --- Housekeeping
