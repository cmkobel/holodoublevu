#!/usr/bin/env Rscript

library(tidyverse)
library(WGCNA)
source("workflow/scripts/utils.R")

# --- Inputs

metadata <- read_tsv(snakemake@input[["metadata"]] %>% as.character())
# metadata <- read_tsv("resources/metadata_v1.6.tsv")

proteome_intensities_raw <- read_rds(snakemake@input[["imputed"]] %>% as.character())
# proteome_intensities_raw <- read_rds("results/luing/imputed/proteome_intensities.rds")

imputation_group = read_tsv(snakemake@input[["imputation_group"]] %>% as.character())
# imputation_group = "read_rds("results/luing/imputed/imputation_group.tsv")"

samples <- snakemake@params[["samples"]]
# samples <- c("D06T6S", "D08T6S", "D16T6S", "D25T6S", "D28T1S", "D28T2S", "D28T3S", "D28T4S", "D28T5S", "D28T6S", "D29T6S", "D35T6S", "D51T1S", "D51T2S", "D51T3S", "D51T4S", "D51T5S", "D51T6S", "D60T6S", "D61T6S", "D76T6S", "L06T6R", "L08T6R", "L16T6R", "L25T6R", "L28T6R", "L29T6R", "L35T6R", "L42T6R", "L51T6R", "L60T6R", "L61T6R", "L76T6R", "W06T6R", "W08T6R", "W16T6R", "W25T6R", "W29T6R", "W35T6R", "W42T6R", "W51T6R", "W60T6R", "W61T6R", "W76T6R") # luing example
