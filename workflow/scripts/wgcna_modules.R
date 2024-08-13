#!/usr/bin/env Rscript

library(tidyverse)
library(WGCNA)
source("workflow/scripts/utils.R")

# --- Inputs

metadata <- read_tsv(snakemake@input[["metadata"]] %>% as.character())
proteome_intensities <- read_rds(snakemake@input[["imputed"]] %>% as.character())
groups <- read_tsv(snakemake@input[["groups"]] %>% as.character())
samples <- snakemake@params[["samples"]]

# For debugging
if (F) {
    metadata <- read_tsv("resources/metadata_v1.6.tsv")
    proteome_intensities <- read_rds("results/ig/luing/imputed/proteome_intensities.rds")
    groups <- read_tsv("results/ig/luing/imputed/groups.tsv")
    samples <- c("D06T6S", "D08T6S", "D16T6S", "D25T6S", "D28T1S", "D28T2S", "D28T3S", "D28T4S", "D28T5S", "D28T6S", "D29T6S", "D35T6S", "D51T1S", "D51T2S", "D51T3S", "D51T4S", "D51T5S", "D51T6S", "D60T6S", "D61T6S", "D76T6S", "L06T6R", "L08T6R", "L16T6R", "L25T6R", "L28T6R", "L29T6R", "L35T6R", "L42T6R", "L51T6R", "L60T6R", "L61T6R", "L76T6R", "W06T6R", "W08T6R", "W16T6R", "W25T6R", "W29T6R", "W35T6R", "W42T6R", "W51T6R", "W60T6R", "W61T6R", "W76T6R") # luing example
}

# --- Housekeeping


metadata %>% handful()
proteome_intensities %>% handful()
groups %>% handful()
samples %>% handful()



# --- Dendrogram


i <- (proteome_intensities %>% # This is a debug setup.
    group_by(group_index) %>%
    group_split())[[1]]

hclust_naive <- lapply(
    (proteome_intensities %>%
        group_by(group_index) %>%
        group_split()),
    function(i) {
        i %>% handful()

        current_group_index <- i$group_index[[1]]
        message("group index ", current_group_index)

        groups %>% filter(group_index == current_group_index)

        wide <- i %>%
            pivot_wider(id_cols = c(sample, group_index), names_from = "protein", values_from = "intensity") %>%
            select(-group_index)

        sample_tree <- wide %>%
            column_to_rownames("sample") %>%
            dist() %>%
            hclust()

        trait_colors_data <- wide %>%
            select(sample) %>%
            supacow_separate_sample(keep_sample = T) %>%
            select(animal, sample) %>%
            left_join(metadata, by = "animal") %>%
            column_to_rownames("sample") %>%
            # select(where(is.numeric))
            mutate(across(everything(), ~ .x %>%
                as.factor() %>%
                as.numeric())) %>%
            identity()

        list(
            sample_tree = sample_tree,
            trait_colors_data = trait_colors_data,
            group_index = current_group_index
        )
    }
)

lapply(
    hclust_naive,
    function(i) {
        plotDendroAndColors(
            i$sample_tree,
            i$trait_colors_data %>% numbers2colors(),
            groupLabels = names(i$trait_colors_data),
            main = paste(
                groups %>% filter(group_index == current_group_index),
                collapse = ", "
            )
        )
    }
)
