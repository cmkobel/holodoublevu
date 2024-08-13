#!/usr/bin/env Rscript

library(tidyverse)
library(missRanger)
source("workflow/scripts/utils.R")

# --- Inputs

threshold_present_in_number_of_samples_to_keep_in_imputation <- 4

metadata <- read_tsv(snakemake@input[["metadata"]] %>% as.character())
# metadata <- read_tsv("resources/metadata_v1.6.tsv")

proteome_intensities_raw <- read_rds(snakemake@input[["filtered"]] %>% as.character())
# proteome_intensities_raw <- read_rds("results/filtered/proteome_intensities.rds")


samples <- snakemake@params[["samples"]]
# samples <- c("D06T6S", "D08T6S", "D16T6S", "D25T6S", "D28T1S", "D28T2S", "D28T3S", "D28T4S", "D28T5S", "D28T6S", "D29T6S", "D35T6S", "D51T1S", "D51T2S", "D51T3S", "D51T4S", "D51T5S", "D51T6S", "D60T6S", "D61T6S", "D76T6S", "L06T6R", "L08T6R", "L16T6R", "L25T6R", "L28T6R", "L29T6R", "L35T6R", "L42T6R", "L51T6R", "L60T6R", "L61T6R", "L76T6R", "W06T6R", "W08T6R", "W16T6R", "W25T6R", "W29T6R", "W35T6R", "W42T6R", "W51T6R", "W60T6R", "W61T6R", "W76T6R") # luing example

paste("number of samples", length(samples))
paste(samples)

imputation_group <- tibble(sample = samples, included = T)
imputation_group %>% handful()
# --- Housekeeping

# Check that the wanted breeds are included.
paste("Breeds within selected samples:")
imputation_group %>%
    supacow_separate_sample() %>%
    left_join(metadata, by = "animal") %>%
    distinct(breed)

message("th")
# Define proteome intensities only with the samples wanted.
proteome_intensities <- proteome_intensities_raw %>%
    left_join(imputation_group, by = "sample") %>%
    filter(included) %>%
    select(-included) %>%
    supacow_separate_sample() %>%
    add_collection_column()

message("nt")
proteome_intensities %>%
    count(collection)



# --- Group by source and collection and remove samples that are not in at least 4

groups = proteome_intensities %>%
    group_by(source, collection) %>%
    summarize() %>%
    ungroup() %>%
    mutate(group_index = 1:n())

message("These are the groups and their indexes which we're going to use when imputing.")
groups

proteome_intensities_indexed <- proteome_intensities %>%
    left_join(groups, by = c("source", "collection")) %>%
    group_by(group_index) %>%
    identity()


# --- impute for each group


imputed <- lapply(
    group_split(proteome_intensities_indexed),
    function(i) {
        title_ <- i$group_index[[1]]
        message("imputing group index ", title_)


        i %>%
            count(protein) %>%
            ggplot(aes(n)) +
            geom_histogram() +
            labs(title = title_, subtitle = "samples per protein")

        worthy <- i %>%
            count(protein, name = "samples_per_protein") %>%
            filter(samples_per_protein >= threshold_present_in_number_of_samples_to_keep_in_imputation)


        i %>%
            supacow_paste_sample() %>%
            select(sample, protein, intensity) %>%
            left_join(worthy) %>%
            drop_na(samples_per_protein) %>%
            pivot_wider(id_cols = "protein", names_from = "sample", values_from = "intensity") %>%
            # head(5000) %>% # DEBUG!!!
            column_to_rownames("protein") %>%
            missRanger() %>%
            as_tibble(rownames = "protein") %>%
            mutate(group_index = title_) %>%
            pivot_longer(-c(protein, group_index), names_to = "sample", values_to = "intensity")

        # Filter based on the number of samples that the protein is in. I think the default that I used to use was 4.
        # Is there a neat way of counting that?
    }
) %>%
    bind_rows()


imputed %>%
    handful()

imputed %>%
    write_rds_and_tsv(snakemake@output[["imputed"]] %>% as.character())

imputation_group %>%
    write_tsv(snakemake@output[["imputation_group"]] %>% as.character())
