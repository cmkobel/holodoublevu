#!/usr/bin/env Rscript

# rm(list = ls())
library(tidyverse)
source("workflow/scripts/utils.R")

# --- Inputs

metadata_file <- snakemake@input[["metadata"]] %>% as.character()
groups_file = snakemake@input[["groups"]] %>% as.character()
net_results_file <- snakemake@input[["wgcna_modules"]] %>% as.character()
proteome2genome_file <- snakemake@input[["proteome2genome"]] %>% as.character()
# annotations_file = snakemake@input[["annotations"]] %>% as.character()


output_species_table_file <- snakemake@output[["species_table"]] %>% as.character()

# For debugging
if (F) {
    metadata_file <- "resources/metadata_v1.6.tsv"
    groups_file <- "results/ig/luing/imputed/groups.tsv"
    net_results_file <- "results/ig/luing/wgcna/modules.rds"
    proteome2genome_file <- "/glittertind/home/carl/PhD/26_proteomics_analysis/proteomic_integration/results/01_maps/01_map_proteome2samples_c6_v4.rds.gz"
    # annotations_file <- "results/annotation/annotation.emapper.annotations"

    output_species_table_file <- "results/ig/aberdeen/wgcna/inspected/species_table/species_table.rds"

    figno_var <<- 1000
}

height <- 12
width <- 12

# pdf(generate_fig_name(output_rds_file, "dendro"), height = height, width = width)


metadata <- read_tsv(metadata_file)
groups <- read_tsv(groups_file)
net_results <- read_rds(net_results_file)
# annotations_raw <- read_tsv(annotations_file, skip = 4, comment = "##", na = "-")
proteome2genome_raw <- read_rds(proteome2genome_file)



# --- Housekeeping


handful(metadata)
handful(groups)
handful(net_results[[1]]$net$colors)
handful(proteome2genome_raw)

proteome2genome_raw %>%
    head() %>%
    view()

proteome2genome = proteome2genome_raw %>%
    select(protein, sample, starts_with("tax_"))

# Does E look alright?
proteome2genome %>% filter(str_detect(protein, "dbE"))


proteome2genome %>%
    count(str_sub(protein, 1, 4)) %>%
    arrange(desc(n))

# Print one random from each.
proteome2genome %>%
    group_by(str_sub(protein, 1, 4)) %>%
    slice_sample(n = 1)


# Let's see if they map at all.



species_table <- lapply( # one group, e.g. "D, slaughter, 1"
    groups %>%
        rowwise() %>%
        group_split(), # Debug i = (groups %>% rowwise() %>% group_split())[[1]]
    function(i) {
        message("group ", paste(i, collapse = ", "))


        group_taxonomized <- net_results[[i$group_index]]$net$colors %>%
            as_tibble(rownames = "protein") %>%
            rename(module = value) %>%
            left_join(proteome2genome) %>%
            mutate(group_index = i$group_index[[1]]) %>%
            relocate(group_index)


        # Count NAs
        group_taxonomized %>%
            count(is.na(tax_kingdom)) # Wtf apparently none?? I'm coding better than I thought.

        group_taxonomized
    }
) %>%
    bind_rows() %>%
    # Clean up binomial name. (As it isn't otherwise uniform throughout the different databases.)
    mutate(
        tax_binomial = case_when(
            is.na(tax_species) ~ paste(tax_genus, "s??p."), # Unknown species
            tax_genus == str_sub(tax_species, 1, str_length(tax_genus)) ~ paste(tax_species), # Coded in genus only?
            TRUE ~ paste(tax_genus, tax_species) # Regular
        )
    ) %>%
    mutate(str_length(tax_genus))

species_table %>%
    select(tax_genus, tax_species, tax_binomial) %>%
    handful()


species_table %>%
    write_rds_and_tsv(output_species_table_file)


# --- Visualize

handful(species_table)


# 1 family
lapply(
    species_table %>%
        group_by(group_index) %>%
        group_split(), # i = (species_table %>% group_by(group_index) %>% group_split())[[1]]
    function(i) {
        j <- i %>%
            count(module, tax_kingdom, tax_family)

        j %>%
            ggplot(aes(module, tax_family, fill = n)) +
            scale_fill_viridis_b(begin = 0, end = .85) +
            ggforce::facet_col(tax_kingdom ~ ., scales = "free_y", space = "free") +
            geom_tile() +
            theme_bw() +
            labs(
                title = paste(
                    groups %>% filter(group_index == i$group_index[[1]]),
                    collapse = ", "
                ),
                subtitle = "Count of proteins per module, for taxonomical groups"
            )

        height_multiplier <- j %>%
            count(tax_family) %>%
            nrow()

        # ggsave(generate_fig_name(output_species_table_file, "tax_family"), height = (height_multiplier / 7) + 2, width = 10)
        ggsave(generate_fig_name(output_species_table_file, paste_("tax_family", filter(groups, group_index == i$group_index[[1]])$presentable)), height = (height_multiplier / 7) + 2, width = 10, limitsize = F)
    }
)


# 2 genus
lapply(
    species_table %>%
        group_by(group_index) %>%
        group_split(), # i = (species_table %>% group_by(group_index) %>% group_split())[[1]]
    function(i) {
        j <- i %>%
            count(module, tax_kingdom, tax_genus)

        j %>%
            ggplot(aes(module, tax_genus, fill = n)) +
            scale_fill_viridis_b(begin = 0, end = .85) +
            ggforce::facet_col(tax_kingdom ~ ., scales = "free_y", space = "free") +
            geom_tile() +
            theme_bw() +
            labs(
                title = paste(
                    groups %>% filter(group_index == i$group_index[[1]]),
                    collapse = ", "
                ),
                subtitle = "Count of proteins per module, for taxonomical groups"
            )

        height_multiplier <- j %>%
            count(tax_genus) %>%
            nrow()


        # ggsave(generate_fig_name(output_species_table_file, "tax_genus"), height = (height_multiplier / 7) + 2, width = 10)
        ggsave(generate_fig_name(output_species_table_file, paste_("tax_genus", filter(groups, group_index == i$group_index[[1]])$presentable)), height = (height_multiplier / 7) + 2, width = 10, limitsize = F)
    }
)


# 3 binomial
lapply(
    species_table %>%
        group_by(group_index) %>%
        group_split(), # i = (species_table %>% group_by(group_index) %>% group_split())[[1]]
    function(i) {
        j <- i %>%
            count(module, tax_intermediate = paste(tax_kingdom, tax_phylum, sep = ", "), tax_binomial)

        j %>%
            ggplot(aes(module, tax_binomial, fill = n)) +
            scale_fill_viridis_b(begin = 0, end = .85) +
            ggforce::facet_col(tax_intermediate ~ ., scales = "free_y", space = "free") +
            geom_tile() +
            theme_bw() +
            labs(
                title = paste(
                    groups %>% filter(group_index == i$group_index[[1]]),
                    collapse = ", "
                ),
                subtitle = "Count of proteins per module, for taxonomical groups"
            )

        height_multiplier <- j %>%
            count(tax_binomial) %>%
            nrow()


        ggsave(generate_fig_name(output_species_table_file, paste_("tax_binomial", filter(groups, group_index == i$group_index[[1]])$presentable)), height = (height_multiplier / 7) + 2, width = 10, limitsize = F)
    }
)
