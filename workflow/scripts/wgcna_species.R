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
#Hvis jeg læser trait_modules_of_interest.tsv ind i pathway or species kan jeg highlighte i plotsene
trait_modules_of_interest_file = snakemake@input[["trait_modules_of_interest"]] %>% as.character()

output_species_table_file <- snakemake@output[["species_table"]] %>% as.character()



# For debugging
if (interactive()) {
    metadata_file <- "resources/metadata_v1.6.tsv"
    groups_file <- "results/ig/both/imputed/groups.tsv"
    net_results_file <- "results/ig/both/wgcna/modules.rds"
    proteome2genome_file <- "~/PhD/26_proteomics_analysis/proteomic_integration/results/01_maps/01_map_proteome2samples_c6_v4.rds.gz"
    # annotations_file <- "results/annotation/annotation.emapper.annotations"
    trait_modules_of_interest_file = "results/ig/both/wgcna/inspected/module_membership_trait_significance/trait_modules_of_interest.tsv"

    output_species_table_file <- "results/ig/aberdeen/wgcna/species_table/species_table.rds"

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
trait_modules_of_interest = read_tsv(trait_modules_of_interest_file)



# --- Housekeeping


handful(metadata)
handful(groups)
handful(net_results[[1]]$net$colors)
handful(proteome2genome_raw)

proteome2genome_raw %>%
    head()

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


trait_modules_of_interest %>% 
    filter(trait == "vsplit") %>%
    handful()


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
    
    
    # This is to fix up some incoherencies that should rather be fixed in the prot2genome mapping.
    # If family is missing, take order, etc.
    mutate(
        tax_family = coalesce(tax_family, tax_order),
        tax_genus = coalesce(tax_genus, tax_family)
        ) %>%
    
    # Another incoherence is that some dbB's are missing the bovine tax id. This should also be fixed outside in the original source file.
    # Eukaryota Animalia Chordata Mammalia Artiodactyla Bovidae Bos taurus Bos taurus
    mutate(
        tax_domain  = case_when(is.na(sample) & str_detect(protein, "^dbB\\|") ~ "Eukaryota", TRUE ~ tax_domain ),
        tax_kingdom = case_when(is.na(sample) & str_detect(protein, "^dbB\\|") ~ "Animalia", TRUE ~ tax_kingdom),
        tax_phylum  = case_when(is.na(sample) & str_detect(protein, "^dbB\\|") ~ "Chordata", TRUE ~ tax_phylum ),
        tax_class   = case_when(is.na(sample) & str_detect(protein, "^dbB\\|") ~ "Mammalia", TRUE ~ tax_class  ),
        tax_order   = case_when(is.na(sample) & str_detect(protein, "^dbB\\|") ~ "Artiodactyla", TRUE ~ tax_order  ),
        tax_family  = case_when(is.na(sample) & str_detect(protein, "^dbB\\|") ~ "Bovidae", TRUE ~ tax_family ),
        tax_genus   = case_when(is.na(sample) & str_detect(protein, "^dbB\\|") ~ "Bos", TRUE ~ tax_genus  ),
        tax_species = case_when(is.na(sample) & str_detect(protein, "^dbB\\|") ~ "taurus", TRUE ~ tax_species)
    ) %>%
    
    # Clean up binomial name. (As it isn't otherwise uniform throughout the different databases.)
    mutate(
        tax_binomial = case_when(
            is.na(tax_species) ~ paste(tax_genus, "sp."), # Unknown species
            tax_genus == str_sub(tax_species, 1, str_length(tax_genus)) ~ paste(tax_species), # Coded in genus only?
            TRUE ~ paste(tax_genus, tax_species) # Regular
        )
    )




species_table %>%
    select(tax_genus, tax_species, tax_binomial) %>%
    handful()


# save this file because now it is clean
if (!interactive()) {
    species_table %>%
        write_rds_and_tsv(output_species_table_file)
}


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
        
        dist_ = j %>%
            pivot_wider(id_cols = module, names_from = tax_family, values_from = n, values_fill = 0) %>%
            drop_na(module) %>%
            column_to_rownames(var = "module") %>%
            dist(method = "binary") %>%
            hclust()

        j %>%
            mutate(module = factor(module, levels = dist_$labels[dist_$order])) %>%
            ggplot(aes(module, tax_family, fill = n)) +
            scale_fill_viridis_b(begin = 0, end = .85, trans = "log") +
            ggforce::facet_col(tax_kingdom ~ ., scales = "free_y", space = "free") +
            geom_tile() +
            theme_bw() +
            labs(
                title = filter(groups, group_index == i$group_index[[1]])$presentable,
                subtitle = "Count of proteins per module, for taxonomical groups",
                caption = "Modules sorted by binary distance."
            ) + 
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
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
        
        dist_ = j %>%
            pivot_wider(id_cols = module, names_from = tax_genus, values_from = n, values_fill = 0) %>%
            drop_na(module) %>%
            column_to_rownames(var = "module") %>%
            dist(method = "binary") %>%
            hclust()
    
        # modules_to_highlight = trait_modules_of_interest %>%
        #     filter(trait == "vsplit" & group_index == i$group_index[[1]])

        j %>%
            mutate(module = factor(module, levels = dist_$labels[dist_$order])) %>%
            ggplot(aes(module, tax_genus, fill = n)) +
            scale_fill_viridis_b(begin = 0, end = .85, trans = "log") +
            ggforce::facet_col(tax_kingdom ~ ., scales = "free_y", space = "free") +
            # geom_tile(
            #     mapping = aes(str_extract(module, "\\d+"), "trait_vsplit", fill = abs(coefficient)),
            #     data = modules_to_highlight
            #     ) +
            geom_tile() + 
            theme_bw() +
            labs(
                title = filter(groups, group_index == i$group_index[[1]])$presentable,
                subtitle = "Count of proteins per module, for taxonomical groups",
                caption = "Modules sorted by binary distance."
            ) +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
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
        
        dist_ = j %>%
            pivot_wider(id_cols = module, names_from = tax_binomial, values_from = n, values_fill = 0) %>%
            #pivot_wider(id_cols = module, names_from = tax_binomial, values_from = n, values_fn = list) %>%
            drop_na(module) %>% # 
            column_to_rownames(var = "module") %>%
            dist(method = "binary") %>% # "euclidean"27 (def), "maximum"28, "manhattan", "canberra"30, "binary"25 or "minkowski"26.
            hclust()
        
        j %>%
            mutate(module = factor(module, levels = dist_$labels[dist_$order])) %>%
            ggplot(aes(module, tax_binomial, fill = n)) +
            scale_fill_viridis_b(begin = 0, end = .85, trans = "log") +
            ggforce::facet_col(tax_intermediate ~ ., scales = "free_y", space = "free") +
            geom_tile() +
            theme_bw() +
            labs(
                title = filter(groups, group_index == i$group_index[[1]])$presentable,
                subtitle = "Count of proteins per module, for taxonomical groups",
                caption = "Modules sorted by binary distance."
            ) +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            )

        height_multiplier <- j %>%
            count(tax_binomial) %>%
            nrow()


        ggsave(generate_fig_name(output_species_table_file, paste_("tax_binomial", filter(groups, group_index == i$group_index[[1]])$presentable)), height = (height_multiplier / 7) + 2, width = 10, limitsize = F)
        
    }
)
