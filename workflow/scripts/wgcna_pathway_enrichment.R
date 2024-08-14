#!/usr/bin/env Rscript
# rm(list = ls())
library(tidyverse)
library(clusterProfiler)
source("workflow/scripts/utils.R")

# --- Inputs

metadata_file <- snakemake@input[["metadata"]] %>% as.character()
groups_file = snakemake@input[["groups"]] %>% as.character()
net_results_file <- snakemake@input[["wgcna_modules"]] %>% as.character()
annotations_file = snakemake@input[["annotations"]] %>% as.character()
kegg_data_file = snakemake@input[["kegg_data"]] %>% as.character()

output_pathway_enrichment_file <- snakemake@output[["pathway_enrichment"]] %>% as.character()

# For debugging
if (F) {
    metadata_file <- "resources/metadata_v1.6.tsv"
    groups_file <- "results/ig/luing/imputed/groups.tsv"
    net_results_file <- "results/ig/luing/wgcna/modules.rds"
    annotations_file <- "results/annotation/annotation.emapper.annotations"
    kegg_data_file <- "resources/kegg_data.tsv"
    figno_var <<- 1000
}

height <- 12
width <- 12

# pdf(generate_fig_name(output_rds_file, "dendro"), height = height, width = width)


metadata <- read_tsv(metadata_file)
groups <- read_tsv(groups_file)
net_results <- read_rds(net_results_file)
annotations_raw <- read_tsv(annotations_file, skip = 4, comment = "##", na = "-")
kegg_data_raw <- read_tsv(kegg_data_file)


# --- Housekeeping

kegg_data <- kegg_data_raw %>%
    rename(ko_id = ortholog_ko)

annotations_raw %>% handful()
annotations_all <- annotations_raw %>%
    rename(query = `#query`) %>%
    # select(`#query`, starts_with("KEGG_"))
    # head(100) %>%
    separate_wider_regex(
        query,
        c(
            database_short = "^db[A-Z]",
            "\\d*\\|", # Crucial to not be greedy, but rather accept only a digit.
            protein_short = "[^ ]+", # Anything but a space.
            ".*$"
        ),
        cols_remove = F
    ) %>%
    mutate(protein = paste0(database_short, "|", protein_short)) %>%
    relocate(protein)

annotations_all %>% handful()

annotations <- annotations_all %>%
    select(-database_short, -protein_short, -query)



if (annotations$protein %>% length() != annotations$protein %>%
    unique() %>%
    length()) {
    stop("Why are the protein names not unique? Shouldn't they be?")
}


# --- module pathways
# Count the number of proteins from each pathway in each module.
# Jeg fik lige en interessant ide. Hvorfor er det egentlig at jeg har fokuseret så meget på at kigge på enrichment? Hvad med at kigge på antal ortologer for en pathway istedet? Er det ikke lige så godt som at køre en masse hypergeometriske test, bare at se hvad der er til stede? det er i hvert fald hurtigere. Og antallet måske også lettere at fortolke end en p-værdi?

# Geneset
term2gene <- kegg_data %>%
    select(term = pathway, gene = ko_id)
term2gene %>% handful()


# Add the annotation to the proteins of a single layer.
# This makes sense, since one protein is only present in one module.



pe_analyses <- lapply( # one group, e.g. "D, slaughter, 1"
    groups %>%
        rowwise() %>%
        group_split(), # Debug i = (groups %>% rowwise() %>% group_split())[[1]]
    function(i) {
        message("group ", paste(i, collapse = ", "))


        group_annotated <- net_results[[i$group_index]]$net$colors %>%
            as_tibble(rownames = "protein") %>%
            rename(module = value) %>%
            left_join(annotations) %>%
            group_by(module)


        lapply( # One module e.g. "0"
            group_annotated %>% group_split(), # Debug: # j = (group_annotated %>% group_split())[[1]]
            function(j) {
                message("group ", paste(i, collapse = ", "), ", module ", j$module[[1]])
                layer_ready <- j %>%
                    select(protein, module, KEGG_ko) %>%
                    drop_na(KEGG_ko) %>%
                    mutate(ko_extraction = str_split(KEGG_ko, ",")) %>%
                    unnest(ko_extraction) %>%
                    mutate(ko = str_extract(ko_extraction, "ko:(K\\d+)", group = 1)) %>%
                    select(-KEGG_ko, -ko_extraction)


                # Ready to perform enrichment analysis.

                gene <- layer_ready %>%
                    drop_na(ko) %>%
                    pull(ko)

                clusterProfiler::enricher(
                    gene,
                    TERM2GENE = term2gene
                ) %>%
                    as_tibble() %>%
                    select(-Description) %>% # Redundant with ID
                    rename(pathway = ID) %>%
                    mutate(
                        group_index = i$group_index,
                        module = layer_ready$module[[1]]
                    ) %>%
                    relocate(group_index, module)
            }
        ) %>% bind_rows()
    }
) %>% bind_rows()

pe_analyses %>% write_rds_and_tsv(output_pathway_enrichment_file)
