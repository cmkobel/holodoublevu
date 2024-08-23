
#!/usr/bin/env Rscript
# rm(list = ls())
library(tidyverse)
library(clusterProfiler)
library(patchwork) # Does clusterProfiler already require this?

source("workflow/scripts/utils.R")

# --- Inputs

metadata_file <- snakemake@input[["metadata"]] %>% as.character()
groups_file = snakemake@input[["groups"]] %>% as.character()
net_results_file <- snakemake@input[["wgcna_modules"]] %>% as.character()
annotations_file = snakemake@input[["annotations"]] %>% as.character()
kegg_data_file = snakemake@input[["kegg_data"]] %>% as.character()
trait_modules_of_interest_file = snakemake@input[["trait_modules_of_interest"]] %>% as.character()


output_pathway_enrichment_file <- snakemake@output[["pathway_enrichment"]] %>% as.character()

# For debugging
if (interactive()) {
    metadata_file <- "resources/metadata_v1.6.tsv"
    groups_file <- "results/ig/both/imputed/groups.tsv"
    net_results_file <- "results/ig/both/wgcna/modules.rds"
    annotations_file <- "results/annotation/annotation.emapper.annotations"
    kegg_data_file <- "resources/kegg_data.tsv"
    trait_modules_of_interest_file = "results/ig/both/wgcna/inspected/module_membership_trait_significance/trait_modules_of_interest.tsv"
    

    output_pathway_enrichment_file <- "results/ig/aberdeen/wgcna/pathway_enrichment/pathway_enrichment.rds"

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
trait_modules_of_interest = read_tsv(trait_modules_of_interest_file)


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


if (!interactive()) {
    annotations %>%
        write_rds_and_tsv(snakemake@output[["annotations_clean"]][[1]])
}


if (annotations$protein %>% length() != annotations$protein %>%
    unique() %>%
    length()) {
    stop("Why are the protein names not unique? Shouldn't they be?")
}


# --- module pathways
# Count the number of proteins from each pathway in each module.
# Jeg fik lige en interessant ide. Hvorfor er det egentlig at jeg har fokuseret så meget på at kigge på enrichment? Hvad med at kigge på antal ortologer for en pathway istedet? Er det ikke lige så godt som at køre en masse hypergeometriske test, bare at se hvad der er til stede? det er i hvert fald hurtigere. Og antallet måske også lettere at fortolke end en p-værdi?

# Hierarchy used for sorting
pathway_hierarchy <- kegg_data %>%
    distinct(pathway_class = class, pathway_group = group, pathway) %>%
    arrange(pathway_class, pathway_group, pathway)
handful(pathway_hierarchy)


# Geneset
term2gene <- kegg_data %>%
    select(term = pathway, gene = ko_id)
term2gene %>% handful()


# Add the annotation to the proteins of a single layer.
# This makes sense, since one protein is only present in one module.

# This is the most fucked up un-debuggable nested lapply on group_split()-ted tables I've ever done.


pe_analyses <- lapply( # one group, e.g. "D, slaughter, 1"
    groups %>%
        rowwise() %>%
        group_split(), # i = (groups %>% rowwise() %>% group_split())[[1]]
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

if (!interactive()) {
    pe_analyses %>% write_rds_and_tsv(output_pathway_enrichment_file)
}


# --- plots

# Simple tile showing pathways and modules


lapply(
    pe_analyses %>%
        group_by(group_index) %>%
        group_split(), # i = (pe_analyses %>% group_by(group_index) %>% group_split())[[1]]
    function(i) {
        j <- i %>%
            left_join(pathway_hierarchy, by = "pathway") %>%
            # slice_sample(n = 100) %>%
            filter(!str_detect(pathway_class, "^09160|09190")) # Remove boring pathways (human diseases and whatnot)
        
        dist_ = j %>%
            pivot_wider(id_cols = module, names_from = pathway, values_from = p.adjust, values_fill = 0) %>%
            drop_na(module) %>%
            column_to_rownames(var = "module") %>%
            dist(method = "binary") %>%
            hclust()

        
        plot_top = j %>%
            mutate(module = factor(module, levels = dist_$labels[dist_$order])) %>% # sort modules on distance
            mutate(pathway = factor(pathway, levels = pathway_hierarchy$pathway)) %>% # sort pathways hierarchically
            ggplot(aes(module, pathway, fill = p.adjust)) +
            # facet_wrap(~pathway_class, space = "free") +
            scale_fill_viridis_b(begin = 0, end = .85) +
            scale_x_discrete(drop = FALSE) +
            ggforce::facet_col(pathway_class ~ ., scales = "free_y", space = "free") +
            geom_tile() +
            theme_bw() +
            labs(
                title = filter(groups, group_index == i$group_index[[1]])$presentable,
                subtitle = "Pathway enrichment analysis on WGCNA modules. Only p.adjust significant values reported."
            ) + 
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            )
        
        
        j
        
        # I wish to add a plot underneath that shows a selected trait. Because that would make it easy to look for patterns that could explain stuff.
        plot_bottom = trait_modules_of_interest %>% 
            #filter(trait == "vsplit") %>% 
            mutate(module = factor(module, levels = paste0("ME", dist_$labels[dist_$order]))) %>% # sort modules on distance
            ggplot(aes(module,  trait, fill = coefficient)) +
            scale_fill_viridis_b(begin = 0, end = .85, option = "H") +
            scale_x_discrete(drop = FALSE) +
            geom_tile() +
            theme(
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
            ) + 
            labs(
                caption = "Modules clustered by binary distance of pathways."
            )
        
        
        plot_top / plot_bottom + 
            patchwork::plot_layout(heights = c(10, 1))

        height_multiplier <- j %>%
            count(pathway) %>%
            nrow()

        ggsave(generate_fig_name(output_pathway_enrichment_file, paste_("pathway", filter(groups, group_index == i$group_index[[1]])$presentable)), height = (height_multiplier / 7) + 2, width = 12)
    }
)
