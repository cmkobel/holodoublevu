#!/usr/bin/env Rscript

# rm(list = ls())
library(tidyverse)
library(clusterProfiler)
library(patchwork) # Does clusterProfiler already require this?

source("workflow/scripts/utils.R")

# --- Inputs

metadata_file <- snakemake@input[["metadata"]] %>% as.character()
groups_file = snakemake@input[["groups"]] %>% as.character()
trait_modules_of_interest_file = snakemake@input[["trait_modules_of_interest"]] %>% as.character()
mm_ts_file = snakemake@input[["mm_ts"]] %>% as.character()

output_figures_file <- snakemake@output[["figures"]] %>% as.character()

# For debugging
if (interactive()) {
    metadata_file <- "resources/metadata_v1.6.tsv"
    groups_file <- "results/ig/both/imputed/groups.tsv"
    trait_modules_of_interest_file = "results/ig/both/wgcna/inspected/module_membership_trait_significance/trait_modules_of_interest.tsv"
    mm_ts_file = "results/ig/both/wgcna/inspected/module_membership_trait_significance/mm_ts_annotated_tax.rds"

    output_figures_file <- "results/ig/both/wgcna/inspected/contrasting_visualizations/plot.flag"

    figno_var <<- 1000
}

height <- 12
width <- 12

# pdf(generate_fig_name(output_rds_file, "dendro"), height = height, width = width)


metadata <- read_tsv(metadata_file)
groups <- read_tsv(groups_file)
trait_modules_of_interest = read_tsv(trait_modules_of_interest_file)
mm_ts = read_rds(mm_ts_file)


# --- Housekeeping



handful(mm_ts)

glimpse(mm_ts)


lapply(
    groups %>% filter(collection != "tube") %>% group_by(group_index) %>% group_split(), # i = (groups %>% group_by(group_index) %>% group_split())[[1]]
    function(i) {
        
        
        significant_species = trait_modules_of_interest %>%
            filter(group_index == i$group_index) %>%
            filter(trait == "vsplit") %>%
            filter(pvalue <= 0.05) %>%
            left_join(
                mm_ts %>%
                    filter(group_index == i$group_index),
                by = c("group_index", "trait", "module")
            ) %>%
            
            count(tax_genus) %>%
            mutate(n >= quantile(n, 1-0.0000005)) %>% # Hmmmm there is an error here??
            distinct(tax_genus)
        
        
        trait_module_order = trait_modules_of_interest %>% 
            filter(group_index == i$group_index) %>%
            filter(trait == "vsplit") %>%
            filter(pvalue <= 0.05) %>%
            arrange(coefficient) %>%
            distinct(module)
        
        if (nrow(trait_module_order) > 0) {
        
            plot_data = mm_ts %>%
                filter(group_index == i$group_index) %>%
                filter(trait == "vsplit") %>%
                count(tax_intermediate = tax_kingdom, tax_genus, module) %>%
                
                mutate(module = factor(module, levels = trait_module_order$module))
            
            plot_data %>%
                ggplot(aes(module, tax_genus, fill = n)) + 
                scale_fill_viridis_b(trans = "log") + # works, but has too many digits
                ggforce::facet_col(tax_intermediate ~ ., scales = "free_y", space = "free") +
                geom_tile() + 
                labs(
                    title = i$presentable,
                    subtitle = "Top biomass species",
                    caption = "Modules sorted by ascending trait (vsplit) correlation."
                )
            
            
        
            height_multiplier <- plot_data %>%
                distinct(tax_genus) %>%
                count(tax_genus) %>%
                pull(n) %>%
                sum()
            
            
            ggsave(generate_fig_name(output_figures_file, paste_("vsplit modules", i$presentable)), height = (height_multiplier / 7) + 2, width = 12)
            
        } else {
            plot_ = NULL
            #return(NULL)
        }
        
    }
)



