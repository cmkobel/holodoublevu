#!/usr/bin/env Rscript

# rm(list = ls())
library(tidyverse)
library(WGCNA)
library(patchwork)
source("workflow/scripts/utils.R")

# --- Inputs

metadata_file <- snakemake@input[["metadata"]] %>% as.character()
net_results_file <- snakemake@input[["wgcna_modules"]] %>% as.character()
groups_file = snakemake@input[["groups"]] %>% as.character()
output_rds_file <- snakemake@output[["mod2mod"]] %>% as.character()
output_module_membership_trait_significance_file = snakemake@output[["module_membership_trait_significance"]] %>% as.character()


# For debugging
if (interactive()) {
    metadata_file <- "resources/metadata_v1.6.tsv"
    net_results_file <- "results/ig/both/wgcna/modules.rds"
    groups_file <- "results/ig/both/imputed/groups.tsv"
    output_rds_file <- "results/ig/both/wgcna/inspected/module2module.rds"
    output_module_membership_trait_significance_file = "results/ig/both/wgcna/inspected/module_membership_trait_significance/module_membership_trait_significance.rds"
    figno_var <<- 1000
}

height = 12
width <- 12

metadata <- read_tsv(metadata_file)
net_results <- read_rds(net_results_file)
groups <- read_tsv(groups_file)


# --- Housekeeping


metadata %>% handful()
net_results[[1]] %>% str()
net_results[[1]]$net$colors %>% table()

net_results[[1]]$net$colors %>%
    enframe(name = "protein_name_long", value = "module") %>%
    group_by(module) %>%
    count()

net_results[[1]]$net$MEs %>%
    as_tibble(rownames = "sample")

# --- Proteins per module


lapply(
    net_results,
    function(i) { # i = net_results[[1]]
        
        
        colors <- i$net$colors %>%
            as_tibble(rownames = "protein") %>%
            rename(module = value)
            
        layer_stats = colors %>%
            count(module)
        
        

        hist_ <- layer_stats %>%
            ggplot(aes(n)) +
            geom_histogram() +
            labs(
                title = paste("Proteins per module:", filter(groups, group_index == i$group_index)$presentable), 
                subtitle = "Histogram", 
                x = "n [proteins]", 
                y = "count [modules]"
            )
                
        

        column_ <- layer_stats %>%
            ggplot(aes(module, n)) +
            geom_col() +
            labs(
                subtitle = "Column plot", x = "module name", y = "n [proteins]",
                caption = paste("Number of proteins:", nrow(colors), "\nCount of modules:", colors$module %>% unique() %>% length())
            )
        
                

        (hist_ / column_)

        ggsave(generate_fig_name(output_rds_file, paste_("count", filter(groups, group_index == i$group_index)$presentable)))
    }
)


# --- Basic dendro viz




lapply(
    net_results, # i <- net_results[[1]]
    function(i) {
        pdf(generate_fig_name(output_rds_file, paste_("dendro", filter(groups, group_index == i$group_index)$presentable)), height = height, width = width)

        plotDendroAndColors(
            i$net$dendrograms[[1]],
            tibble(index = i$net$blockGenes[[1]]) %>%
                left_join(
                    i$net$colors %>%
                        enframe() %>%
                        mutate(index = 1:n())
                ) %>%
                select(-index) %>%
                deframe(),
            main = paste(filter(groups, group_index == i$group_index)$presentable, "(block 1 only)"),
            # main = "Dendro (block 1 only)",
            dendroLabels = F
        )
        dev.off()

        # ggsave(generate_fig_name(output_rds_file, "dendro"))
    }
)




# --- pheatmaps



# --- ## Visualizing module eigengenes


lapply(
    net_results, # [1], # select one only with [n] where in is in 1:6.
    function(i) {
        pdf(generate_fig_name(output_rds_file, paste_("me", filter(groups, group_index == i$group_index)$presentable)), height = height, width = width)


        pheatmap::pheatmap(
            i$net$MEs,
            main = paste(
                groups %>% filter(group_index == i$group_index),
                collapse = ", "
            )
        )
        dev.off()
    }
)



# --- Cross correlation of modules.

# Pick the layers that represent:
# - digesta x wall
# - digesta x liver
# - liver x wall

layer_digesta <- filter(groups, source == "D" & collection == "slaughter") %>% pull(group_index)
layer_wall <- filter(groups, source == "W") %>% pull(group_index)
layer_liver <- filter(groups, source == "L") %>% pull(group_index)



wanted_keys <- list(
    # Aberdeen Angus X
    list(A = layer_digesta, B = layer_wall, plot_title = "digestaXwall"), # digesta x wall
    list(A = layer_digesta, B = layer_liver, plot_title = "digestaXliver"), # digesta x liver
    list(A = layer_liver, B = layer_wall, plot_title = "liverXwall") # liver x wall
)


bicor_module2module <- lapply(
    wanted_keys,
    function(i) { # i = list(A = 1, B = 7) # DEBUG
        message(i$plot_title)
        pdf(generate_fig_name(output_rds_file, paste_("axis", groups$imputation_group[[1]], i$plot_title)), height = height, width = width)



        return_value <- list()

        message("For: ", i$A, " x ", i$B)

        setHook("grid.newpage", function() grid::pushViewport(grid::viewport(x = 1, y = 1, width = 0.9, height = 0.9, name = "vp", just = c("right", "top"))), action = "prepend")


        # Due to possible outlier removal, we need to find the overlapping samples and make sure that they are on the correct rows.
        inner_animals <- inner_join(
            tibble(
                sample = rownames(net_results[[i$A]][["net"]][["MEs"]])
            ) %>%
                mutate(animal = str_sub(sample, 2, 3)),
            tibble(
                sample = rownames(net_results[[i$B]][["net"]][["MEs"]])
            ) %>%
                mutate(animal = str_sub(sample, 2, 3)),
            by = "animal",
            suffix = c("_A", "_B")
        )


        # This is based on the index positions of a list which is dangerous.
        # We're using the inner_animals table to get only the rows that we need for comparison.

        left_hand_side <- left_join(
            inner_animals %>%
                select(animal, sample = sample_A),
            net_results[[i$A]][["net"]][["MEs"]] %>%
                rownames_to_column("sample")
        ) %>%
            select(-animal) %>%
            column_to_rownames("sample")

        right_hand_side <- left_join(
            inner_animals %>%
                select(animal, sample = sample_B),
            net_results[[i$B]][["net"]][["MEs"]] %>%
                rownames_to_column("sample")
        ) %>%
            select(-animal) %>%
            column_to_rownames("sample")


        message("left dimension: ", nrow(left_hand_side), " ", ncol(left_hand_side))
        message("right dimension: ", nrow(right_hand_side), " ", ncol(right_hand_side))


        bicor_results <- bicorAndPvalue(
            left_hand_side,
            right_hand_side
        )


        presentable_A <- filter(groups, group_index == i$A)$presentable
        presentable_B <- filter(groups, group_index == i$B)$presentable

        pheatmap::pheatmap(
            bicor_results$bicor,
            display_numbers = bicor_results$p %>% gtools::stars.pval(),
            # main = paste0(
            #     paste(keys[c(i$A, i$B),] %>% pull(breed) %>% unique(), collapse = " / "), # breed
            #
            #     "\n bicor(",
            #     keys[i$B,]$source,
            #     ", ",
            #     keys[i$A,]$source,
            #     ")"
            # )
            main = paste0(
                "bicor(", i$A, ", ", i$B, ")"
            )
        )

        message("4")
        setHook("grid.newpage", NULL, "replace")
        grid::grid.text(presentable_B, y = -0.07, gp = grid::gpar(fontsize = 12))
        grid::grid.text(presentable_A, x = -0.07, rot = 90, gp = grid::gpar(fontsize = 12))


        dev.off()
        list(
            "A" = i$A,
            "B" = i$B,
            "presentable_A" = presentable_A,
            "presentable_B" = presentable_B,
            "bicor" = bicor_results
        )
    }
)


if (!interactive()) {
    bicor_module2module %>%
        write_rds(output_rds_file)
}


# We can extract significant module pairs from the bicor_module2module (axis) table holo module pairs or axis module pairs? axis couple? holo module couple?

    
axis_couples = base::lapply(
    bicor_module2module,
    function(x) { # x = bicor_module2module[1]
    
        x$bicor$p %>% # p, the Student p-values corresponding to the calculated correlations
            as_tibble(rownames = "module_A") %>% 
            
            pivot_longer(-module_A, names_to = "module_B", values_to = "pvalue") %>%  
            
            filter(pvalue < 0.05) %>% 
            # This is stupid, but it works. "Tidy data", right? And then I can't fuck it up by having the wrong "wanted_keys" table laying around in a corner.
            mutate(
                module_A_index = x$A,
                module_B_index = x$B,
                presentable_A = x$bicor$presentable_A,
                presentable_B = x$bicor$presentable_B
            ) %>% 
            
            # Add the bicor values. This was the simplest way I could think of. 
            left_join(
                x$bicor$bicor %>% # bicor, the calculated correlations
                    as_tibble(rownames = "module_A") %>% 
                    pivot_longer(-module_A, names_to = "module_B", values_to = "bicor"),
                by = c("module_A", "module_B")
            ) %>% 
            identity()
    }
)  %>% 
bind_rows() %>% 

# Since the names (i.e. "ME3", "ME4") are overlapping between sources, we must include the group index in the module name.
mutate(
    module_A_complete = paste(module_A_index, module_A, sep = "|"),
    module_B_complete = paste(module_B_index, module_B, sep = "|"),
    
) 



# --- Correlating to phenotypes




# By converting this before the for loop of each layer, I know that the factors will be comparable between layers. Makes visual interpretation easier.
metadata_numeric <- metadata %>%
    mutate(across(
        !is.numeric & !matches("animal"),
        ~ .x %>%
            as.factor() %>%
            as.numeric()
    ))

stopifnot(all(metadata$animal == metadata_numeric$animal))

bind_cols(
    metadata,
    metadata_numeric %>%
        rename_at(vars(-animal), function(x) {paste0(x, "_numeric")}) %>%
        select(-animal)
) %>%
    pivot_longer(-animal, values_transform = as.character) %>%
    mutate("_numeric" = str_detect(name, "_numeric$")) %>%
    mutate(name = str_remove(name, "_numeric$")) %>% 
    pivot_wider(id_cols = c(animal, name), names_from = `_numeric`, values_from = value) %>%
    select(-animal) %>%
    filter(`TRUE` != `FALSE`) %>% # The ones that were already numeric
    distinct() %>%
    arrange(name) %>% 
    write_tsv(paste0(dirname(output_rds_file), "/trait_key.tsv"))




# You will get problems with duplicate keys if you run this with the tube samples.

pheno = lapply(
    net_results[filter(groups, collection != "tube")$group_index], # x = net_results[filter(groups, collection != "tube")$group_index][[1]]
    function(x) {
        pdf(generate_fig_name(output_rds_file, paste_("pheno", filter(groups, group_index == x$group_index)$presentable)), height = height, width = width)

        bicor_result <- bicorAndPvalue(
            x = x$net$MEs,
            y = x$net$MEs %>%
                rownames_to_column("sample") %>%
                select(sample) %>%
                transmute(animal = str_sub(sample, 2, 3)) %>%
                left_join(
                    metadata_numeric,
                    by = "animal"
                ) %>%
                select(where(~ var(.x, na.rm = T) > 0)) %>%
                column_to_rownames("animal"),
            use = "pairwise.complete.obs"
        )
        
        pheatmap::pheatmap(
            mat = bicor_result$bicor,
            display_numbers = bicor_result$p %>%
                gtools::stars.pval(),
            main = filter(groups, group_index == x$group_index)$presentable
        ) %>%
            print()

        dev.off()
        
        list(
            bicor = bicor_result$bicor,
            p = bicor_result$p,
            group_index = x$group_index
        )
    }
)

# I want to extract the modules that are correlated with a trait


trait_modules_of_interest = lapply(
    pheno,
    function(i) { # i = pheno[[1]]
    a = i$bicor %>% as_tibble(rownames = "module") %>%
        pivot_longer(-module, names_to = "trait", values_to = "coefficient")
    b = i$p %>% as_tibble(rownames = "module") %>%
        pivot_longer(-module, names_to = "trait", values_to = "pvalue")
    
    stopifnot(all(a$module == b$module))
    stopifnot(all(a$trait == b$trait))
    
    bind_cols(a, b %>% select(pvalue)) %>%
        filter(pvalue <= 0.05) %>%
        mutate(group_index = i$group_index) %>%
        relocate(group_index, trait) %>%
        arrange(group_index, trait)
    
    }
) %>% bind_rows()

if (!interactive()) {
    trait_modules_of_interest %>%
        write_tsv(paste0(dirname(output_module_membership_trait_significance_file), "/trait_modules_of_interest.tsv"))
}


# Module membership and gene significance.
# Formerly traits_modules
mm_gs = lapply(
    net_results, # x = net_results[[1]]
    function(x) {
        
        # Do not process the tube samples.
        if (filter(groups, group_index == x$group_index)$collection == "tube") {
            return(NULL)
        }
        
        # Correlate proteins and modules (kME) including pvalues. I really need a way to be able to pick out which exact proteins are of significance. Otherwise there is just too much data.
        protein_module_raw = corAndPvalue(
            x$datExpr,
            x$net$MEs
        )
        
        a = protein_module_raw$cor %>%
            as_tibble(rownames = "protein") %>%
            pivot_longer(-protein, names_to = "module", values_to = "coefficient")
        
        b = protein_module_raw$p %>%
            as_tibble(rownames = "protein") %>%
            pivot_longer(-protein, names_to = "module", values_to = "pvalue")
        
        stopifnot(all(a$protein == b$protein))
        stopifnot(all(a$module == b$module))
        
        module_membership = bind_cols(a, b %>% select(pvalue)) # Eigengene-based connectivity, also known as module membership.
        
        # Correlate proteins and traits
        protein_trait_raw = corAndPvalue(
            x$datExpr, # samples by proteins
            x$datExpr %>% # samples by traits
                as_tibble(rownames = "sample") %>%
                select(sample) %>% 
                transmute(animal = str_sub(sample, 2, 3)) %>%
                left_join(
                    metadata_numeric,
                    by = "animal"
                )
            
        )
        
        a = protein_trait_raw$cor %>%
            as_tibble(rownames = "protein") %>%
            pivot_longer(-protein, names_to = "trait", values_to = "coefficient")
        
        b = protein_trait_raw$p %>%
            as_tibble(rownames = "protein") %>%
            pivot_longer(-protein, names_to = "trait", values_to = "pvalue")
        
        stopifnot(all(a$protein == b$protein))
        stopifnot(all(a$trait == b$trait))
        
        gene_significance = bind_cols(a, b %>% select(pvalue))
        
        list(
            group_index = x$group_index,
            module_membership = module_membership,
            gene_significance = gene_significance
            
        )
    }
)


# Plot each module that is significantly correlated to a trait.
mm_gs_filtered = lapply(
    trait_modules_of_interest %>% rowwise() %>% group_split(), # i = (trait_modules_of_interest %>% rowwise() %>% group_split)[[1]]
    function(i) { 
        
        message(paste(i, collapse = ", "))
        
        presentable_ = filter(groups, group_index == i$group_index)$presentable
        
        
        df_module = mm_gs[[i$group_index]]$module_membership %>%
            filter(module == i$module) %>%
            rename(
                module_membership = coefficient, # formerly module_connectivity
                module_pvalue = pvalue
            )
        
        
        df_trait = mm_gs[[i$group_index]]$gene_significance %>%
            filter(trait == i$trait) %>%
            rename(
                trait_correlation = coefficient, # formerly trait_correlation
                trait_pvalue = pvalue
            )
        
        mm_gs_joined = left_join(
            df_module,
            df_trait,
            by = "protein"
        ) %>%
            mutate(significance = case_when(
                module_pvalue < 0.05 & trait_pvalue < 0.05 ~ "module & trait",
                module_pvalue < 0.05 ~ "module",
                trait_pvalue < 0.05 ~ "trait",
                TRUE ~ "none"
                )
            )
        
        mm_gs_joined %>%
            ggplot(aes(module_membership, trait_correlation, color = significance)) + 
            geom_point() +
            geom_smooth(mapping = aes(group = 1), method = "lm") +
            labs(
                title = filter(groups, group_index == i$group_index)$presentable,
                subtitle = "Gene significance: Module connectivity and trait correlation.",
                x = paste(i$module, "(module membership)"),
                y = paste(i$trait, "(trait correlation)")
            )
        
        ggsave(generate_fig_name(output_module_membership_trait_significance_file, paste_("moduletrait", i$trait, presentable_, i$module)), height = height, width = width)
        
        # Save the significant ones to disk.
        mm_gs_joined %>%
            filter(significance == "module & trait") %>%
            arrange(module_membership, trait_correlation) %>%
            mutate(group_index = i$group_index) %>%
            relocate(group_index)
        
    }
    
) %>%
    bind_rows()


if (!interactive()) {
    mm_gs_filtered %>%
        write_rds_and_tsv(output_module_membership_trait_significance_file)
        
}

