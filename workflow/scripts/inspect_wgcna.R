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


# For debugging
if (interactive()) {
    metadata_file <- "resources/metadata_v1.6.tsv"
    net_results_file <- "results/ig/both/wgcna/modules.rds"
    groups_file <- "results/ig/both/imputed/groups.tsv"
    output_rds_file <- "results/ig/both/wgcna/inspected/module2module.rds"
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


        list(
            "A" = i$A,
            "B" = i$B,
            "presentable_A" = presentable_A,
            "presentable_B" = presentable_B,
            "bicor" = bicor_results
        )
        dev.off()
    }
)


bicor_module2module %>%
    write_rds(output_rds_file)


# --- Correlating to phenotypes




# By converting this before the for loop of each layer, I know that the factors will be comparable between layers. Makes visual interpretation easier.
metadata_numeric <- metadata %>%
    mutate(across(
        # c(everything(), -animal),
        -animal,
        ~ .x %>%
            as.factor() %>%
            as.numeric()
    ))






# You will get problems with duplicate keys if you run this with the tube samples.


lapply(
    net_results[groups %>%
        filter(collection != "tube") %>%
        pull(group_index)],
    function(x) {
        pdf(generate_fig_name(output_rds_file, paste_("pheno", filter(groups, group_index == x$group_index)$presentable)), height = height, width = width)

        presentable <- paste(
            groups %>% filter(group_index == x$group_index),
            collapse = ", "
        )

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
            main = presentable
        ) %>%
            print()

        dev.off()
    }
)
