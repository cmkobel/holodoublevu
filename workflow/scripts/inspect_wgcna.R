#!/usr/bin/env Rscript
# rm(list = ls())
library(tidyverse)
library(WGCNA)
source("workflow/scripts/utils.R")

# --- Inputs

metadata_file <- snakemake@input[["metadata"]] %>% as.character()
net_results_file <- snakemake@input[["wgcna_modules"]] %>% as.character()
groups_file = snakemake@input[["groups"]] %>% as.character()
output_rds_file <- snakemake@output[["mod2mod"]] %>% as.character()

# For debugging
if (F) {
    metadata_file <- "resources/metadata_v1.6.tsv"
    net_results_file <- "results/ig/luing/wgcna/modules.rds"
    groups_file <- "results/ig/luing/imputed/groups.tsv"
    output_rds_file <- "results/ig/luing/wgcna/inspected/module2module.rds"
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



# --- Basic dendro viz


i <- net_results[[1]]

lapply(
    net_results,
    function(i) {
        pdf(generate_fig_name(output_rds_file, "dendro"), height = height, width = width)
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
            # main = paste(i$title, "(block 1 only)"),
            main = "Dendro",
            dendroLabels = F
        )
        dev.off()

        # ggsave(generate_fig_name(output_rds_file, "dendro"))
    }
)




# --- pheatmaps
p <- lapply(
    net_results, # [1], # select one only with [n] where in is in 1:6.
    function(x) {
        pdf(generate_fig_name(output_rds_file, "pheatmap"), height = height, width = width)
        pheatmap::pheatmap(
            x$net$MEs,
            # main = paste0("Module Eigengenes\n", keys_presentable %>% filter(group_index == x$group_index) %>% pull(presentable) %>% unique())
            main = paste(
                groups %>% filter(group_index == i$group_index),
                collapse = ", "
            )
        )
        dev.off()
        # ggsave(generate_fig_name(output_rds_file, "pheatmap"))
    }
)



# --- ## Visualizing module eigengenes


lapply(
    net_results, # [1], # select one only with [n] where in is in 1:6.
    function(i) {
        pdf(generate_fig_name(output_rds_file, "me_pheatmap"), height = height, width = width)
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



wanted_keys = list(
    # Aberdeen Angus X
    list(A = layer_digesta, B = layer_wall), # digesta x wall
    list(A = layer_digesta, B = layer_liver), # digesta x liver
    list(A = layer_liver, B = layer_wall) # liver x wall
)


bicor_module2module <- lapply(
    wanted_keys,
    function(x) {
        pdf(generate_fig_name(output_rds_file, "mod2mod"), height = height, width = width)

        i <- x # i = list(A = 1, B = 7) # DEBUG

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
        message("this far")




        presentable_A <- paste(
            groups %>% filter(group_index == i$A),
            collapse = ", "
        )
        presentable_B <- paste(
            groups %>% filter(group_index == i$B),
            collapse = ", "
        )



        message("3")
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

net_results[groups %>%
    filter(collection != "tube") %>%
    pull(group_index)] %>%
    lapply(
        function(x) {
            pdf(generate_fig_name(output_rds_file, "pheno"), height = height, width = width)

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
