#!/usr/bin/env Rscript

library(tidyverse)
library(WGCNA)
source("workflow/scripts/utils.R")

# --- Inputs

metadata_file <- snakemake@input[["metadata"]] %>% as.character()
proteome_intensities_file <- snakemake@input[["imputed"]] %>% as.character()
groups_file <- snakemake@input[["groups"]] %>% as.character()
samples <- snakemake@params[["samples"]]

output_wgcna_modules_file = snakemake@output[["wgcna_modules"]] %>% as.character()


# For debugging
if (F) {
    metadata_file <- "resources/metadata_v1.6.tsv"
    proteome_intensities_file <- "results/ig/luing/imputed/proteome_intensities.rds"
    groups_file <- "results/ig/luing/imputed/groups.tsv"
    samples <- c("D06T6S", "D08T6S", "D16T6S", "D25T6S", "D28T1S", "D28T2S", "D28T3S", "D28T4S", "D28T5S", "D28T6S", "D29T6S", "D35T6S", "D51T1S", "D51T2S", "D51T3S", "D51T4S", "D51T5S", "D51T6S", "D60T6S", "D61T6S", "D76T6S", "L06T6R", "L08T6R", "L16T6R", "L25T6R", "L28T6R", "L29T6R", "L35T6R", "L42T6R", "L51T6R", "L60T6R", "L61T6R", "L76T6R", "W06T6R", "W08T6R", "W16T6R", "W25T6R", "W29T6R", "W35T6R", "W42T6R", "W51T6R", "W60T6R", "W61T6R", "W76T6R") # luing example

    output_wgcna_modules_file <- "results/ig/{imputation_group}/wgcna/modules.rds"
}


metadata <- read_tsv(metadata_file)
proteome_intensities <- read_rds(proteome_intensities_file)
groups <- read_tsv(groups_file)




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

        current_group <- groups %>% filter(group_index == current_group_index)
        current_group

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
        pdf(generate_fig_name(output_wgcna_modules_file, paste_("hclust", filter(groups, group_index == i$group_index)$presentable)))

        plotDendroAndColors(
            i$sample_tree,
            i$trait_colors_data %>% numbers2colors(),
            groupLabels = names(i$trait_colors_data),
            main = filter(groups, group_index == i$group_index)$presentable
        )

        dev.off()
    }
)


# --- Pick scalefree threshold
# These will be added to the groups table



i <- (proteome_intensities %>% # Debug friendly setup.
    group_by(group_index) %>%
    group_split())[[1]]

thresholds <- lapply(
    (proteome_intensities %>%
        group_by(group_index) %>%
        group_split()),
    function(i) {
        i %>% handful()

        current_group_index <- i$group_index[[1]]
        message("group index ", current_group_index)

        current_group <- groups %>% filter(group_index == current_group_index)
        current_group

        wide <- i %>%
            pivot_wider(id_cols = c(sample, group_index), names_from = "protein", values_from = "intensity") %>%
            select(-group_index)

        pst <- pickSoftThreshold(
            wide %>% column_to_rownames("sample"),
            RsquaredCut = .8, # Default is .85. Maybe set it to .8
            powerVector = c(seq(1, 10, by = 2), seq(12, 30, by = 3))
        )


        list(
            group_index = current_group_index,
            pst = pst,
            wide = wide
        )
    }
)


# pst plots
pdf(generate_fig_name(output_wgcna_modules_file, "pst"))

lapply(
    thresholds,
    function(i) {
        i$pst$fitIndices %>%
            pivot_longer(-Power) %>%
            ggplot(aes(Power, value)) +
            geom_point() +
            geom_line() +
            facet_wrap(~name, scales = "free") +
            geom_vline(xintercept = i$pst$powerEstimate, linetype = "dashed", color = "red") +
            geom_vline(xintercept = 12, linetype = "dashed", color = "green3", alpha = 0.7) +
            geom_hline(yintercept = .8, linetype = "dashed", color = "grey", alpha = 0.5) +
            theme() +
            labs(
                title = filter(groups, group_index == i$group_index)$presentable,
                subtitle = "WGCNA::pickSoftThreshold",
                caption = paste0("powerEstimate = ", i$pst$powerEstimate)
            )
    }
)

dev.off()

groups_estimate <- groups %>%
    rowwise() %>%
    mutate(
        threshold_estimate = (function(x) {
            thresholds[[x]]$pst$powerEstimate
        })(group_index)
    )


# Since I really want to be able to compare across layers, I'm gonna pick the same threshold.
groups_final = groups_estimate %>%
    mutate(threshold_final = 12)


# --- Make modules



i <- (proteome_intensities %>% # Debug friendly setup.
    group_by(group_index) %>%
    group_split())[[1]]

net_results <- lapply(
    (proteome_intensities %>%
        group_by(group_index) %>%
        group_split()),
    function(i) {
        i %>% handful()

        current_group_index <- i$group_index[[1]]
        message("group index ", current_group_index)

        current_group <- groups_final %>% filter(group_index == current_group_index)
        current_group

        wide <- i %>%
            pivot_wider(id_cols = c(sample, group_index), names_from = "protein", values_from = "intensity") %>%
            select(-group_index)


        current_power <- current_group$threshold_final
        message("current power ", current_power)

        # Should probably try this with corType = "bicor" instead of the default corType = "pearson"
        net <- blockwiseModules(
            datExpr = wide %>% column_to_rownames("sample"),
            power = current_power,
            networkType = "signed",
            TOMType = "signed",
            # corType = "bicor", # Wanxin's idea.
            numericLabels = T # T if you want to manually recolor.
            # reassignThreshold = 0 # Torgeir's suggestion
        )

        list(
            group_index = current_group_index,
            net = net

            # numeric_net = numeric_net,
            # new_words = new_words
        )
    }
)


# Save net results to disk
net_results %>%
    write_rds(snakemake@output[["wgcna_modules"]] %>% as.character())
