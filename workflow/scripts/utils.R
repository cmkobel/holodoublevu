supacow_separate_sample <- function(df, keep_debug_columns = FALSE, keep_sample = FALSE) {
    # separates sample column into source-animal-timepoint and solution columns.
    rv <- df %>%
        separate_wider_regex(
            sample,
            patterns = c(
                source = "^[A-Z]",
                animal = "\\d{2}",
                timepoint = "T\\d",
                solution = "[A-Z]"
            ),
            too_few = "debug"
        )

    if (keep_sample) {
        rv <- rv
    } else {
        rv <- select(rv, -sample)
    }

    if (keep_debug_columns) {
        rv <- rv
    } else {
        rv <- select(rv, -sample_ok, -sample_matches, -sample_remainder)
    }
}

supacow_paste_sample <- function(df, keep_columns = FALSE) {
    rv <- df %>%
        mutate(sample = paste0(source, animal, timepoint, solution))

    if (keep_columns) {
        rv <- rv
    } else {
        rv <- select(rv, -source, -animal, -timepoint, -solution)
    }

    rv
}


handful <- function(df, .return = F, n = 7, seed = NULL) {
    # Shows a random handful op items. For safety it doesn't return anything (by default).

    message("Showing ", n, " random rows:")

    if (!is.null(seed)) {
        message("Using specific seed for reproducibility.")
        set.seed(seed)
    }

    if ("data.frame" %in% class(df)) {
        slice_sample(df, n = n) %>%
            print()

        message(nrow(df), " rows in total.")
    } else {
        sample(df, n) %>%
            print()
        message(length(df), " items in total.")
    }

    if (.return) {
        df
    } else {
        NULL # This is for safety, so noone leaves a handful() call in the middle of production code while debugging.
    }
}

add_collection_column <- function(df) {
    df %>%
        mutate(collection = case_when(
            timepoint == "T6" ~ "slaughter",
            TRUE ~ "tube"
        ))
}



write_rds_and_tsv <- function(x, path_to_rds) {
    nrows <- nrow(x)
    ncols <- ncol(x)

    if (!str_detect(path_to_rds, "\\.rds$")) {
        stop("path_to_rds must end on \".rds\" but is \"", path_to_rds, "\"")
    }

    message("Writing ", nrows, " x ", ncols, " (rows x cols) to ", path_to_rds)
    write_rds(x, path_to_rds)

    path_to_tsv <- str_replace(path_to_rds, "\\.rds$", ".tsv")

    message("Writing ", nrows, " x ", ncols, " (rows x cols) to ", path_to_tsv)
    write_tsv(x, path_to_tsv)
}




# A function that returns a new global figure enumerator for each call

figno_var <- 0
figno <- function() {
    figno_var <<- figno_var + 1
    figno_var
}


generate_fig_name <- function(adjacent_file, title = "untitled", extension = "pdf") {
    rv <- paste0(dirname(adjacent_file), "/", title, "_fig", figno(), ".", extension)

    message("Plotting to file ", rv)
    dir.create(dirname(rv), showWarnings = F)

    rv
}


# height <- 12
# width <- 12

# pdf(generate_fig_name(output_file, paste("hclust", pastecomma(filter( groups, group_index == i$group_index)))))
# ggsave(generate_fig_name(output_pathway_enrichment_file, paste("title", pastecomma(filter(groups, group_index == i$group_index))), height = (height_multiplier / 7) + 2, width = width, create.dir = T))



pastecomma <- function(...) {
    paste(
        ...,
        collapse = ", ",
        sep = "_"
    )
}

paste_ <- function(...) {
    paste(
        ...,
        collapse = "_",
        sep = "_"
    )
}
