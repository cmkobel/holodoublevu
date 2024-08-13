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


handful <- function(df, .return = F, n = 7) {
    # Shows a random handful op items. For safety it doesn't return anything (by default).

    message("Showing ", n, " random rows:")

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
