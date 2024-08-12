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


show_some <- function(df, n = 10) {
    # Like head, but with random rows.
    df %>%
        slice_sample(n = n) %>%
        print()

    df
}

add_collection_column <- function(df) {
    df %>%
        mutate(collection = case_when(
            timepoint == "T6" ~ "slaughter",
            TRUE ~ "tube"
        ))
}
