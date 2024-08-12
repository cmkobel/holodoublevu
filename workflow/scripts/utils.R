supacow_separate_sample <- function(df, keep_debug_columns = FALSE) {
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
        ) %>%
        select(-sample)

    if (keep_debug_columns) {
        rv
    } else {
        rv %>%
            select(-sample_ok, -sample_matches, -sample_remainder)
    }
}


random_head <- function(df, n = 10) {
    # Like head, but with random rows.
    df %>%
        slice_sample(n = n)
}
