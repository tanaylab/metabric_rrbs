get_mut_df <- function() {
    return(fread(here("data/mutations.tsv")))
}
