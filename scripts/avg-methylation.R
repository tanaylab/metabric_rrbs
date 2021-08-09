get_all_meth <- function(){
    genomic_avg_meth <- fread(here("data/genomic_msp1_avg_meth_filt.csv")) %>% as_tibble()
    prom_avg_meth <- fread(here("data/promoter_avg_meth_filt.csv")) %>% 
        select(-name, -name3.chr) %>% 
        distinct(chrom, start, end, .keep_all = TRUE) %>% 
        as_tibble()

    all_mat_raw <- bind_rows(
        prom_avg_meth, 
        genomic_avg_meth
    )

    return(all_mat_raw)
}

# Get average of ER+/ER-/normal samples of each locus
get_all_summary_meth <- function(){
    get_all_meth() %>% 
        mutate(
            `ER+` = rowMeans(select(., any_of(ER_positive_samples)), na.rm = TRUE), 
            `ER-` = rowMeans(select(., any_of(ER_negative_samples)), na.rm = TRUE), 
            normal = rowMeans(select(., any_of(normal_samples)), na.rm = TRUE)
        ) %>% 
        select(chrom, start, end, `ER+`, `ER-`, normal) %cache_df% 
        here("data/all_meth_summary.tsv") %>% 
        as_tibble()
}


get_all_summary_meth_sd <- function(){
    get_all_meth() %>% 
        mutate(
            `ER+` = matrixStats::rowSds(as.matrix(select(., any_of(ER_positive_samples))), na.rm = TRUE), 
            `ER-` = matrixStats::rowSds(as.matrix(select(., any_of(ER_negative_samples))), na.rm = TRUE), 
            normal = matrixStats::rowSds(as.matrix(select(., any_of(normal_samples))), na.rm = TRUE)
        ) %>% 
        select(chrom, start, end, `ER+`, `ER-`, normal) %cache_df% 
        here("data/all_meth_summary_sd.tsv") %>% 
        as_tibble()
}