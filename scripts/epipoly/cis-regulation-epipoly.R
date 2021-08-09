get_promoter_cis_reg_epipoly <- function() {       
    {
        cands <- fread(here("data/promoter_cis_cands.tsv")) %>% as_tibble()

        loci <- cands %>% group_by(chrom, start, end, ER, name) %>% summarise(type = ifelse(any(r == 1), "reg", "bg")) %>% ungroup()
        meth_clust_epi <- get_intervals_epipoly(loci %>% select(chrom:end, everything()), max_dist = 0)
        prom_epipoly <- meth_clust_epi %>% left_join(samp_data %>% select(samp, ER_samp = ER1)) %>% filter(ER == ER_samp)

        prom_epipoly
    } %cache_df% here("data/promoter_cis_reg_epipoly.tsv") %>% as_tibble()
}

get_genomic_cis_reg_epipoly <- function() {    
    {
        cands <- bind_rows(
            fread(here("data/genomic_cis_cands_ER_positive.tsv")),
            fread(here("data/genomic_cis_cands_ER_negative.tsv")),
            fread(here("data/genomic_cis_cands_normal.tsv"))
        )

        loci <- cands %>%
            filter(type == "obs", !is.na(dist), abs(dist) >= 2000) %>%
            mutate(type = ifelse(rank == 1 & abs(dist) <= 1e6, "reg", "bg")) %>%
            group_by(chrom, start, end, ER) %>%
            summarise(type = ifelse(any(type == "reg"), "reg", "bg")) %>%
            ungroup()

        meth_clust_epi <- get_intervals_epipoly(loci %>% select(chrom:end, everything()), max_dist = 0)
        epipoly <- meth_clust_epi %>% left_join(samp_data %>% select(samp, ER_samp = ER1)) %>% filter(ER == ER_samp)
        
     } %cache_df% here("data/genomic_cis_reg_epipoly.tsv") %>% as_tibble()
}