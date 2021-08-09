get_samples_coverage_dist <- function(tracks = samp_data$track, breaks = c(0:100, 1e5), intervals = gintervals.all()) {    

    # Calculate CpG coverage distribution per sample
    samp_covs <- plyr::alply(tracks, 1, function(x) gdist(glue("{x}.cov"), breaks, include.lowest = FALSE, intervals = intervals) %>% tibble::enframe(), .parallel = TRUE)

    samp_covs <- map2_dfr(samp_covs, samp_data$samp, ~ .x %>%
        rename(cov = name, cpgs = value) %>%
        mutate(samp = .y))

    n_cpgs <- nrow(gintervals.intersect("intervs.global.seq_CG", intervals))

    zero_covs <- samp_covs %>%
        group_by(samp) %>%
        summarise(s_cpgs = sum(cpgs)) %>%
        mutate(cpgs = n_cpgs - s_cpgs) %>%
        mutate(cov = "(-1,0]") %>%
        select(cov, cpgs, samp)

    samp_covs <- bind_rows(samp_covs, zero_covs) %>% arrange(samp)

    samp_covs <- samp_covs %>% 
        mutate(cov = gsub("\\(-*\\d+,", "", cov), cov = gsub("]", "", cov), cov = as.numeric(cov))

    return(samp_covs)
}

get_promoter_expr_covs <- function(){
    mat_expr <- get_gene_expression_mat()
    mat_expr <- mat_expr %>%
            select(any_of(c("chrom", "start", "end", "name", "name3.chr", samp_data$samp))) 
    prom_expr <- mat_expr %>%
            mutate(
                max_expr = matrixStats::rowMaxs(as.matrix(mat_expr[, -1:-5]), na.rm = TRUE),
                n_expr = rowSums(!is.na(as.matrix(mat_expr[, -1:-5])), na.rm = TRUE),
            ) %>%
            select(chrom, start, end, name, name3.chr, max_expr, n_expr)

    promoter_covs <- fread(here(main_config$promoter_methylation$cov_file)) %>%
            select(any_of(c("chrom", "start", "end", paste0(samp_data$samp, ".cov")))) %>%
            mutate(cov = rowMeans(.[, -1:-3], na.rm = TRUE)) %>%
            select(chrom, start, end, cov) %>% 
            distinct()

    prom_expr_cov <- promoter_covs %>% inner_join(prom_expr, by = c("chrom", "start", "end"))

    return(prom_expr_cov)
}

get_sample_tot_meth_calls <- function(intervals = gintervals.all()) {
    samp_covs <- plyr::adply(samp_data %>% select(samp, track), 1, function(x) {
        message(x$samp)
        enframe(gsummary(glue("{x$track}.cov"), intervals))
    }, .parallel = TRUE)

    samp_covs <- samp_covs %>%
        select(-track) %>%
        spread(name, value) %>%
        as_tibble()

    return(samp_covs)

}
