get_pat_stats <- function(recalc = FALSE) {
    fn <- main_config$data_files$pats_stats

    if (!file.exists(fn) || recalc) {
        pats_d <- plyr::adply(samp_data, 1, function(x) {
            print(x$samp)
            return(gpatterns.extract_all(x$track, samples = x$samp, dsn = 30, na.rm = TRUE))
        }, .parallel = TRUE)
        pats_d <- pats_d %>%
            select(chrom, start, end, samp, fid, ncpg:epipoly) %>%
            as_tibble()
        fwrite(pats_d, fn)
    } else {
        pats_d <- fread(fn) %>% as_tibble()
    }

    return(pats_d)
}

get_intervals_epipoly <- function(intervals, max_dist = 0) {
    pat_stats <- get_pat_stats()
    intervals_fids <- intervals %>%
        gintervals.neighbors1(pat_stats %>% distinct(chrom, start, end, fid)) %>%
        filter(abs(dist) <= max_dist) %>%
        select(-(chrom1:end1))
    intervals_epi <- intervals_fids %>% left_join(pat_stats %>% select(samp, fid, pat_meth, epipoly))
    return(intervals_epi)
}