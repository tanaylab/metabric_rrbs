## Screen for coverage ###
gscreen_coverage <- function(tracks, min_cov, min_samples) {
    add_covs_chrom <- function(tracks, chr, min_cov, min_samples) {
        message(glue("starting {chr}"))
        chrom_len <- gintervals.all() %>%
            filter(chrom == chr) %>%
            .$end
        poss <- rep(0, chrom_len)
        i <- 1
        for (track in tracks) {
            if ((i / length(tracks) * 100) %% 5 == 0) {
                message(sprintf("%s: %s%%", chr, i / length(tracks) * 100))
            }

            i <- i + 1
            a <- gextract(paste0(track, ".cov"), intervals = gintervals.all() %>% filter(chrom == chr), colnames = "cov") %>% tbl_df()
            poss[a$start] <- poss[a$start] + (a$cov >= min_cov)
        }
        covs <- tibble(start = which(poss != 0), chrom = chr) %>%
            mutate(samples = poss[start]) %>%
            select(chrom, start, samples) %>%
            filter(samples >= min_samples)
        message(glue("finished {chr}"))
        return(covs)
    }

    res <- as.character(gintervals.all()$chrom) %>% plyr::alply(1, function(x) add_covs_chrom(tracks, x, min_cov = min_cov, min_samples = min_samples), .parallel = TRUE)
    covs <- res %>%
        map_df(~.x) %>%
        mutate(end = start + 1) %>%
        select(chrom, start, end, samples)
    return(covs)
}

## Marginal coverage ###
covs_marginal <- function(tracks) {
    add_covs_chrom <- function(tracks, chr) {
        message(glue("starting {chr}"))
        chrom_len <- gintervals.all() %>%
            filter(chrom == chr) %>%
            pull(end)
        poss <- rep(0, chrom_len)
        for (track in tracks) {
            tr_data <- gextract(paste0(track, ".cov"), intervals = gintervals.all() %>% filter(chrom == chr), colnames = "cov") %>% as_tibble()
            poss[tr_data$start] <- poss[tr_data$start] + tr_data$cov
        }
        message(glue("finished {chr}"))
        tibble(start = which(poss != 0), chrom = chr) %>%
            mutate(cov = poss[start]) %>%
            select(chrom, start, cov)
    }

    res <- as.character(gintervals.all()$chrom) %>% plyr::alply(1, function(x) add_covs_chrom(tracks, x), .parallel = TRUE)

    covs <- res %>%
        map_df(~.x) %>%
        mutate(end = start + 1) %>%
        select(chrom, start, end, cov)
    return(covs)
}
