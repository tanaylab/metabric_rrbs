get_msp1_meth <- function(recalc = FALSE, min_cov = 20) {
    cov_file <- here(main_config$msp1_methylation$cov_file)
    avg_file <- here(main_config$msp1_methylation$avg_meth_file)
    if (!file.exists(avg_file) || !file.exists(cov_file) || recalc) {
        meth_list <- gextract_meth_parallel(samp_data$track, samp_data$samp, intervals = "intervs.msp1.fid", iterator = "intervs.msp1.fid", min_cov = min_cov)
        fwrite(meth_list$cov, cov_file)
        fwrite(meth_list$avg, avg_file)
    } else {
        meth_list <- list(avg = as_tibble(fread(avg_file)), cov = as_tibble(fread(cov_file)))
    }
    invisible(meth_list)
}

get_promoter_meth <- function(recalc = FALSE, min_cov = 20) {
    cov_file <- here(main_config$promoter_methylation$cov_file)
    avg_file <- here(main_config$promoter_methylation$avg_meth_file)
    if (!file.exists(avg_file) || !file.exists(cov_file) || recalc) {
        meth_list <- gextract_meth_parallel(samp_data$track, samp_data$samp, intervals = promoter_intervs, iterator = promoter_intervs, min_cov = min_cov)
        fwrite(meth_list$cov, cov_file)
        fwrite(meth_list$avg, avg_file)
    } else {
        meth_list <- list(avg = as_tibble(fread(avg_file)), cov = as_tibble(fread(cov_file)))
    }
    invisible(meth_list)
}

get_promoter_avg_meth <- function(add_names = FALSE, normal_fraction = 0.7, tumor_fraction = 0.7, recalc = FALSE) {
    avg_meth_file <- here(main_config$promoter_methylation$filtered_avg_meth_file)

    if (!file.exists(avg_meth_file) || recalc) {
        promoter_meth <- get_promoter_meth()
        promoter_avg_meth <- filter_avg_meth(promoter_meth$avg, normal_fraction = normal_fraction, tumor_fraction = tumor_fraction)

        promoter_avg_meth <- promoter_avg_meth %>%
            left_join(promoter_intervs %>% select(chrom, start, end, name, name3.chr), by = c("chrom", "start", "end")) %>%
            select(chrom, start, end, name, name3.chr, everything())
        fwrite(promoter_avg_meth, avg_meth_file)
    } else {
        promoter_avg_meth <- fread(avg_meth_file)
    }

    if (!add_names) {
        promoter_avg_meth <- promoter_avg_meth %>% select(-name, -name3.chr) %>% distinct(chrom, start, end, .keep_all = TRUE)
    }
    return(promoter_avg_meth)
}

get_genomic_avg_meth <- function(normal_fraction = 0.7, tumor_fraction = 0.7, recalc = FALSE) {
    cov_file <- here(main_config$genomic_msp1_methylation$cov_file)
    avg_meth_file <- here(main_config$genomic_msp1_methylation$avg_meth_file)

    if (!file.exists(avg_meth_file) || recalc) {
        msp1_meth <- get_msp1_meth()
        msp1_genomic_intervs <- gintervals.load("intervs.msp1.fid") %>%
            gintervals.neighbors1(promoter_intervs) %>%
            filter(abs(dist) > 0) %>%
            select(chrom, start, end) %>%
            gintervals.neighbors1("intervs.global.exon") %>%
            filter(abs(dist) > 10) %>%
            select(chrom, start, end)
            
        genomic_meth <- msp1_meth$avg %>% inner_join(msp1_genomic_intervs, by = c("chrom", "start", "end"))
        genomic_meth <- filter_avg_meth(genomic_meth, normal_fraction = normal_fraction, tumor_fraction = tumor_fraction)              
        genomic_meth_cov <- msp1_meth$cov %>% inner_join(genomic_meth %>% select(chrom, start, end), by = c("chrom", "start", "end"))
        
        fwrite(genomic_meth_cov, cov_file)
        fwrite(genomic_meth, avg_meth_file)
    } else {
        genomic_meth <- fread(avg_meth_file)
    }

    return(genomic_meth)
}

get_genomic_avg_meth_no_filt <- function(recalc = FALSE) {
    cov_file <- here("data/genomic_msp1_cov_no_filt.csv")
    avg_file <- here("data/genomic_msp1_avg_no_filt.csv")
    if (!file.exists(avg_file) || !file.exists(cov_file) || recalc) {
        msp1_meth <- get_msp1_meth()
        msp1_genomic_intervs <- gintervals.load("intervs.msp1.fid") %>%
            gintervals.neighbors1(promoter_intervs) %>%
            filter(abs(dist) > 0) %>%
            select(chrom, start, end) %>%
            gintervals.neighbors1("intervs.global.exon") %>%
            filter(abs(dist) > 10) %>%
            select(chrom, start, end)
            
        genomic_meth <- msp1_meth$avg %>% inner_join(msp1_genomic_intervs, by = c("chrom", "start", "end"))
        genomic_meth_cov <- msp1_meth$cov %>% inner_join(msp1_genomic_intervs %>% select(chrom, start, end), by = c("chrom", "start", "end"))

        fwrite(genomic_meth_cov, cov_file)
        fwrite(genomic_meth, avg_file)
    } else {
        meth_list <- list(avg = as_tibble(fread(avg_file)), cov = as_tibble(fread(cov_file)))
    }
    invisible(meth_list)
}

get_exon_meth <- function(recalc = FALSE, min_cov = 20) {
    cov_file <- here(main_config$exon_methylation$cov_file)
    avg_file <- here(main_config$exon_methylation$avg_meth_file)
    if (!file.exists(avg_file) || !file.exists(cov_file) || recalc) {
        meth_list <- gextract_meth_parallel(samp_data$track, samp_data$samp, intervals = exon_intervs, iterator = exon_intervs, min_cov = min_cov)
        fwrite(meth_list$cov, cov_file)
        fwrite(meth_list$avg, avg_file)
    } else {
        meth_list <- list(avg = as_tibble(fread(avg_file)), cov = as_tibble(fread(cov_file)))
    }
    invisible(meth_list)
}

get_exon_avg_meth <- function(add_names = FALSE, normal_fraction = 0.7, tumor_fraction = 0.7, recalc = FALSE) {
    avg_meth_file <- here(main_config$exon_methylation$filtered_avg_meth_file)

    if (!file.exists(avg_meth_file) || recalc) {
        exon_meth <- get_exon_meth()
        exon_avg_meth <- filter_avg_meth(exon_meth$avg, normal_fraction = normal_fraction, tumor_fraction = tumor_fraction)

        exon_avg_meth <- exon_avg_meth %>%
            left_join(exon_intervs %>% select(chrom, start, end, name), by = c("chrom", "start", "end")) %>%
            select(chrom, start, end, name, everything())
        fwrite(exon_avg_meth, avg_meth_file)
    } else {
        exon_avg_meth <- fread(avg_meth_file)
    }

    if (!add_names) {
        exon_avg_meth <- exon_avg_meth %>% select(-name) %>% distinct(chrom, start, end, .keep_all = TRUE)
    }
    return(exon_avg_meth)
}

get_all_meth <- function(){
    bind_rows(
        get_promoter_avg_meth(), 
        get_genomic_avg_meth()
    )
}

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

get_promoter_summary_stats <- function(recalc = FALSE) {
    promoter_stats_fn <- here(main_config$promoter_methylation$stats_file)

    if (!file.exists(promoter_stats_fn) || recalc) {
        promoter_meth <- get_promoter_avg_meth()
        promoter_stats <- loci_meth_summary_stats(promoter_meth)
        fwrite(promoter_stats, promoter_stats_fn)
    } else {
        promoter_stats <- fread(promoter_stats_fn)
    }

    return(promoter_stats)
}


