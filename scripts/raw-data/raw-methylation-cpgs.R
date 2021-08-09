get_promoter_cpgs_meth_no_filt <- function(recalc = FALSE) {
    cov_file <- here("data/promoter_cpgs_cov_no_filt.csv")
    avg_file <- here("data/promoter_cpgs_avg_meth_no_filt.csv")    
    if (!file.exists(avg_file) || !file.exists(cov_file) || recalc) {
        meth_list <- gextract_meth_parallel(samp_data$track, samp_data$samp, intervals = promoter_intervs, iterator = "intervs.global.seq_CG", min_cov = 1)
        fwrite(meth_list$cov, cov_file)
        fwrite(meth_list$avg, avg_file)
    } else {
        meth_list <- list(avg = as_tibble(fread(avg_file)), cov = as_tibble(fread(cov_file)))
    }
    invisible(meth_list)
}


get_promoter_cpgs_meth <- function(recalc = FALSE, min_cov = 20) {
    cov_file <- main_config$promoter_cpgs_methylation$cov_file
    avg_file <- main_config$promoter_cpgs_methylation$avg_meth_file
    if (!file.exists(avg_file) || !file.exists(cov_file) || recalc) {
        meth_list <- gextract_meth_parallel(samp_data$track, samp_data$samp, intervals = promoter_intervs, iterator = "intervs.global.seq_CG", min_cov = min_cov)
        fwrite(meth_list$cov, cov_file)
        fwrite(meth_list$avg, avg_file)
    } else {
        meth_list <- list(avg = as_tibble(fread(avg_file)), cov = as_tibble(fread(cov_file)))
    }
    invisible(meth_list)
}

get_promoter_cpgs_avg_meth <- function(add_names = FALSE, normal_fraction = 0.7, tumor_fraction = 0.7, recalc = FALSE) {
    avg_meth_file <- main_config$promoter_cpgs_methylation$filtered_avg_meth_file

    if (!file.exists(avg_meth_file) || recalc) {
        promoter_meth <- get_promoter_cpgs_meth()
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

get_genomic_msp1_cpgs_meth_no_filt <- function(recalc = FALSE) {
    cov_file <- here("data/genomic_msp1_cpgs_cov_no_filt.csv")
    avg_file <- here("data/genomic_msp1_cpgs_avg_no_filt.csv")
    if (!file.exists(avg_file) || !file.exists(cov_file) || recalc) {
        intervs <- get_genomic_avg_meth() %>% select(chrom, start, end)
        intervs <- gintervals.intersect("intervs.global.seq_CG", intervs)
        meth_list <- gextract_meth_parallel(samp_data$track, samp_data$samp, intervals = intervs, iterator = "intervs.global.seq_CG", min_cov = 1)
        fwrite(meth_list$cov, cov_file)
        fwrite(meth_list$avg, avg_file)
    } else {
        meth_list <- list(avg = as_tibble(fread(avg_file)), cov = as_tibble(fread(cov_file)))
    }
    invisible(meth_list)
}

get_msp1_cpgs_meth <- function(recalc = FALSE, min_cov = 20) {
    cov_file <- main_config$msp1_cpgs_methylation$cov_file
    avg_file <- main_config$msp1_cpgs_methylation$avg_meth_file
    if (!file.exists(avg_file) || !file.exists(cov_file) || recalc) {
        meth_list <- gextract_meth_parallel(samp_data$track, samp_data$samp, intervals = "intervs.msp1.fid", iterator = "intervs.global.seq_CG", min_cov = min_cov)
        fwrite(meth_list$cov, cov_file)
        fwrite(meth_list$avg, avg_file)
    } else {
        meth_list <- list(avg = as_tibble(fread(avg_file)), cov = as_tibble(fread(cov_file)))
    }
    invisible(meth_list)
}

get_genomic_cpgs_avg_meth <- function(normal_fraction = 0.7, tumor_fraction = 0.7, recalc = FALSE) {
    cov_file <- main_config$genomic_msp1_cpgs_methylation$cov_file
    avg_meth_file <- main_config$genomic_msp1_cpgs_methylation$avg_meth_file

    if (!file.exists(cov_file) || !file.exists(avg_meth_file) || recalc) {
        msp1_meth <- get_msp1_cpgs_meth()
        msp1_genomic_intervs <- gintervals.intersect("intervs.global.seq_CG", "intervs.msp1.fid") %>%
            gintervals.neighbors1(promoter_intervs) %>%
            filter(abs(dist) > 0) %>%
            select(chrom, start, end) %>%
            gintervals.neighbors1("intervs.global.exon") %>%
            filter(abs(dist) > 10) %>%
            select(chrom, start, end)
            
        genomic_meth <- msp1_meth$avg %>% inner_join(msp1_genomic_intervs, by = c("chrom", "start", "end"))
        genomic_meth <- filter_avg_meth(genomic_meth, normal_fraction = normal_fraction, tumor_fraction = tumor_fraction)              
        genomic_meth_cov <- msp1_meth$cov %>% inner_join(genomic_meth %>% select(chrom, start, end), by = c("chrom", "start", "end"))
        
        fwrite(genomic_meth_cov, main_config$genomic_msp1_cpgs_methylation$cov_file)
        fwrite(genomic_meth, main_config$genomic_msp1_cpgs_methylation$avg_meth_file)
    } else {
        genomic_meth <- fread(main_config$genomic_msp1_cpgs_methylation$avg_meth_file)
    }

    return(genomic_meth)
}


export_samp_raw_data <- function(samp, track, hg19_ofn, hg38_ofn){
    message(samp)
    gsetroot(here("db/trackdb"))
    data <- gextract(
            paste0(track, c(".cov", ".meth")), 
            iterator=glue("{track}.cov"), 
            intervals=glue("{track}.cov"), 
            colnames=c("cov", "meth")) %>% 
        select(-intervalID) %>% 
        arrange(chrom, start, end) %>% 
        mutate(avg = meth / cov)

    fwrite(data, hg19_ofn, sep="\t", quote=FALSE)

    # liftover to hg38
    on.exit(gsetroot(here("db/trackdb")))
    gsetroot("/home/aviezerl/hg38")
    hg38_data <- data %>% 
        mutate(intervalID = 1:n()) %>% 
        select(intervalID, cov, meth, avg) %>% 
        inner_join(gintervals.liftover(intervals = data, chain = here("data/hg19ToHg38.over.chain")), by = "intervalID") %>% 
        select(chrom, start, end, cov, meth, avg) %>% 
        arrange(chrom, start, end)

    fwrite(hg38_data, hg38_ofn, sep="\t", quote=FALSE)
}


export_all_raw_data <- function(){
    dir.create(here("data/raw_exported/hg19"), showWarnings = FALSE)
    dir.create(here("data/raw_exported/hg38"), showWarnings = FALSE)
    doMC::registerDoMC(40)
    plyr::a_ply(samp_data, 1, function(x) export_samp_raw_data(
        x$samp, 
        x$track, 
        glue(here("data/raw_exported/hg19/{x$samp}.tsv.gz")),
        glue(here("data/raw_exported/hg38/{x$samp}.tsv.gz"))
        ), 
        .parallel=TRUE)
}