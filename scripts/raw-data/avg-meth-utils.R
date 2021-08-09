
gextract_meth_parallel <- function(tracks, names, intervals = gintervals.all(),
    iterator = "intervs.global.seq_CG",  min_cov  = NULL, threads = 30, io_saturation = TRUE, ...){

    smd <- tibble(track = tracks, name = names) 
    nchunks <- round(nrow(smd) / threads)
    smd <- smd %>%
        mutate(id = 1:n()) %>%
        mutate(chunk = ntile(id, nchunks)) %>%
        select(-id)    

    smd <- smd %>% 
        group_by(chunk) %>% 
        summarise(cmd = 
            glue("gextract_meth(tracks = c({trs}), names = c({nms}), intervals = intervals, iterator = iterator, min_cov = min_cov) %>% arrange(intervalID) %>% select(-intervalID)", 
                trs = paste(glue("'{track}'"), collapse=","),
                nms = paste(glue("'{name}'"), collapse=",")
                ))

    
    opt <- options(gmultitasking = FALSE)
    on.exit(options(opt))
    res <- misha.ext::gcluster.run2(command_list = smd$cmd, io_saturation = io_saturation, threads = threads)
    
    avg_tab <- map_dfc(res, ~ .x$retv %>% select(-ends_with(".cov"), -(chrom:end)))
    avg_tab <- bind_cols(res[[1]]$retv %>% select(chrom, start, end), avg_tab)
    
    cov_tab <- map_dfc(res, ~ .x$retv %>% select(ends_with(".cov"), -(chrom:end)))
    cov_tab <- bind_cols(res[[1]]$retv %>% select(chrom, start, end), cov_tab)
    
    return(list(avg = avg_tab, cov = cov_tab))
}

filter_avg_meth <- function(raw_avg_meth, normal_fraction = 0.7, tumor_fraction = 0.7){    
    n_normals <- sum(samp_data$type == "ADJNORMAL")
    f_normals <- rowSums(!is.na(raw_avg_meth[, normal_samples])) >= n_normals * normal_fraction
    n_tumors <- sum(samp_data$type == "TUMOUR")
    f_tumors <- rowSums(!is.na(raw_avg_meth[, tumor_samples])) >= n_tumors * tumor_fraction

    raw_avg_meth <- raw_avg_meth[f_tumors & f_normals, ]
    raw_avg_meth <- raw_avg_meth %>% distinct(chrom, start, end, .keep_all = TRUE)
    message(glue("{scales::comma(nrow(raw_avg_meth))} loci"))
    
    return(raw_avg_meth)
}

loci_meth_summary_stats <- function(avg_meth) {
    normal_mat <- avg_meth %>%
        select(one_of(normal_samples)) %>%
        as.matrix()
    tumor_mat <- avg_meth %>%
        select(one_of(tumor_samples)) %>%
        as.matrix()

    locus_data <- tibble(
        tumor = rowMeans(tumor_mat, na.rm = TRUE),
        n_tumor = rowSums(!is.na(tumor_mat)),
        v_tumor = matrixStats::rowVars(tumor_mat, na.rm = TRUE),
        normal = rowMeans(normal_mat, na.rm = TRUE),
        n_normal = rowSums(!is.na(normal_mat)),
        v_normal = matrixStats::rowVars(normal_mat, na.rm = TRUE)
    )
    locus_data <- bind_cols(avg_meth %>% select(-one_of(c(normal_samples, tumor_samples))), locus_data)
    return(locus_data)
}


filter_avg_matrix <- function(b, min_samples = NULL, var_quant = NULL, min_cov = NULL) {
    f <- rep(TRUE, nrow(b))
    if (!is.null(min_samples)) {
        if (is.null(min_cov)) {
            min_cov <- 1
        }
        message(qq("taking loci that have more than @{min_samples} samples covered by more than @{min_cov}"))
        f <- apply(b[, -c(1:3)], 1, function(x) sum(x >= min_cov, na.rm = T) >= min_samples)
    }

    if (!is.null(var_quant)) {
        sds <- matrixStats::rowSds(as.matrix(b[, -c(1:3)]), na.rm = T)
        f <- f & (sds >= quantile(sds, 1 - var_quant))
    }

    message(qq("@{sum(f)} loci"))
    return(f)
}