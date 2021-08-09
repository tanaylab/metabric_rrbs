get_gene_expression_mat <- function(recalc = FALSE) {
    expr_mat_fn <- here(main_config$data_files$gene_expression_matrix)

    if (!file.exists(expr_mat_fn) || recalc) {
        rna_mat <- fread(main_config$metadata_files$rna_matrix_csv)
        rna_feat <- fread(main_config$metadata_files$rna_feature_csv)
        list2env(main_config$genomic_regions, envir = environment())
        rna_feat_promoters <- as_tibble(define_promoters(main_config$metadata_files$rna_feature_csv, upstream = promoters$upstream, promoters$downstream))

        locus_stats <- get_promoter_summary_stats()
        rna_feat_promoters <- rna_feat_promoters %>% left_join(locus_stats %>% select(chrom, start, end, normal_meth = normal, tumor_meth = tumor, n_tumor_meth = n_tumor, n_normal_meth = n_normal))
        rna_feat_promoters <- rna_feat_promoters %>% mutate(n_samp_expr = rowSums(as.matrix(!is.na(rna_mat[, -1]))))

        chosen_promoters <- rna_feat_promoters %>%
            mutate(
                prev_ord = 1:n(),
                has_expr_data = n_samp_expr > 0,
                has_meth_normal_data = n_normal_meth >= length(normal_samples) * 0.7,
                has_meth_tumor_data = n_tumor_meth >= length(tumor_samples) * 0.7
            ) %>%
            arrange(
                name,
                has_expr_data,
                has_meth_normal_data,
                has_meth_tumor_data,
                normal_meth,
                tumor_meth
            ) %>%
            distinct(name, .keep_all = TRUE)

        rna_mat_f <- bind_cols(chosen_promoters %>% select(chrom, start, end, name, name3.chr), rna_mat[chosen_promoters$prev_ord, ] %>% select(-name.chr))
        colnames(rna_mat_f) <- gsub("MB-", "MB_", colnames(rna_mat_f))
        fwrite(rna_mat_f, expr_mat_fn)
    } else {
        rna_mat_f <- fread(expr_mat_fn)
    }

    return(rna_mat_f)
}

get_gene_tss <- function() {
    # rna_feat <- fread(here(main_config$metadata_files$rna_feature_csv)) %>%
    #     select(chrom = seqnames, start, end, strand, name = name2, name3.chr) %>%
    #     mutate(strand = ifelse(strand == "+", 1, -1)) %>%
    #     select(chrom_expr1 = chrom, start_expr1 = start, end_expr1 = end, strand, full_name = name3.chr)
    # rna_feat <- rna_feat %>%
    #     mutate(tss = ifelse(strand == 1, start_expr1, end_expr1)) %>%
    #     select(full_name, tss)
    return(fread(here("data/gene_tss.tsv")))
}

get_gene_tss_coord <- function() {
    # fread(here(main_config$metadata_files$rna_feature_csv)) %>%
    #     select(chrom = seqnames, start, end, strand, name = name2, full_name = name3.chr) %>%
    #     left_join(get_gene_tss()) %>%
    #     mutate(strand = ifelse(strand == "+", 1, -1)) %>%
    #     mutate(start = tss, end = tss + 1) %>%
    #     select(chrom, start, end, strand, name)
    return(fread(here("data/gene_tss_coords.tsv")))
}

get_mean_expression <- function() {
    {
        expr_mat <- get_gene_expression_mat()
        mean_expr <- expr_mat %>%
            distinct(chrom, start, end, name, name3.chr, .keep_all = TRUE) %>%
            mutate(
                `ER+` = rowMeans(select(., any_of(ER_positive_samples))),
                `ER-` = rowMeans(select(., any_of(ER_negative_samples))),
                normal = rowMeans(select(., any_of(normal_samples)))
            ) %>%
            select(chrom, start, end, name, name3.chr, normal, `ER+`, `ER-`)
     } %cache_df% here("data/mean_expression.tsv") %>% as_tibble()
}

expr_intervs_to_mat <- function(df){
    df %>% select(-(chrom:end), -name3.chr) %>% column_to_rownames("name") %>% as.matrix()
    
}
