get_xist_loci <- function() {        
    {
        xist_expr <- get_gene_expression_mat() %>%
            filter(name == "XIST") %>%
            gather("samp", "expr", -(chrom:name3.chr)) %>% as_tibble()

        all_norm_meth <- fread(here("data/all_norm_meth.tsv")) %>% filter(chrom == "chrX") %>% as_tibble() 
        ER_positive_mat <- all_norm_meth %>% select(chrom:end, any_of(ER_positive_samples)) %>% intervs_to_mat()
        ER_negative_mat <- all_norm_meth %>% select(chrom:end, any_of(ER_negative_samples)) %>% intervs_to_mat()

        xist_expr_vec <- xist_expr %>% select(samp, expr)  %>% deframe()
        samples <- intersect(names(xist_expr_vec), colnames(ER_positive_mat))
        ER_pos_cors <- tgs_cor(t(ER_positive_mat[, samples]), as.matrix(xist_expr_vec[samples]), pairwise.complete.obs = TRUE)

        samples <- intersect(names(xist_expr_vec), colnames(ER_negative_mat))
        ER_neg_cors <- tgs_cor(t(ER_negative_mat[, samples]), as.matrix(xist_expr_vec[samples]), pairwise.complete.obs = TRUE)

        xist_loci <- cbind(ER_pos_cors, ER_neg_cors) %>% mat_to_intervs() %>% filter(V1 >= 0.2 |  V2 >= 0.2) %>% distinct(chrom, start, end)
    } %cache_df% here("data/xist_loci.tsv") %>% as_tibble()
}

get_xist_loci_meth <- function() {
    get_xist_loci() %>%
        inner_join(get_promoter_avg_meth()) %>%
        gather("samp", "meth", -(chrom:end))
}

get_xist_loci_expr <- function() {
    get_gene_expression_mat() %>%
        # mutate(start = start + 1, end = end + 1) %>%
        inner_join(get_xist_loci()) %>%
        gather("samp", "expr", -(chrom:name3.chr))
}

get_xist_cna <- function() {
    xist_cna <- get_xist_loci() %>%
        gintervals.neighbors1(cna %>%
            mutate(end = ifelse(start == end, start + 1, end)) %>%
            filter(chrom == "chrX"), maxneighbors = nrow(samp_data)) %>%
        filter(dist == 0) %>%
        # select(chrom, start, end, samp, cna = cna_val)
        select(chrom, start, end, samp, cna = cna_round)

    xist_cna <- xist_cna %>%
        mutate(cna_grp = cut(cna, breaks=c(-1, 0,1,2,10), labels=c("0N", "1N", "2N", ">=3N"))) %>% 
        # mutate(cna_grp = cut(cna, breaks = c(0, 0.25, 0.75, 1.25, 10), labels = c("0N", "1N", "2N", ">=3N"), include.lowest = T)) %>%
        # mutate(cna_grp = cut(cna, breaks = c(0, 0.25, 0.75, 1.25, 10), labels = c("0N", "1N", "2N", ">=3N"), include.lowest = T)) %>%
        filter(cna_grp != "0N") %>%
        left_join(samp_data %>% select(samp, ER = ER1))

    return(xist_cna)
}

get_xist_expr_cna <- function() {
    get_xist_cna() %>%
        inner_join(get_xist_loci_expr()) %>%
        group_by(chrom, start, end, ER, cna_grp) %>%
        summarise(expr = mean(expr, na.rm = TRUE)) %>%
        ungroup() %>%
        spread(cna_grp, expr)
}

get_xist_meth_cna <- function() {
    get_xist_cna() %>%
        inner_join(get_xist_loci_meth()) %>%
        group_by(chrom, start, end, ER, cna_grp) %>%
        summarise(meth = mean(meth, na.rm = TRUE)) %>%
        ungroup() %>%
        spread(cna_grp, meth)
}