annotate_cis_cands <- function(cands, sigma_meth = 2, sigma_expr = 2, meth_diff =0.2, expr_diff = 1){
    norm_meth <- fread(here("data/all_norm_meth.tsv")) %>% as_tibble()  %>% intervs_to_mat()  
    raw_meth <- get_all_meth() %>% intervs_to_mat()
    expr_mat_all <- get_gene_expression_mat() %>% expr_intervs_to_mat()    

    res <- map_dfr(unique(cands$ER), function(ER) {
        cands <- cands %>% filter(ER == !!ER)
        ER_samples <- samp_data$samp[samp_data$ER1 == ER]

        samples <- intersect(ER_samples, colnames(norm_meth))
        samples <- samples[!is.na(samples)]     
        cands <- cands  %>% unite("coord", chrom:end, remove = FALSE) 
        coord <- cands$coord        

        meth_mat <- norm_meth[unique(coord), samples]        

        samples <- intersect(ER_samples, colnames(expr_mat_all)) 
        samples <- samples[!is.na(samples)]     
        expr_mat <- as.matrix(expr_mat_all[unique(cands$gene), samples])
        
        m <- rowMeans(meth_mat, na.rm=TRUE)
        sd <- sigma_meth*matrixStats::rowSds(meth_mat, na.rm=TRUE)
        names(sd) <- names(m)
        m_meth <- enframe(m, "coord", "mean_meth")
        sd_meth <- enframe(sd, "coord", "sd_meth")
        n_hypo <- enframe(rowSums(meth_mat <= (m - sd), na.rm=TRUE), "coord", "n_hypometh")
        n_hyper <- enframe(rowSums(meth_mat >= (m + sd), na.rm=TRUE), "coord", "n_hypermeth")

        m <- rowMeans(expr_mat, na.rm=TRUE)
        sd <- sigma_meth*matrixStats::rowSds(expr_mat, na.rm=TRUE)
        names(sd) <- names(m)
        m_expr <- enframe(m, "gene", "mean_expr")
        sd_expr <- enframe(sd, "gene", "sd_expr")
        n_induced <- enframe(rowSums(expr_mat >= (m + sd), na.rm=TRUE), "gene", "n_induced")
        n_repressed <- enframe(rowSums(expr_mat <= (m - sd), na.rm=TRUE), "gene", "n_repressed")
        
        N_mat <- (1*!is.na(meth_mat)) %*% (1*t(!is.na(expr_mat)))
        
        
        N <- map_dbl(1:nrow(cands), ~ N_mat[coord[.x], cands$gene[.x]])

        cands$N_considered <- N
        
        cands <- cands %>% left_join(n_hypo) %>% left_join(n_induced) %>% left_join(n_hyper) %>% left_join(n_repressed) %>% left_join(m_meth) %>% left_join(sd_meth) %>% left_join(m_expr) %>% left_join(sd_expr) 


        # Annotate vs normal
        raw_meth <- raw_meth[unique(coord), ]
        normal_meth <- raw_meth[, normal_samples]
        normal_mean <- rowMeans(normal_meth, na.rm=TRUE)
        normal_sd <- matrixStats::rowSds(normal_meth, na.rm=TRUE)

        samples <- intersect(ER_samples, colnames(raw_meth))
        raw_meth_mat <- raw_meth[, samples]
        n_hypermeth_vs_normal <- rowSums((raw_meth_mat >= normal_mean + meth_diff) & (raw_meth_mat >= normal_mean +sigma_meth*normal_sd), na.rm=TRUE)
        n_hypometh_vs_normal<- rowSums((raw_meth_mat <= normal_mean - meth_diff) & (raw_meth_mat <= normal_mean -sigma_meth*normal_sd), na.rm=TRUE)
        n_stable_vs_normal <- rowSums(!is.na(raw_meth_mat)) - n_hypermeth_vs_normal - n_hypometh_vs_normal

        normal_meth_df <- tibble(coord = rownames(raw_meth_mat), normal_meth = normal_mean, normal_meth_sd = normal_sd, n_hypermeth_vs_normal = n_hypermeth_vs_normal, n_hypometh_vs_normal = n_hypometh_vs_normal, n_stable_vs_normal = n_stable_vs_normal)
        
        cands <- cands %>% left_join(normal_meth_df, by = "coord")

        samples <- intersect(normal_samples, colnames(expr_mat_all)) 
        normal_expr <- expr_mat_all[unique(cands$gene), samples]
        normal_expr_mean <- rowMeans(normal_expr, na.rm=TRUE)
        normal_expr_sd <- matrixStats::rowSds(normal_expr, na.rm=TRUE)
        

        n_induced_vs_normal <- rowSums((expr_mat >= normal_expr_mean + expr_diff) & (expr_mat >= normal_expr_mean +sigma_expr*normal_expr_sd), na.rm=TRUE)
        n_repressed_vs_normal <- rowSums((expr_mat <= normal_expr_mean - expr_diff) & (expr_mat >= normal_expr_mean -sigma_expr*normal_expr_sd), na.rm=TRUE)
        n_stable_expr_vs_normal <- rowSums(!is.na(expr_mat)) - n_induced_vs_normal - n_repressed_vs_normal

        normal_expr_df <- tibble(gene = rownames(expr_mat), normal_expr = normal_expr_mean, normal_expr_sd = normal_expr_sd, n_induced_vs_normal = n_induced_vs_normal, n_repressed_vs_normal = n_repressed_vs_normal, n_stable_expr_vs_normal = n_stable_expr_vs_normal)

        cands <- cands %>% left_join(normal_expr_df, by = "gene")

        cands <- cands %>% select(-coord)

        

        return(cands)
    })
}

