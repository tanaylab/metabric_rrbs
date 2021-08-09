get_cis_genomic_examples_cg_meth <- function(genes = c("DNMT3A", "GATA3", "TBX1", "FGFR4", "PAX8"), scope = 2e4, min_dist = 1e5, min_tss_dist = 2e3) {    
    {
        cands <- bind_rows(
            fread(here("data/genomic_cis_cands_ER_positive.tsv")),
            fread(here("data/genomic_cis_cands_ER_negative.tsv")),
            fread(here("data/genomic_cis_cands_normal.tsv"))) %>% as_tibble()

        top_cands <- cands %>% filter(gene %in% genes, type == "obs", rank == 1,  !is.na(dist), abs(dist) <= min_dist, abs(dist) >= min_tss_dist) %>% arrange(gene, cor) %>% group_by(gene) %>% slice(1)

        rm(cands)    
        gc()

        intervs <- top_cands %>% mutate(start = pmin(start, start_expr), end = pmax(end, start_expr)) %>% select(chrom, start, end, everything()) %>% gpatterns:::gintervals.centers() %>% gpatterns:::gintervals.expand(scope) %>% select(chrom, start, end)

        meth_df <- gextract_meth_parallel(samp_data$track, samp_data$samp, intervals = intervs, iterator = "intervs.global.seq_CG", min_cov = 1)        
          
        cg_meth <- meth_df$cov %>% gather("samp", "cov", -(chrom:end)) %>% mutate(samp = gsub("\\.cov", "", samp))
        cg_meth1 <- meth_df$avg %>% gather("samp", "avg", -(chrom:end))
        stopifnot(all(cg_meth$chrom == cg_meth1$chrom))
        stopifnot(all(cg_meth$start == cg_meth1$start))
        stopifnot(all(cg_meth$end == cg_meth1$end))
        stopifnot(all(cg_meth$samp == cg_meth1$samp))
        
        cg_meth <- bind_cols(cg_meth, cg_meth1 %>% select(avg))
        cg_meth <- cg_meth %>% filter(!is.na(cov), !is.na(avg)) %>% mutate(meth = avg * cov)
        cg_meth <- cg_meth %>% left_join(samp_data %>% select(samp, ER = ER1))    

        cg_meth

    }  %cache_df% here("data/cis_genomic_examples_cg_meth.tsv") %>% as_tibble()        
}

plot_cis_genomic_example <- function(df, gene, expr_df, meth_df, ofn=NULL, k_smooth = 40, ER = 'ER+', scope_start = 1e3, scope_end = 1e4, min_cov = 10, add_pval = TRUE){

    df <- df %>% filter(gene == !!gene, ER == !!ER) %>% arrange(cor) %>% slice(1)
    samples <- samp_data %>% filter(ER1 == !!ER) %>% pull(samp)
    scope <- data.frame(chrom = df$chrom[1], start = min(df$start - scope_start, df$start_expr - scope_start), end = max(df$end + scope_end, df$start_expr + scope_end))
    expr_df <- expr_df %>% filter(gene == !! gene)
    samp_ord <-  expr_df %>% filter(gene == !!gene)  %>% arrange(expr) %>% pull(samp)    

    meth_df <- meth_df %>% 
        filter(samp %in% samples) %>% #, cov >= min_cov) %>%
        mutate(meth = meth / cov) %>% 
        # select(-cov) %>%
        mutate(samp = factor(samp, levels=rev(samp_ord))) %>% 
        filter(!is.na(samp)) %>%         
        select(chrom, start, end, everything()) %>% 
        gintervals.neighbors1(scope) %>%
        filter(dist == 0) %>%
        select(-(chrom1:dist)) 
        
    meth_df_f <- meth_df %>%            
        gintervals.neighbors1(df %>% select(chrom, start, end)) %>%
        filter(dist == 0) %>% 
        select(-(chrom1:dist))
    
    # Select the CpG with maximum correlation
    meth_df_f <- meth_df_f %>% 
        left_join(expr_df) %>% 
        group_by(chrom, start, end) %>% 
        mutate(cor = cor(meth, expr, use = 'pairwise.complete.obs', method="spearman")) %>% 
        ungroup() %>% 
        filter(abs(cor) == max(abs(cor)))
        
    locus_cor <- meth_df_f$cor[1]

    meth_df_f <- meth_df_f %>%     
        filter(cov >= min_cov) %>%
        mutate(samp = factor(samp, levels = rev(samp_ord))) %>%
        arrange(samp) %>%
        mutate(meth_smooth = c(rep(mean(meth[1:(k_smooth - 1)], na.rm = TRUE), k_smooth - 1), zoo::rollmean(meth, k = k_smooth, na.rm = TRUE)))

    enh_loc <- meth_df_f$start     


    p_line <- meth_df_f %>%
        ggplot(aes(x = expr, y = meth_smooth, group = 1)) +
        # ggrastr::geom_point_rast(aes(y = meth, color=ER), size = 0.1) +
        geom_point(aes(y = meth, color=ER), size = 0.1) +
        scale_color_manual(values=annot_colors$ER1) + 
        guides(color=FALSE) + 
        geom_line() +        
        ylab('Avg. methylation') + 
        xlab('Gene expression') + 
        ylim(0, 1) + 
        ggtitle(glue("{scope$chrom} ({scope$start}-{scope$end})")) + 
        theme(plot.title = element_text(size = 5), aspect.ratio=1)


    if (add_pval) {        
        p_line <- p_line + annotate('text', label = sprintf("~ rho == %0.2f", locus_cor), parse = TRUE, x = 6.5, y = 0.2, size = 3, family = "Arial")        
    }

    p_genome <- meth_df %>%
        distinct(start) %>%
        ggplot(aes(x = start, xend = start, y = 1, yend = 0)) +
            geom_segment(size=0.1, alpha=0.5) +
            geom_point(size=0.1, alpha=0.5) +            
            ylim(-0.2, 5) +
            scale_x_continuous(limits = c(scope$start, scope$end)) +
            xlab('') +
            ylab('') +
            geom_hline(yintercept = 0) +
            ggpubr::theme_pubr() +
            theme(text = element_text(family = 'Arial', size = 6),
                  axis.ticks.y = element_blank(),
                  axis.text.y = element_blank(),
                  axis.line.y = element_blank(),
                  axis.line.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank()) +
            annotate('segment',
                     x = enh_loc,
                     xend = enh_loc,
                     y = 3.5,
                     yend = 2,
                     color = 'red',
                     arrow = arrow(length = unit(0.05, 'inches'))) +
            annotate('segment',
                     x = df$start_expr[1],
                     xend = df$start_expr[1],
                     y = 1.5,
                     yend = 2.7,
                     color = 'blue') +
            annotate('segment',
                    x = df$start_expr[1],
                    xend = df$start_expr[1] + df$strand_expr * (scope$end - scope$start) * 0.1,
                    y = 2.7,
                    yend = 2.7,
                    color = 'blue',
                    arrow = arrow(length = unit(0.05, 'inches'))) +
            annotate('text', label = df$gene[1], x = df$start_expr[1], y = 3.5, size = 2)

    

    margin <- margin(t = 0, b = 0, r = 3, l = 3, unit = "pt")
    p <- plot_grid(
        p_genome + theme(plot.margin=margin), 
        p_line + theme(plot.margin=margin), 
        ncol=1, 
        rel_heights=c(0.3,0.7), axis="l")

    if (!is.null(ofn)) {
        p <- p + ggsave(ofn, width = 7, height = 6)
    }
    return(p)
    

}
