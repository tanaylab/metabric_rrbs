get_promoter_cand_interval <- function(cands, name, ER) {
    cands <- cands %>% left_join(promoter_intervs %>% select(chrom, start, end, name, strand, full_name = name3.chr)) %>% left_join(get_gene_tss())
    cands %>%
        filter(name == !!name, ER == !!ER) %>%
        mutate(min_cor = min(cor, na.rm = TRUE)) %>%
        filter(cor == min_cor) %>%        
        select(chrom, tss, name, strand) %>%
        mutate(start = tss, end = start + 1) %>%
        select(chrom, start, end, strand, name)
}


get_cis_promoter_examples_cg_meth <- function(genes = c("KRT7", "CABP4", "BRCA1"), scope = 1e4) {        
    {
        cands <- fread(here("data/promoter_cis_cands.tsv")) %>% as_tibble()
        example_cands <- map_dfr(genes, ~ get_promoter_cand_interval(cands, .x, "ER+"))

        intervs <- example_cands %>%
            mutate(start = start - scope, end = end + scope) %>%
            select(chrom, start, end)
        meth_df <- gextract_meth_parallel(samp_data$track, samp_data$samp, intervals = intervs, iterator = "intervs.global.seq_CG", min_cov = 1)

        cg_meth <- meth_df$cov %>%
            gather("samp", "cov", -chrom, -start, -end) %>%
            mutate(samp = gsub("\\.cov", "", samp)) %>% 
            left_join(meth_df$avg %>% 
                gather("samp", "avg", -chrom, -start, -end) 
            ) %>% 
            filter(!is.na(cov), !is.na(avg)) %>%
            mutate(meth = cov * avg) %>% 
            left_join(samp_data %>% select(samp, ER = ER1))
        cg_meth
    } %cache_df% here("data/cis_promoter_examples_cg_meth.tsv") %>% as_tibble()
}

plot_cis_promoter_example <- function(df, cg_meth, ER, resolution = 1e3, min_n = 50, k = 5, min_cov = 10, plot_all_tss_lines = FALSE) {
    scope <- df %>% mutate(start = start - resolution, end = end + resolution)

    x <- cg_meth %>%
        gintervals.neighbors1(scope) %>%
        filter(!is.na(chrom1)) %>%
        filter(dist == 0) %>%
        select(-(chrom1:end1))

    expr_df <- get_gene_expression_mat() %>%
        filter(name == df$name) %>%
        gather("samp", "expr", -(chrom:name3.chr)) %>%
        select(name, samp, expr)

    x <- x %>%
        left_join(expr_df) %>%
        filter(!is.na(expr))

    gene_cors <- x %>%
        filter(ER == !!ER, cov >= min_cov) %>%
        mutate(avg_m = meth / cov) %>%
        group_by(chrom, start, end) %>%
        summarise(cor = cor(expr, avg_m, use = "pairwise.complete.obs"), avg_m = sum(meth, na.rm = TRUE) / sum(cov, na.rm = TRUE), n = n()) %>%
        ungroup() %>%
        right_join(x %>% distinct(chrom, start, end))

    gene_cors1 <- gene_cors %>%
        filter(n >= min_n, !is.na(cor)) %>%
        mutate(cor_m = c(cor[1:(k - 1)], zoo::rollmean(cor, k = k)))
    
    p_scatter <- gene_cors1 %>%
        ggplot(aes(x = start, y = cor, color = avg_m)) +
        geom_hline(color = "darkgray", yintercept = 0) +
        geom_vline(color = "darkblue", xintercept = df$start, linetype = "dashed") +
        geom_line(inherit.aes = FALSE, aes(x = start, y = cor_m)) +
        geom_point(size = 1) +
        xlim(scope$start, scope$end) +
        scale_color_gradientn("Average
Methylation", colors = rev(c("darkred", "gray", "cyan")), na.value = "transparent", limits = c(0, 1)) +
        xlab("") +
        ylab("Expression-Methylation
correlation") +
        ylim(-1, 1) +
        theme(legend.key.width = unit(0.3, "cm"))

    cpg_intervs <- giterator.intervals(intervals=scope, iterator="intervs.global.seq_CG")

    p_genome <- cpg_intervs %>%
    # p_genome <- gene_cors1 %>%
        distinct(start) %>%
        ggplot(aes(x = start, xend = start, y = 1, yend = 0)) +
        geom_segment() +
        geom_point(size = 0.5) +
        ylim(0, 5) +
        scale_x_continuous(limits = c(scope$start, scope$end)) +
        xlab("") +
        ylab("") +
        ggpubr::theme_pubr() +
        theme(
            text = element_text(family = "Arial", size = 6),
            axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            axis.line.y = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank()
        )
    
    genes <- get_gene_tss_coord() %>%
        gintervals.neighbors1(scope) %>%
        filter(dist == 0) %>%
        filter(name != df$name) %>%
        distinct(name, .keep_all = TRUE) %>%
        select(chrom, start, end, strand, name) %>%
        bind_rows(df)



    for (i in 1:nrow(genes)) {
        p_genome <- p_genome +
            annotate("segment",
                x = genes$start[i],
                xend = genes$start[i],
                y = 1.5,
                yend = 2.7,
                color = "blue"
            ) +
            annotate("segment",
                x = genes$start[i],
                xend = genes$start[i] + genes$strand[i] * (scope$end - scope$start) * 0.05,
                y = 2.7,
                yend = 2.7,
                color = "blue",
                arrow = arrow(length = unit(0.05, "inches"))
            ) +
            annotate("text", label = genes$name[i], x = genes$start[i], y = 3.9, size = 2)
        if (plot_all_tss_lines && genes$name[i] != df$name) {
            p_scatter <- p_scatter +
                geom_vline(color = "darkblue", xintercept = genes$start[i], linetype = "dashed")
        }
    }

    margin <- margin(t = 0, b = 0, r = 3, l = 3, unit = "pt")
    p <- plot_grid(
        plot_grid(
            p_genome + theme(plot.margin = margin), 
            p_scatter + guides(color = FALSE) + theme(plot.margin = margin), ncol = 1, align = "v", rel_heights = c(0.3, 0.7)
        ), 
        get_legend(p_scatter), 
        nrow = 1, 
        rel_widths = c(0.8, 0.2))

    return(p)
}
