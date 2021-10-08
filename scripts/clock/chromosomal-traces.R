
get_tor_clock_chrom_trace <- function(chrom, ER, iterator = 1e4) {
    {
        intervals <- gintervals.all() %>% filter(chrom == !! chrom)
        feats_df <- get_all_features()
        chr_bins <- gextract(main_config$genomic_regions$tor_track, intervals, iterator = iterator, colnames = "tor") %>%
            select(-intervalID) %>%
            as_tibble()
            
        loss_meth_grp <- feats_df %>%
            filter(ER == !!ER) %>%
            mutate(loss_grp = cut(clock, quantile(clock, c(0, 0.3, 0.7, 1)), include.lowest = TRUE, labels = c("low", "mid", "high")))
        low_tumors <- loss_meth_grp %>%
            filter(ER == !!ER, loss_grp == "low") %>%
            pull(samp)
        low_tumors_tracks <- samp_data %>% filter(samp %in% low_tumors) %>% pull(track)
        
        high_tumors <- loss_meth_grp %>%
            filter(ER == !!ER, loss_grp == "high") %>%
            pull(samp)
        high_tumors_tracks <- samp_data %>% filter(samp %in% high_tumors) %>% pull(track)
        
        normal_tracks <- samp_data %>% filter(samp %in% normal_samples) %>% pull(track)        
        
        meth_df <- gextract_meth_parallel(c(low_tumors_tracks, high_tumors_tracks, normal_tracks), c(low_tumors, high_tumors, normal_samples), intervals = intervals, iterator = iterator, min_cov = 30)      
        
        meth_df_sum <- meth_df$avg %>% mutate(
                low_loss = rowMeans(select(., any_of(low_tumors)), na.rm=TRUE),
                high_loss = rowMeans(select(., any_of(high_tumors)), na.rm=TRUE),
                normal = rowMeans(select(., any_of(normal_samples)), na.rm=TRUE),
                low_loss_n = rowSums(!is.na(select(., any_of(low_tumors)))),
                high_loss_n = rowSums(!is.na(select(., any_of(high_tumors)))),
                normal_n = rowSums(!is.na(select(., any_of(normal_samples))))
            ) %>% 
            select(chrom, start, end, low_loss, high_loss, normal,  low_loss_n, high_loss_n, normal_n)
            
        df <- meth_df_sum %>% left_join(chr_bins)
        
        df              
    } %cache_df% here(glue("data/tor_clock_chrom_trace_{chrom}_{iterator}.tsv", iterator = format(iterator, scientific = FALSE))) %>% as_tibble()
}

plot_tor_clock_chrom_track <- function(chrom, ER, iterator = 1e4, min_n = 50, trace_df = NULL, tor_breaks = main_config$genomic_regions$tor_high_low){
    if (is.null(trace_df)){
        trace_df <- get_tor_clock_chrom_trace(chrom, ER, iterator)            
    }
    
    trace_df <- trace_df %>% 
        mutate(
            low_loss = ifelse(low_loss_n >= min_n, low_loss, NA), 
            high_loss = ifelse(high_loss_n >= min_n, high_loss, NA), 
            normal = ifelse(normal_n >= min_n, normal, NA)
        ) %>% 
        filter(!is.na(tor))
        
    trace_df <- trace_df %>% 
        select(chrom, start, end, tor, low_loss, high_loss, normal) %>% 
        gather("loss_grp", "avg", -(chrom:end), -tor)
            

    chrom_len <- gintervals.all() %>% filter(chrom == !!chrom) %>% pull(end)
    
    margin <- margin(t = 0, b = 0, r = 1, l = 1, unit = "pt")
    p_chr_meth <- trace_df %>%
        group_by(loss_grp) %>%
        mutate(avg = zoo::rollapply(avg, width = 50, FUN = function(x) mean(x, na.rm = TRUE), fill = NA, align = "right")) %>%
        ungroup() %>%
        mutate(loss_grp = factor(loss_grp, levels = c("normal", "low_loss", "high_loss"))) %>%
        arrange(loss_grp) %>%
        ggplot(aes(x = start, y = avg, color = loss_grp, alpha = loss_grp)) + 
            geom_line(size=0.2) + 
            scale_color_manual(values = c("high_loss" = "red", "low_loss" = "darkblue", normal = "black")) +
            scale_alpha_manual(values = c("high_loss" = 1, "low_loss" = 1, "normal" = 0.5)) + 
            xlab("") + 
            scale_x_continuous(expand = c(0, 0), limits=c(0, chrom_len)) +
            ylab("Methylation") + 
            xlab(gsub("chr", "Chromosome ", chrom)) + 
            theme(
                panel.grid.major.x = element_blank(),
                panel.grid.minor.x = element_blank()
            )

    p_chr_tor <- trace_df %>%
        mutate(tor_grp = cut(tor, breaks = tor_breaks, labels = c("late", "early"))) %>%
        ggplot(aes(x = start, y = 1, fill = tor_grp)) + 
            theme_void() + 
            geom_raster() + 
            scale_fill_manual(values = c("late" = "black", "early" = "darkgreen")) + 
            scale_x_continuous(expand = c(0, 0), limits=c(0, chrom_len))
            
    p_chr <- plot_grid(
        p_chr_tor + 
            theme(plot.margin = margin) + 
            guides(fill = FALSE), 
        p_chr_meth + 
            theme(plot.margin = margin) + 
            guides(color = FALSE, alpha = FALSE), ncol = 1, align = "hv", rel_heights = c(0.4, 0.6))
    return(list(p_meth = p_chr_meth, p_tor = p_chr_tor, p = p_chr))
}