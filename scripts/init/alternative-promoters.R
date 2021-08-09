#' Resolve alternative promoters by taking the promoter with lowest methylation in normal samples
resolve_alt_promoters <- function(intervals){   
    intervals %>% 
        left_join(get_all_summary_meth(), by = c("chrom", "start", "end")) %>% 
        left_join(promoter_intervs %>% select(chrom, start, end, name), by = c("chrom", "start", "end")) %>% 
        filter(!is.na(name)) %>% 
        arrange(name, normal) %>% 
        group_by(name) %>% 
        slice(1) %>% 
        select(chrom, start, end, name) %>% 
        ungroup()
}