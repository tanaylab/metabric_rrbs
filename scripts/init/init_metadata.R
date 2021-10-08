define_promoters <- function(refseq_file, upstream = 500, downstream = 50) {
    refseq <- fread(refseq_file) %>%
        select(chrom = seqnames, start, end, strand, name = name2, name3.chr) %>%
        mutate(start = ifelse(strand == "+", start, end)) %>%
        mutate(end = start + 1, strand = case_when(strand == "+" ~ 1, strand == "-" ~ -1)) %>%
        select(chrom, start, end, strand, name, name3.chr) %>%
        filter(chrom %in% gintervals.all()$chrom)
    # account for misha being 0 based
    refseq <- refseq %>% mutate(start = start - 1, end = end - 1)
    promoters <- refseq %>%
        as.data.frame() %>%
        mutate(start = ifelse(strand == 1, start - upstream, start - downstream), end = ifelse(strand == 1, end + downstream, end + upstream)) %>%
        gintervals.force_range()
    return(promoters)
}


define_genomic_regions <- function(config_file) {
    conf <- yaml::read_yaml(config_file)
    list2env(conf$genomic_regions, envir = environment())

    k4me3_peaks <<- gpatterns.putative_enhancers(promoters$H3K4me3_tracks, quant_thresh = promoters$H3K4me3_thresh, normalize = NULL, min_tss_dist = 0) %cache_df% here(promoters$k4me3_file)
   
    
    promoter_intervs <<- define_promoters(refseq_file = promoters$refseq_file, upstream = promoters$upstream, promoters$downstream) %>% 
        gintervals.neighbors1(k4me3_peaks) %>%
        mutate(k4me3 = dist == 0) %>%
        select(-(chrom1:dist))  %cache_df% here(promoters$file)

    promoter_intervs_extended <<- promoter_intervs %>% mutate(start = ifelse(strand == 1, start - 500, start), end = ifelse(strand == 1, end, end + 500))

    k27ac_peaks <<- gpatterns.putative_enhancers(enhancers$H3K27ac_tracks, quant_thresh = enhancers$H3K27ac_thresh, normalize = enhancers$H3K27ac_size, min_tss_dist = enhancers$H3K27ac_tss_dist) %cache_df% here(enhancers$k27ac_file)
    

    enh_intervs <<- gpatterns.putative_enhancers(enhancers$H3K4me1_tracks, quant_thresh = enhancers$H3K4me1_thresh, normalize = enhancers$H3K4me1_size, min_tss_dist = enhancers$H3K4me1_tss_dist) %cache_df% here(enhancers$file)

    enh_intervs_tumors <<- gpatterns.putative_enhancers(enhancers$H3K4me1_tumor_tracks, quant_thresh = enhancers$H3K4me1_thresh, normalize = enhancers$H3K4me1_size, min_tss_dist = enhancers$H3K4me1_tss_dist) %cache_df% here(enhancers$file_tumors)
    

    k27_peaks <<- gpatterns.putative_enhancers(polycomb$H3K27me3_track, quant_thresh = polycomb$H3K27me3_thresh, normalize = NULL) %cache_df% here(polycomb$file)  
    
        
    exon_intervs <<- gintervals.load("intervs.global.exon")
        
}


init_metadata <- function(config_file, recalc = FALSE) {
    conf <- yaml::read_yaml(config_file)

    samp_data <<- fread(here("data/samp_data.csv")) %>% as_tibble()    

    tumor_samples <<- samp_data$samp[samp_data$type == "TUMOUR"]
    normal_samples <<- samp_data$samp[samp_data$type == "ADJNORMAL"]
    ER_positive_samples <<- samp_data$samp[samp_data$ER == "positive" & !is.na(samp_data$ER)]
    ER_negative_samples <<- samp_data$samp[samp_data$ER == "negative" & !is.na(samp_data$ER)]

    init_colors(conf$colors)
    
    survival <<- samp_data %>%         
        mutate(y = T / 365) %>%
        select(samp, ER, ER1, age, death, T, y) %>% 
        as_tibble()

    cna <<- fread(here("data/cna.tsv")) %>% as_tibble() 
}

init_colors <- function(color_config) {
    list2env(color_config, envir = environment())
    color_config <<- color_config
    annot_colors <<- map(color_config, unlist)
}
