annotate_loci <- function(loci) {
    k4me1_names <- c("k4me1_luminal", "k4me1_myo1", "k4me1_myo2", "k4me1_hmec")
    k4me1_tracks <- main_config$genomic_regions$enhancers$H3K4me1_tracks
    walk2(k4me1_names, k4me1_tracks, ~ gvtrack.create(.x, .y, func = "global.percentile.max"))

    gvtrack.create("k27me3", main_config$genomic_regions$polycomb$H3K27me3_track, func = "global.percentile.max")
    gvtrack.create("tss_d", "intervs.global.tss", func = "distance")
    loci_annot <- gextract.left_join(c("seq.CG_500_mean", main_config$genomic_regions$tor_track, "tss_d", "k27me3", k4me1_names), iterator = loci, intervals = loci, colnames = c("cg_cont", "tor", "tss_d", "k27me3", k4me1_names)) %>%
        select(-(chrom1:end1)) %>%
        as_tibble()

    return(loci_annot)
}

annotate_loci_tumors_datasets <- function(loci){
    opt <- options(gmax.data.size = 1e9, gmultitasking = FALSE)
    on.exit(options(opt))
    k4me1_tracks <- main_config$genomic_regions$enhancers$H3K4me1_tumor_tracks
    k4me1_names <- gsub("GSE85158.", "", k4me1_tracks)
    walk2(k4me1_names, k4me1_tracks, ~ gvtrack.create(.x, .y, func = "global.percentile.max"))
    
    k27me3_tracks <- main_config$genomic_regions$polycomb$H3K27me3_track_tumors
    k27me3_names <- gsub("GSE85158.", "", k27me3_tracks)
    walk2(k27me3_names, k27me3_tracks, ~ gvtrack.create(.x, .y, func = "global.percentile.max"))
    
    loci_annot <- gextract.left_join(c(k27me3_names, k4me1_names), iterator = loci, intervals = loci, colnames = c(k27me3_names, k4me1_names)) %>%
        select(-(chrom1:end1)) %>%
        as_tibble()
}

get_loci_annot <- function(){
	get_all_meth() %>% 
		select(chrom, start, end) %>% 
		annotate_loci() %cache_df% 
		here("data/loci_annotation.tsv")
}

get_loci_annot_tumor_datasets <- function(){
    get_all_meth() %>% 
        select(chrom, start, end) %>% 
        annotate_loci_tumors_datasets() %cache_df% 
        here("data/loci_annotation_tumor_datasets.tsv")
}