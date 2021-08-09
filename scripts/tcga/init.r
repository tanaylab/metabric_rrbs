get_tcga_brca_expr_mat <- function(){
    mat <- readr::read_rds(here("data/TCGA-BRCA/TCGA_BRCA_expr.rds"))
    mat <- log10(mat + 1e-9)
    return(mat)
}


get_tcga_brca_meth <- function(intervals = gintervals.all(), track = "TCGA.BRCA_450k"){
    opt <- options(gmax.data.size=1e9, gmultitasking=FALSE)
    on.exit(options(opt))
    gvtrack.array_extract(track, intervals, iterator=intervals) %>% arrange(intervalID) %>% select(-intervalID) %>% as_tibble() 
}

get_all_tcga_brca_meth <- function(){
    gtrack.array.extract("TCGA.BRCA_450k", gtrack.array.get_colnames("TCGA.BRCA_450k"), gintervals.all()) %cache_df% here("data/TCGA_BRCA_450k_all_meth.tsv") %>% as_tibble() 
}

get_tcga_brca_prom_meth <- function(){
    get_tcga_brca_meth(promoter_intervs %>% distinct(chrom, start, end)) %cache_df% here("data/TCGA-BRCA/TCGA_BRCA_450k_promoters.tsv") %>% intervs_to_mat() 
}

get_tcga_brca_prom_meth_TME_norm <- function(){
    get_tcga_brca_meth(promoter_intervs %>% distinct(chrom, start, end), track = "TCGA.BRCA_450k_norm") %cache_df% here("data/TCGA_BRCA_450k_promoters_TME_norm.tsv") %>% intervs_to_mat() 
}

get_tcga_brca_genomic_meth <- function(track = "TCGA.BRCA_450k"){
    genomic_intervs <- gintervals.diff(gintervals.all(), "intervs.global.exon") %>% gintervals.diff(get_promoters(upstream = 2000, downstream = 50))
    tcga_meth_genomic <- gtrack.array.extract(track, intervals = genomic_intervs) %>% arrange(intervalID) %>% select(-intervalID) %>% intervs_to_mat()
    return(tcga_meth_genomic)
}

# get_tcga_brca_genomic_meth_immune_caf_norm <- function(){
#     return(cbind(
#         readr::read_rds(here("data/tcga_brca_genomic_er_positive_norm.rds"),
#         readr::read_rds(here("data/tcga_brca_genomic_er_negative_norm.rds"),
#         readr::read_rds(here("data/tcga_brca_genomic_er_normal_norm.rds")
#         ))    
# }


init_tcga_samp_data <- function(){
    tcga_samp_data <<- fread(here("data/TCGA-BRCA/tcga_samp_data.csv")) %>% as_tibble()
    tcga_ER_positive_samples <<- tcga_samp_data %>% filter(gender == "FEMALE", ER == "ER+") %>% pull(samp_id)
    tcga_ER_negative_samples <<- tcga_samp_data %>% filter(gender == "FEMALE", ER == "ER-") %>% pull(samp_id)
    tcga_normal_samples <<- tcga_samp_data %>% filter(gender == "FEMALE", ER == "normal") %>% pull(samp_id)
}


load_all_tcga_brca_data <- function(){
    # Load expression data
    tcga_expr <<- get_tcga_brca_expr_mat()
    tcga_expr_positive <<- tcga_expr[, intersect(colnames(tcga_expr), tcga_ER_positive_samples)]
    tcga_expr_negative <<- tcga_expr[, intersect(colnames(tcga_expr), tcga_ER_negative_samples)]
    tcga_expr_normal <<- tcga_expr[, intersect(colnames(tcga_expr), tcga_normal_samples)]
    
    # Load promoter methylation data
    tcga_prom_meth <<- get_tcga_brca_prom_meth()
    tcga_prom_meth_positive <<- tcga_prom_meth[, intersect(colnames(tcga_prom_meth), tcga_ER_positive_samples)]
    tcga_prom_meth_negative <<- tcga_prom_meth[, intersect(colnames(tcga_prom_meth), tcga_ER_negative_samples)]
    tcga_prom_meth_normal <<- tcga_prom_meth[, intersect(colnames(tcga_prom_meth), tcga_normal_samples)]
    
    # Load genomic methylation data
    tcga_genomic_meth <<- get_tcga_brca_genomic_meth()
    tcga_genomic_meth_positive <<- tcga_genomic_meth[, intersect(colnames(tcga_genomic_meth), tcga_ER_positive_samples)]
    tcga_genomic_meth_negative <<- tcga_genomic_meth[, intersect(colnames(tcga_genomic_meth), tcga_ER_negative_samples)]
    tcga_genomic_meth_normal <<- tcga_genomic_meth[, intersect(colnames(tcga_genomic_meth), tcga_normal_samples)]
    
    # Combine promoters and genomic methylation
    all_meth <<- get_all_tcga_brca_meth() %>% intervs_to_mat()
    all_meth_positive <<- all_meth[, intersect(colnames(all_meth), tcga_ER_positive_samples)]
    all_meth_negative <<- all_meth[, intersect(colnames(all_meth), tcga_ER_negative_samples)]
    all_meth_normal <<- all_meth[, intersect(colnames(all_meth), tcga_normal_samples)]
}
