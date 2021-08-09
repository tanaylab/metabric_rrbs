get_gene_cna_df <- function() {
    prom.gene <- fread(here("data/pheno.prom.gene.csv"))
    promoter.peaks <- fread(here("data/driver_gene_list.csv"))
    gene_cna_df <- prom.gene %>%
        mutate(samp = gsub("-", "_", METABRIC_ID)) %>%
        select(samp, name, name3.chr, cna = cna.cat.zip) %>%
        left_join(samp_data %>% select(samp, digpath_lymph), by = "samp") %>%
        left_join(promoter.peaks %>% select(name = gene, type), by = "name") %>%
        left_join(get_all_features(), by = c("samp"))
    return(gene_cna_df)
}

cna_import <- function(){
    load(here("metadata/cna.final.RData"))
    as.data.frame(cna.chr) %>% 
        select(chrom = seqnames, start, end, CNA, samp=METABRIC_ID, nMajor, nMinor) %>% 
        mutate(samp = gsub('-', '_', samp)) %>% 
        left_join(samp_data %>% select(samp, cna_ploidy)) %>% 
        filter(!is.na(cna_ploidy)) %>% 
        as_tibble() %>% 
        mutate(cna_val = 2*CNA / cna_ploidy) %>% 
        fwrite(here('metadata/cna_all.tsv'))
}