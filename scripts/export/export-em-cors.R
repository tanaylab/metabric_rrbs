parse_em_cors <- function(l){
    meth_tab <- tgs_matrix_tapply(l$em_cross, l$expr_clust$clust, mean, na.rm=TRUE) %>% 
        t() %>% 
        mat_to_intervs() %>% 
        right_join(l$meth_clust %>% column_to_rownames("locus") %>% mat_to_intervs(), by = c("chrom", "start", "end")) %>% 
        select(chrom, start, end, clust, everything()) %>% 
        as_tibble()
        
    colnames(meth_tab)[-1:-4] <- paste0(colnames(meth_tab)[-1:-4], "_EXPR")
    
    expr_tab <- tgs_matrix_tapply(t(l$em_cross), l$meth_clust$clust, mean, na.rm=TRUE) %>% 
        t() %>% 
        as.data.frame() %>% 
        rownames_to_column("name") %>% 
        right_join(l$expr_clust, by = "name") %>%
        select(name, clust, everything()) %>% 
        as_tibble()
        
    colnames(expr_tab)[-1:-2] <- paste0(colnames(expr_tab)[-1:-2], "_METH")
        
    return(list(expr_tab = expr_tab, meth_tab = meth_tab))
}
