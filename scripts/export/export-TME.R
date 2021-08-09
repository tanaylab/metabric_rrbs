export_em_clust_expression <- function(){
    pos_expr_clust <-  readr::read_rds(here("data/ER_positive_norm_meth.rds"))$em_cross_clust$expr_clust
    neg_expr_clust <-  readr::read_rds(here("data/ER_negative_norm_meth.rds"))$em_cross_clust$expr_clust
    normal_expr_clust <-  readr::read_rds(here("data/normal_norm_meth.rds"))$em_cross_clust$expr_clust

    sheets <- list(
        `ER+` = pos_expr_clust %>% select(CE = clust, gene = name) %>% arrange(CE, gene),
        `ER-` = neg_expr_clust %>% select(CE = clust, gene = name) %>% arrange(CE, gene),
        `normal` = normal_expr_clust %>% select(CE = clust, gene = name) %>% arrange(CE, gene)
    )


    writexl::write_xlsx(x = sheets, path = here("export/Expression_Methylation_Correlations_TME.xlsx"))
}