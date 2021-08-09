gen_dosage_comp_excel <- function(thresh = 0.1){    
    autosome_cna_meth_expr_4N <- get_autosome_dosage_comp_cands("4N")
    autosome_cna_meth_expr_3N <- get_autosome_dosage_comp_cands("3N")

    dosage_cands_gain_4N <- autosome_cna_meth_expr_4N %>% 
        filter(`>=4N` >= (`2N` + thresh) ) %>%         
        filter(`n_>=4N` >= 5, `n_2N` >= 5,  !is.na(name))  %>%
        select(chrom, start, end, name, ER, `meth_>=4N` = `>=4N`, `meth_2N` = `2N`, `n_>=4N`, n_2N, `expr_>=4N`, `expr_2N`, `n_expr_>=4N`, `n_expr_2N`) %>%
        arrange(-`n_>=4N`, ER, name)

    dosage_cands_gain_3N <- autosome_cna_meth_expr_3N %>% 
        filter(`>=3N` >= (`2N` + thresh) ) %>%         
        filter(`n_>=3N` >= 5, `n_2N` >= 5,  !is.na(name))  %>%
        select(chrom, start, end, name, ER, `meth_>=3N` = `>=3N`, `meth_2N` = `2N`, `n_>=3N`, n_2N, `expr_>=3N`, `expr_2N`, `n_expr_>=3N`, `n_expr_2N`) %>%
        arrange(-`n_>=3N`, ER, name)

    dosage_cands_loss <- autosome_cna_meth_expr_4N %>% 
        filter(`1N` <= (`2N` - thresh) ) %>%         
        filter(`n_1N` >= 5, `n_2N` >= 5,  !is.na(name))  %>%
        select(chrom, start, end, name, ER, `meth_1N` = `1N`, `meth_2N` = `2N`, `n_1N`, n_2N, `expr_1N`, `expr_2N`, `n_expr_1N`, `n_expr_2N`) %>%
        arrange(-`n_1N`, ER, name)


    writexl::write_xlsx(x = list(`Gain 3N` = dosage_cands_gain_3N, `Amplification 4N` = dosage_cands_gain_4N, Loss = dosage_cands_loss), path = here("export/S9 - Dosage Compensation Candidates.xlsx"))
}