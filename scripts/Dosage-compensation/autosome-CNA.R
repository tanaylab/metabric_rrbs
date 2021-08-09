get_autosome_loci <- function() {
    get_promoter_avg_meth() %>%
        select(chrom, start, end) %>%
        filter(chrom != "chrX", chrom != "chrY") %>%
        anti_join(get_xist_loci(), by = c("chrom", "start", "end"))
}

get_autosome_loci_meth <- function() {
    get_autosome_loci() %>%
        inner_join(get_promoter_avg_meth()) %>%
        gather("samp", "meth", -(chrom:end))
}

get_autosome_loci_expr <- function() {
    get_gene_expression_mat() %>%        
        inner_join(get_autosome_loci()) %>%
        gather("samp", "expr", -(chrom:name3.chr))
}


get_autosome_cna <- function() {        
        {
            autosome_cna <- get_autosome_loci() %>%
                gintervals.neighbors1(cna %>%
                    mutate(end = ifelse(start == end, start + 1, end)) %>%
                    filter(chrom != "chrX"), maxneighbors = nrow(samp_data)) %>%
                filter(dist == 0) %>%
                select(chrom, start, end, samp, cna = cna_round)

            autosome_cna <- autosome_cna %>%
                mutate(cna_grp = cut(cna, breaks=c(-1, 0,1,2,10), labels=c("0N", "1N", "2N", ">=3N"))) %>% 
                # mutate(cna_grp = cut(cna, breaks = c(0, 0.25, 0.75, 1.25, 10), labels = c("0N", "1N", "2N", ">=3N"), include.lowest = T)) %>%
                filter(cna_grp != "0N") %>%
                left_join(samp_data %>% select(samp, ER = ER1))
            autosome_cna 
        } %cache_df% here("data/autosome_cna.tsv") %>% as_tibble()
}

get_autosome_expr_cna <- function() {        
        {
            auto_df <- get_autosome_cna() %>%
                inner_join(get_autosome_loci_expr()) %>%
                group_by(chrom, start, end, ER, cna_grp) %>%
                summarise(expr = mean(expr, na.rm = TRUE)) %>%
                ungroup() %>%
                spread(cna_grp, expr)
            auto_df 
        } %cache_df% here("data/autosome_expr_cna.tsv") %>% as_tibble()
}

get_autosome_meth_cna <- function() {
    {    
        auto_df <- get_autosome_cna() %>%
            inner_join(get_autosome_loci_meth()) %>%
            group_by(chrom, start, end, ER, cna_grp) %>%
            summarise(meth = mean(meth, na.rm = TRUE)) %>%
            ungroup() %>%
            spread(cna_grp, meth)
        auto_df 
    } %cache_df% here("data/autosome_meth_cna.tsv") %>% as_tibble()
}

get_autosome_num_cna_meth <- function(){    
    {
        num_cna_meth <- get_autosome_cna() %>% 
            inner_join(get_autosome_loci_meth()) %>% 
            filter(!is.na(meth)) %>% 
            count(chrom, start, end, ER, cna_grp) %>% 
            mutate(cna_grp = paste0("n_", cna_grp)) %>% 
            spread(cna_grp, n)
        num_cna_meth
    } %cache_df% here("data/autosome_num_cna_meth.tsv") %>% as_tibble()
}

get_autosome_num_cna_expr <- function(){
    
    {
        num_cna_expr <- get_autosome_cna() %>% 
            inner_join(get_autosome_loci_expr()) %>% 
            filter(!is.na(expr)) %>% 
            count(chrom, start, end, ER, cna_grp) %>% 
            mutate(cna_grp = paste0("n_expr_", cna_grp)) %>% 
            spread(cna_grp, n)

        num_cna_expr  
    } %cache_df% here("data/autosome_num_cna_expr.tsv") %>% as_tibble()
}

get_autosome_cna_meth_expr <- function(){
    {
        num_cna_meth <- get_autosome_num_cna_meth()
        num_cna_expr <- get_autosome_num_cna_expr()
        cna_expr <- get_autosome_expr_cna() %>% gather("cna", "expr", -(chrom:ER)) %>% mutate(cna = paste0("expr_", cna)) %>% spread(cna, expr)
    
        autosome_cna_meth_expr <- get_autosome_meth_cna() %>% 
            left_join(num_cna_meth) %>% 
            left_join(cna_expr) %>%
            left_join(num_cna_expr) %>% 
            left_join(get_gene_expression_mat() %>% select(chrom:end, name))

        autosome_cna_meth_expr
    } %cache_df% here("data/autosome_cna_meth_expr.tsv")
}

get_autosome_dosage_comp_cands <- function(cna_type){    
    {
        autosome_cna <- get_autosome_loci() %>%
            gintervals.neighbors1(cna %>%
                mutate(end = ifelse(start == end, start + 1, end)) %>%
                filter(chrom != "chrX"), maxneighbors = nrow(samp_data)) %>%
            filter(dist == 0) %>%
            select(chrom, start, end, samp, cna = cna_round)

        
        if (cna_type == "4N"){            
            autosome_cna <- autosome_cna %>%
                mutate(cna_grp = cut(cna, breaks=c(-1, 0,1,2,3,10), labels=c("0N", "1N", "2N", "3N", ">=4N"))) %>%             
                filter(cna_grp != "0N") %>%
                left_join(samp_data %>% select(samp, ER = ER1))
        } else if (cna_type == "3N"){            
            autosome_cna <- autosome_cna %>%
                mutate(cna_grp = cut(cna, breaks=c(-1, 0,1,2,10), labels=c("0N", "1N", "2N", ">=3N"))) %>%             
                filter(cna_grp != "0N") %>%
                left_join(samp_data %>% select(samp, ER = ER1))
        } else {
            stop("cna_type can be only 4N or 3N")
        }
        

        num_cna_meth <- autosome_cna %>% 
            inner_join(get_autosome_loci_meth()) %>% 
            filter(!is.na(meth)) %>% 
            count(chrom, start, end, ER, cna_grp) %>% 
            mutate(cna_grp = paste0("n_", cna_grp)) %>% 
            spread(cna_grp, n)

        num_cna_expr <- autosome_cna %>% 
            inner_join(get_autosome_loci_expr()) %>% 
            filter(!is.na(expr)) %>% 
            count(chrom, start, end, ER, cna_grp) %>% 
            mutate(cna_grp = paste0("n_expr_", cna_grp)) %>% 
            spread(cna_grp, n)

        cna_expr <- autosome_cna %>%
            inner_join(get_autosome_loci_expr()) %>%
            group_by(chrom, start, end, ER, cna_grp) %>%
            summarise(expr = mean(expr, na.rm = TRUE)) %>%
            ungroup() %>%
            spread(cna_grp, expr)
        cna_expr <-  cna_expr %>% gather("cna", "expr", -(chrom:ER)) %>% mutate(cna = paste0("expr_", cna)) %>% spread(cna, expr)

        cna_meth <- autosome_cna %>%
            inner_join(get_autosome_loci_meth()) %>%
            group_by(chrom, start, end, ER, cna_grp) %>%
            summarise(meth = mean(meth, na.rm = TRUE)) %>%
            ungroup() %>%
            spread(cna_grp, meth)

        autosome_cna_meth_expr <- cna_meth %>% 
            left_join(num_cna_meth) %>% 
            left_join(cna_expr) %>%
            left_join(num_cna_expr) %>% 
            left_join(get_gene_expression_mat() %>% select(chrom:end, name))

        autosome_cna_meth_expr
    } %cache_df% here(glue("data/autosome_cna_meth_cands_{cna_type}.tsv")) %>% as_tibble()
}