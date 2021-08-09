get_expression_features_cors <- function() {    
    {
        expr_m <- fread(here("data/expression_matrix.csv")) %>% select(-any_of(c("chrom", "start", "end", "name3.chr")))

        expr_mat <- expr_m %>%
            as.data.frame() %>%
            column_to_rownames("name")
        f <- rowSums(!is.na(expr_mat)) > 0
        expr_mat <- expr_mat[f, ]

        feats <- get_all_features()

        cor_df <- plyr::ddply(feats, "ER", function(x) {
            samples <- reduce(list(samp_data$samp, colnames(expr_mat), x$samp), intersect)
            feats_mat <- x %>%
                select(-ER) %>%
                as.data.frame() %>%
                column_to_rownames("samp") %>%
                as.matrix()
            cm <- tgs_cor(t(expr_mat[, samples]), feats_mat[samples, ], pairwise.complete.obs = TRUE)
            cm <- cm %>%
                as.data.frame() %>%
                rownames_to_column("name") %>%
                as_tibble()
            return(cm %>% mutate(ER = x$ER[1]))
        })
    } %cache_df% here("data/expr_feat_cors.tsv") %>% as_tibble()
}

get_gene_annots <- function(){
    gene_annots <- fread(here("data/genes_annot.csv")) %>%
        select(name = gene, type) %>%
        mutate(type = case_when(type == "CC" ~ "Cell Cycle", type == "EMB_TF" ~ "Embryonic TF", TRUE ~ "Other"))
    colors <- tribble(~type, ~color,
                 "Cell Cycle", "#008B45FF",
                 "Embryonic TF", "#3B4992FF",
                 "Other", "purple")
    gene_annots <- gene_annots %>% left_join(colors) %>% as_tibble()
    return(gene_annots)
}

get_gene_features_df <- function(genes){
    expr_df <- fread(here("data/expression_matrix.csv")) %>% select(-(chrom:end), -name3.chr) %>% filter(name %in% genes) %>% gather("samp", "expr", -name)
    feats <- get_all_features()
    return(expr_df %>% left_join(feats, by = "samp") %>% as_tibble())
}