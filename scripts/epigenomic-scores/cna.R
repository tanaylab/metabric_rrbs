calc_features_cna_pvals <- function(feats = c("immune", "caf", "clock", "MG", "ML")) {
        cna_df <- get_gene_cna_df() %>%
            filter(ER != "normal") %>%
            mutate(type = forcats::fct_explicit_na(type)) %>%
            gather("feat", "score", -(samp:ER)) %>% 
            filter(feat %in% feats)

        pvals_er <- cna_df %>%
            group_by(name, ER, type, feat) %>%
            summarise(
                pval_loss = wilcox.test(score[cna == "LOSSLOH"], score[cna %in% c("NEUT", "GAINAMPL")])$p.value,
                pval_gain = wilcox.test(score[cna == "GAINAMPL"], score[cna %in% c("NEUT", "LOSSLOH")])$p.value,
                n_loss = sum(cna == "LOSSLOH", na.rm = TRUE),
                n_gain = sum(cna == "GAINAMPL", na.rm = TRUE),
                n_neut = sum(cna == "NEUT", na.rm = TRUE)
            ) %>%
            ungroup() %>%
            mutate(n = n_loss + n_gain + n_neut, p_loss = n_loss / n, p_gain = n_gain / n, p_neut = n_neut / n)

        pvals_all <- cna_df %>%
            group_by(name, type, feat) %>%
            summarise(
                pval_loss = wilcox.test(score[cna == "LOSSLOH"], score[cna %in% c("NEUT", "GAINAMPL")])$p.value,
                pval_gain = wilcox.test(score[cna == "GAINAMPL"], score[cna %in% c("NEUT", "LOSSLOH")])$p.value,
                n_loss = sum(cna == "LOSSLOH", na.rm = TRUE),
                n_gain = sum(cna == "GAINAMPL", na.rm = TRUE),
                n_neut = sum(cna == "NEUT", na.rm = TRUE)
            ) %>%
            ungroup() %>%
            mutate(n = n_loss + n_gain + n_neut, p_loss = n_loss / n, p_gain = n_gain / n, p_neut = n_neut / n) %>%
            mutate(ER = "all")

        pvals <- bind_rows(pvals_er, pvals_all)
        
        pvals     
}
