plot_score_clinical_annotation <- function(feat, label, feats = get_all_features(), samp_md=samp_data, ER="ER+", width = 1) {    
    p_grade <- feats %>%
        left_join(samp_md %>% select(samp, grade), by = "samp") %>%
        mutate(grade = ifelse(ER == "normal", "normal", grade)) %>%
        filter(!is.na(grade), ER %in% !! ER) %>%
        mutate(            
            grade = factor(grade, levels = c("normal", "1", "2", "3"))
        ) %>%
        ggplot(aes_string(x = feat, color = "grade")) +
        geom_density(size = width) +
        scale_color_manual(values = c("normal" = "gray", "1" = "darkblue", "2" = "red", "3" = "orange")) +
        theme(aspect.ratio = 1) +
        xlab(label) +
        ylab("Density") +
        guides(color = "none")
        
    p_stage <- feats %>% 
        left_join(samp_md %>% select(samp, stage), by = "samp") %>% 
        mutate(stage = ifelse(stage %in% c(0, "DCIS", 1), "0-1", stage)) %>% 
        mutate(stage = ifelse(ER == "normal", "normal", stage)) %>% 
        filter(!is.na(stage)) %>% 
        mutate(stage = factor(stage)) %>% 
        filter(ER %in% !!ER) %>%
        ggplot(aes_string(x = feat, color = "stage")) +
        geom_density(size = width) +
        scale_color_manual(values = c("normal" = "gray", "0-1" = "black", "2" = "blue", "3" = "red", "4" = "orange")) +
        theme(aspect.ratio = 0.7) +
        xlab(label) +
        ylab("Density") +
        guides(color = "none")

    p_cell <- feats %>%
        left_join(samp_md %>% select(samp, cna_cellularity), by = "samp") %>%
        filter(!is.na(cna_cellularity), ER %in% !! ER) %>%        
        ggplot(aes_string(x = feat, color = "cna_cellularity")) +
        geom_density(size = width) +
        scale_color_manual(values = c("black", "blue", "red", "orange", "yellow")) +
        theme(aspect.ratio = 0.7) +
        xlab(label) +
        ylab("Density") +
        guides(color = "none")    

    return(list(grade = p_grade, stage = p_stage, cellularity = p_cell))
}

plot_score_feats_boxp <- function(data, xlab, feats_tidy, nbins = 5) {
    p <- data %>%
        mutate(clin_feat = as.numeric(clin_feat), clin_feat = cut(clin_feat, quantile(clin_feat, 0:nbins / nbins, na.rm = TRUE), labels = 1:nbins, include.lowest = TRUE)) %>%
        left_join(feats_tidy, by = "samp") %>%
        na.omit() %>%
        filter(ER != "normal") %>%
        mutate(feat = factor(feat, levels = c("caf", "immune", "clock", "MG", "ML"))) %>% 
        ggplot(aes(x = clin_feat, y = score, fill = ER, group = clin_feat)) +        
        geom_boxplot(width = 0.5, fatten=0.5, outlier.size = 0.1) +
        scale_fill_manual(values = annot_colors$ER1, guide = FALSE) +
        facet_grid(feat ~ ER, scales = "free_y") +
        ylab("") +
        xlab(xlab) +               
        scale_x_discrete(breaks = c(1, nbins), labels = c("low", "high")) +
        theme(aspect.ratio = 1)        
    return(p)
}