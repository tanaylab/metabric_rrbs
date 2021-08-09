install_and_load_deps <- function(packages) {
    if (!require("BiocManager", character.only = TRUE)) {
        install.packages("BiocManager")
    }

    for (pkg in packages) {
        pkg_name <- gsub("tanaylab/", "", pkg)
        if (!require(pkg_name, character.only = TRUE)) {
            BiocManager::install(pkg)
            library(pkg_name, character.only = TRUE)
        }
    }
}


source_files <- function(scripts_dir) {
    src_files <- list.files(scripts_dir, pattern = "*.[rR]$", full.names = TRUE, recursive = TRUE)
    src_files <- src_files[src_files != glue("{scripts_dir}/init.R")]
    walk(src_files, source)
}


###############################################################
expr_intervs_to_mat <- function(df){
    df %>% select(-(chrom:end), -name3.chr) %>% column_to_rownames("name") %>% as.matrix()    
}

add_ER <- function(df){
    df %>% 
        left_join(samp_data %>% select(samp, ER=ER1), by="samp") %>% 
        mutate(ER = factor(ER, levels=c("ER+", "ER-", "normal")))
}



