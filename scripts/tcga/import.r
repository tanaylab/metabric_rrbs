import_breast_tcga <- function(){
    library(curatedTCGAData)
    library(MultiAssayExperiment)
    library(TCGAutils)

    tcga <- curatedTCGAData("BRCA", dry.run = FALSE)
    
    expr_list <- as.list(assays(experiments(tcga)) )

    import_meth_array_TCGA_data(expr_list$`BRCA_Methylation_methyl27-20160128`, here("data/TCGA_BRCA_27k_mat.tsv"), "TCGA.BRCA_27k", "27k data from BRCA 20160128")

    import_meth_array_TCGA_data(expr_list$`BRCA_Methylation_methyl450-20160128`, here("data/TCGA_BRCA_450k_mat.tsv"), "TCGA.BRCA_450k", "450k data from BRCA 20160128")
    
    import_TCGA_expr_data(expr_list$`BRCA_RNASeq2GeneNorm-20160128`, here("data/TGCA_BRCA_expr.rds"))
}

import_TCGA_expr_data <- function(expr_mat, ofn){
    samp_md <- tibble(barcode = colnames(expr_mat)) %>% mutate(biospec = map(barcode, TCGAbiospec)) %>% unnest(biospec) 
    samp_md <- samp_md %>% mutate(type = case_when(sample %in% c("01", "02") ~ "T", sample == "06" ~ "M", sample %in% "11" ~ "N")) %>% mutate(samp_id = gsub("-", "_", submitter_id), samp_id = paste0(samp_id, "_", type))
    m <- expr_mat
    colnames(m) <- samp_md$samp_id
    readr::write_rds(m, ofn)        
}

import_meth_array_TCGA_data <- function(meth_mat, file, track, description, min_n = 100){
    meth_mat <- as.matrix(meth_mat)
    class(meth_mat) <- 'numeric'
    loci <- tibble(id = rownames(meth_mat)) %>% left_join(gintervals.load("intervs.450k_27k.cpgs")) %>% unite('coord', chrom:end, remove = FALSE)
    meth_mat <- meth_mat[!is.na(loci$start), ]
    loci <- loci %>% filter(!is.na(start))
    rownames(meth_mat) <- loci$coord

    samp_md <- tibble(barcode = colnames(meth_mat)) %>% mutate(biospec = map(barcode, TCGAbiospec)) %>% unnest(biospec) 
    samp_md <- samp_md %>% mutate(type = case_when(sample %in% c("01", "02") ~ "T", sample == "06" ~ "M", sample %in% "11" ~ "N")) %>% mutate(samp_id = gsub("-", "_", submitter_id), samp_id = paste0(samp_id, "_", type))
    colnames(meth_mat) <- samp_md$samp_id

    f <- rowSums(!is.na(meth_mat)) >= min_n
    meth_mat <- meth_mat[f, ]
    meth_df <- as.data.frame(meth_mat) %>% rownames_to_column("coord") %>% separate(coord, c("chrom", "start", "end"), sep="_") %>% mutate(start = as.numeric(start), end = as.numeric(end))
    data.table::fwrite(meth_df, file, na = "nan", row.names = FALSE, quote = FALSE, sep = "\t", scipen=50)
    if (gtrack.exists(track)) {
        gtrack.rm(track, force=TRUE)
    }
    gtrack.array.import(track, description, file)
    gtrack.var.set(track=track, var="samp_md", value=as.data.frame(samp_md))
}

import_TCGA_surv_data <- function(){
    # From: "An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics"
    download.file("https://ars.els-cdn.com/content/image/1-s2.0-S0092867418302290-mmc1.xlsx", destfile = "data/Liu2018_survivalData.xlsx")
    surv_df <- readxl::read_xlsx("data/Liu2018_survivalData.xlsx", sheet = "TCGA-CDR")[, -1] %>% filter(type == "BRCA")
    coad_surv <- surv_df %>% mutate(id = gsub("-", "_", bcr_patient_barcode), y = OS.time / 365) %>% select(patient = id, y, death = OS)
    fwrite(coad_surv, here("data/TCGA_BRCA_survival.tsv"), sep = "\t")
}


import_TCGA_samp_data <- function(tcga){
    tcga_patient_data <- colData(tcga) %>% as_tibble() %>% select(patientID, ER = ER.Status, HER2 = HER2.Final.Status, PR = PR.Status, gender = Gender, PAM50 = PAM50.mRNA, age =  Age.at.Initial.Pathologic.Diagnosis, stage = AJCC.Stage)
    tcga_patient_data <- tcga_patient_data %>% 
        mutate(PR = case_when(PR == "Positive" ~ "PR+", PR == "Negative" ~ "PR-")) %>%
        mutate(ER = case_when(ER == "Positive" ~ "ER+", ER == "Negative" ~ "ER-")) %>%
        mutate(HER2 = case_when(HER2 == "Positive" ~ "HER2+", HER2 == "Negative" ~ "HER2-")) %>%
        mutate(IHC = case_when(
            PR == "PR-" & ER == "ER-" & HER2 == "HER2-" ~ "TNBC", 
            ER == "ER-" & HER2 == "HER2+" ~ "ER-HER2+", 
            ER == "ER+" & HER2 == "HER2-" ~ "ER+HER2-", 
            ER == "ER+" & HER2 == "HER2+" ~ "ER+HER2+"
        )) %>% 
        dplyr::rename(submitter_id = patientID)
    
    tcga_samp_data <- bind_rows(
        gtrack.var.get("TCGA.BRCA_450k", "samp_md"),
        gtrack.var.get("TCGA.BRCA_27k", "samp_md")
    ) %>% as_tibble()
    
    tcga_samp_data <- tcga_samp_data %>% left_join(tcga_patient_data) %>% select(samp_id, everything()) %>% mutate(ER = ifelse(type == "N", "normal", ER)) %>% mutate(patient = gsub("-", "_", submitter_id))
    fwrite(tcga_samp_data, here("data/tcga_samp_data.csv"))
}