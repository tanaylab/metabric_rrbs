---
jupyter:
  jupytext:
    formats: ipynb,Rmd
    text_representation:
      extension: .Rmd
      format_name: rmarkdown
      format_version: '1.2'
      jupytext_version: 1.11.4
  kernelspec:
    display_name: R
    language: R
    name: ir
---

# Export supplementary tables 


```r
source(here::here("scripts/init.R"))
```

## Table 1: METABRIC samples metadata


METABRIC samples profiled in this study.  

Metadata fields are:


- samp
- patient
- batch
- track
- total_reads
- mapped_reads
- mapped_frac
- cg_num
- meth_calls
- global_avg_meth
- type
- age
- grade
- stage
- ER
- IHC
- iC10
- PAM50
- matched_normal
- matched_tumor
- digpath_lymph
- digpath_stromal
- IMC_Fibroblasts
- IMC_Lymphocytes
- ASCAT_cellularity
- ASCAT_ploidy
- MathScore
- log10_global_epm


```r
samp_qc <- fread(here("data/sample_qc.csv")) %>% as_tibble()
```


```r
tab <- samp_data %>% left_join(samp_qc) %>% select(samp, patient, batch, track, total_reads, mapped_reads, mapped_frac, cg_num, meth_calls, global_avg_meth, type, age, grade, stage, ER, IHC, iC10, PAM50, matched_normal, matched_tumor, digpath_lymph, digpath_stromal, IMC_Fibroblasts, IMC_Lymphocytes, ASCAT_cellularity, MathScore, log10_global_epm)
```

```
## Joining, by = "track"
```


```r
writexl::write_xlsx(tab, here("export/S1 - Sample Information.xlsx"))
```

## Table 2: Immune expression signature 


Table of genes that were used to define the Immune expression signature, separated by ER status.


```r
immune_genes_tab <- map2_dfr(c("ER_positive", "ER_negative", "normal"), c("ER+", "ER-", "normal"), ~
        tibble(gene = get_TME_genes(readr::read_rds(here(glue("data/{.x}_norm_meth.rds")))$em_cross_clust, caf_gene = NULL)) %>% mutate(ER = .y) ) %cache_df% here("data/immune_genes_by_er.tsv")

immune_genes_tab %>% count(ER)
```

```
##       ER   n
## 1    ER- 345
## 2    ER+ 195
## 3 normal 864
```

```r
head(immune_genes_tab)
```

```
##      gene  ER
## 1   ACAP1 ER+
## 2     ADA ER+
## 3   ADAM7 ER+
## 4    AIM2 ER+
## 5    AOAH ER+
## 6 APBB1IP ER+
```


```r
writexl::write_xlsx(immune_genes_tab, here("export/S2 - Immune expression signature.xlsx"))
```

## Table 3: CAF expression signature 


Table of genes that were used to define the CAF expression signature, separated by ER status.


```r
caf_genes_tab <- map2_dfr(c("ER_positive", "ER_negative", "normal"), c("ER+", "ER-", "normal"), ~
        tibble(gene = get_TME_genes(readr::read_rds(here(glue("data/{.x}_norm_meth.rds")))$em_cross_clust, immune_gene = NULL)) %>% mutate(ER = .y) ) %cache_df% here("data/caf_genes_by_er.tsv")

caf_genes_tab %>% count(ER)
```

```
##       ER   n
## 1    ER- 360
## 2    ER+ 207
## 3 normal 592
```

```r
head(caf_genes_tab)
```

```
##      gene  ER
## 1   ABCA6 ER+
## 2   ABCA8 ER+
## 3   ABCB1 ER+
## 4  ACVRL1 ER+
## 5 ADAMTS9 ER+
## 6 ALDH1A2 ER+
```


```r
writexl::write_xlsx(caf_genes_tab, here("export/S3 - CAF expression signature.xlsx"))
```

## Table 4: Expression-Methylation correlation


Expression-Methylation correlation tables. 
Shown are pairs of tables of the following datasets:
1. Raw promoter methylation (Extended Data Fig. 2a,b).
2. Immune-CAF normalized methylation for loci that were correlated with MG/ML epigenomic instability (Fig. 2f).
3. X chromosome immune-CAF normalized methylation (Extended Data Fig 9a). 

Each dataset is represented by an expression table, showing the mean correlation of each one of 30/32 methylation clusters in each gene, and a methylation table showing the mean correlation of every locus to 30/32 expression clusters. 


```r
em_raw <- list(
        `ER+` = readr::read_rds(here("data/ER_positive_norm_meth.rds"))$em_cross_clust,
        `ER-` = readr::read_rds(here("data/ER_negative_norm_meth.rds"))$em_cross_clust,
        `normal` = readr::read_rds(here("data/normal_norm_meth.rds"))$em_cross_clust
)
```


```r
em_raw <- map(em_raw, parse_em_cors)
```


```r
norm_em <- readr::read_rds(here("data/MG_ML_em_cross_cor_clust.rds"))
norm_em <- parse_em_cors(norm_em)
```


```r
x_em <- readr::read_rds(here("data/X_er_positive_em_cross_cor_clust.rds"))
x_em <- parse_em_cors(x_em)
```


```r
sheets <- list(
    "Raw promoters ER+ (Expression)" = em_raw[["ER+"]]$expr_tab, 
    "Raw promoters ER+ (Methylation)" = em_raw[["ER+"]]$meth_tab,
    "Raw promoters ER- (Expression)" = em_raw[["ER-"]]$expr_tab, 
    "Raw promoters ER- (Methylation)" = em_raw[["ER-"]]$meth_tab, 
    "MG,ML loci ER+ (Expression)" = norm_em$expr_tab,
    "MG,ML loci ER+ (Methylation)" = norm_em$meth_tab,
    "X chromosome ER+ (Expression)" = x_em$expr_tab, 
    "X chromosome ER+ (Methylation)" = x_em$meth_tab
)
```


```r
writexl::write_xlsx(sheets, here("export/S4 - Expression-methylation correlations.xlsx"))
```

## Table 5: Methylation scores


Table with CAF, Immune, Clock, MG and ML methylation scores per METABRIC sample.


```r
tab <- get_all_features()
```


```r
writexl::write_xlsx(tab, here("export/S5 - Methylation scores.xlsx"))
```

## Table 6: Gene expression correlation to epigenomic instability


Tables of genes that have an absolute expression correlation higher than 0.3 to MG/ML/Clock scores, together with their respective correlations.


```r
feat_gene_cors <- get_expression_features_cors()
cor_thresh <- 0.3
```


```r
tab <- feat_gene_cors %>% filter(abs(MG) >= cor_thresh | abs(ML) >= cor_thresh | abs(clock) > cor_thresh) %>% mutate(ER = factor(ER, levels = c("ER+", "ER-", "normal"))) %>% select(name, ER,  clock.cor = clock, MG.cor = MG, ML.cor = ML, caf.cor = caf, immune.cor = immune) %>% arrange(ER, clock.cor, MG.cor, ML.cor)
```


```r
writexl::write_xlsx(tab, here("export/S6 - Gene expression correlation to methylation scores.xlsx"))
```

## Table 7: Methylation layers loci


Coordinates of loci that were used in order to calculate the CAF, Immune, Clock, MG and ML scores.  


```r
clust_df <- fread(here("data/ER_positive_loci_clust.tsv") ) %>% as_tibble()
```


```r
tab <- clust_df %>% rename(layer = clust) %>% filter(layer %in% c("ML", "MG", "clock")) %>% arrange(layer, chrom, start, end) 
head(tab)
```

```
## # A tibble: 6 x 4
##   chrom   start     end layer
## 1  chr1  134998  135215 clock
## 2  chr1  837884  838076 clock
## 3  chr1  907866  908002 clock
## 4  chr1  996514  996647 clock
## 5  chr1 1086481 1086732 clock
## 6  chr1 1381236 1381347 clock
```


```r
writexl::write_xlsx(tab, here("export/S7 - Methylation layers loci.xlsx"))
```

## Table 8: Cis regulation candidates 


Pairs of loci and genes that are candidates for cis regulation in different FDR thresholds. Promoters and genomic loci are shown in separate tables. 


```r
min_dist <- 5e5
min_tss_dist <- 200
```


```r
genomic_cis_cands <- bind_rows(
    fread(here("data/genomic_cis_cands_ER_positive.tsv")),
    fread(here("data/genomic_cis_cands_ER_negative.tsv")),
    fread(here("data/genomic_cis_cands_normal.tsv")))
head(genomic_cis_cands)
```

```
##   chrom  start    end type rank   gene        cor chrom_expr start_expr
## 1  chr1  10496  10587  obs    1 CT45A3 -0.2446525       chrX  134883487
## 2  chr1  10588  10639  obs    1  DSCR8 -0.1658337      chr21   39493544
## 3  chr1 134998 135215  obs    1 MAGEC2 -0.2381125       chrX  141293076
## 4  chr1 546168 546310  obs    1 MAGEA8 -0.2542376       chrX  149009940
## 5  chr1 565396 565791  obs    1 RAD51C -0.2153869      chr17   56769933
## 6  chr1 567121 567237  obs    1 TIMM23 -0.2271249      chr10   51623386
##    end_expr strand_expr dist n_obs n_shuff        fdr  ER
## 1 134883488           1   NA  2680      92 0.03432836 ER+
## 2  39493545           1   NA  2680      92 0.03432836 ER+
## 3 141293077          -1   NA  2680      92 0.03432836 ER+
## 4 149009941           1   NA  2680      92 0.03432836 ER+
## 5  56769934           1   NA  2680      92 0.03432836 ER+
## 6  51623387          -1   NA  2680      92 0.03432836 ER+
```

```r
dim(genomic_cis_cands)
```

```
## [1] 55616700       16
```


```r
cands_genomic <- genomic_cis_cands %>%
    filter(type == "obs", !is.na(dist), abs(dist) <= min_dist, abs(dist) >= min_tss_dist) %>%
    arrange(cor) %>% 
    as_tibble()
```


```r
source(here::here("scripts/init.R"))
```


```r
genomic_list <- list(
    `Genomic (best correlation)` = cands_genomic %>% filter(rank == 1),
    `Genomic (FDR < 0.05)` = cands_genomic %>% filter(fdr < 0.05), 
    `Genomic (FDR < 0.1)` = cands_genomic %>% filter(fdr < 0.1)
)

genomic_list <- map(genomic_list, ~ annotate_cis_cands(.x, sigma_meth = 2, sigma_expr = 2, meth_diff = 0.2, expr_diff = 1))
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```


```r
genomic_list <- map(genomic_list, ~ .x %>% 
    select(gene, 
           ER, 
           chrom, 
           start, 
           end,            
           distance = dist, 
           fdr, 
           rank, 
           cor, 
           mean_meth, 
           sd_meth, 
           mean_expr, 
           sd_expr, 
           normal_meth, 
           normal_meth_sd, 
           normal_expr, 
           normal_expr_sd, 
           n_hypometh, 
           n_induced, 
           n_hypermeth, 
           n_repressed, 
           N_considered, 
           n_hypometh_vs_normal,
           n_stable_vs_normal, 
           n_hypermeth_vs_normal, 
           n_repressed_vs_normal, 
           n_stable_expr_vs_normal, 
           n_induced_vs_normal) %>% 
                    arrange(cor) 
)
```


```r
cands_prom <- fread(here("data/promoter_cis_cands.tsv")) %>% as_tibble() 
```


```r
prom_list <- list(
    `Promoters (best correlation)` = cands_prom %>% filter(r == 1),
    `Promoters (FDR < 0.005)` = cands_prom %>% filter(fdr < 0.005), 
    `Promoters (FDR < 0.01)` = cands_prom %>% filter(fdr < 0.01),
    `Promoters (FDR < 0.05)` = cands_prom %>% filter(fdr < 0.05)
)
```


```r
prom_list <- map(prom_list, ~ annotate_cis_cands(.x %>% rename(gene = name), sigma_meth = 2, sigma_expr = 2, meth_diff = 0.2, expr_diff = 1))
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
```

```
## Joining, by = "gene"
```

```
## Joining, by = "coord"
## Joining, by = "coord"
```

```
## Joining, by = "gene"
## Joining, by = "gene"
```


```r
prom_list <- map(prom_list, ~ 
                .x %>% 
        mutate(diff = abs(best - kth)) %>% 
        select(gene, ER, chrom, start, end, fdr, rank = r, cor, diff_from_second_best = diff, mean_meth, sd_meth, mean_expr, sd_expr, normal_meth, normal_meth_sd, normal_expr, normal_expr_sd, n_hypometh, n_induced, n_hypermeth, n_repressed, N_considered, n_hypometh_vs_normal, n_stable_vs_normal, n_hypermeth_vs_normal, n_repressed_vs_normal, n_stable_expr_vs_normal, n_induced_vs_normal) %>%
                arrange(cor) )
```


```r
cis_cands <- c(prom_list, genomic_list)    
```


```r
writexl::write_xlsx(x = cis_cands, path = here("export/S8 - Cis Regulation Candidates.xlsx"))
```

## Table 9: Dosage compensation in autosomes 


Table of autosome promoters that show increase of at least 10% in methylation when amplified to at least 3N vs 2N (‘Gain 3N’ table), at least 4N vs 2N (‘Amplification 4N’). The table ‘Loss’ shows loci that showed decreased methylation of at least 10% when losing a copy (1N vs 2N). 


```r
gen_dosage_comp_excel()
```


```r
gc()
```

```
##              used    (Mb) gc trigger    (Mb)   max used    (Mb)
## Ncells    4173577   222.9    8029078   428.8    8029078   428.8
## Vcells 2635840254 20109.9 5129628008 39136.0 5101547905 38921.8
```
