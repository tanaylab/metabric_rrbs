# Files needed for the analysis

Due to the size of the METABRIC-RRBS dataset (~2.2TB full, 55GB only pileup), we generated a few smaller processed files to help reproduce the analysis. 
See scripts that generate those files at `raw-data.Rmd` and `pipeline.Rmd`.  

In addition, during the analysis some files are heavy to compute and therefore are cached for convenience. See below a list of all the files. 

In general - in order to run the analysis notebooks you would need to first download the processed files from https://metabric-rrbs.s3.eu-west-1.amazonaws.com/analysis_files.tar.gz. 

The analysis files bundle contains a `misha` db and additional processed files: 

## Misha DB

A [misha](https://github.com/tanaylab/misha) database of `hg19` is needed for the analysis. It contains the basic intervals such as gene annotations, together with additional tracks that are used in the analysis. You can see those tracks at the [configuration file](scripts/config/config.yaml).

## Additional files


### METABRIC data

- Promoter methylation: 
    all: `data/promoter_avg_meth.csv`
    filtered (coverage >= 20 in at least 70% of tumor samples and 70% of normal samples): `data/promoter_avg_meth_filt.csv`
    coverage: `data/promoter_cov.csv`    
- Non-promoter methylation (MSP1 fragments, not on promoter and at least 10 bp frmom exon):
    all: `data/genomic_msp1_avg_meth.csv`
    filtered (coverage >= 20 in at least 70% of tumor samples and 70% of normal samples): `data/genomic_msp1_avg_meth_filt.csv`
- Expression `data/expression_matrix.csv`
    See raw-data.Rmd for alternative promoter choice. 
- Mean Epipolymorhism per locus: `data/loci_epipoly_mean.tsv`. See the Epipolymorhism notebook for the code that generated it. 
- Copy Number Abberations `data/cna.tsv`
- Mutations: `data/mutations.tsv`
- Survival: `data/survival.tsv`
- Samples metadata: `data/samp_data.csv`
### QC files

- `data/sample_qc.csv`: per sample coverage statistics.
- `data/sample_coverage_dist.tsv`: distribution of CpG coverage per sample. 
- `data/cov_cpgs.tsv`: CpGs that are covered by at least 5 reads in half or more of the samples. 
- `data/cpg_cov_marginal.tsv`: Marginal coverage per CpG. 
- `data/sample_tot_meth_calls.tsv`: Total methylation calls per sample. 
- `data/sample_tot_meth_calls_promoters.tsv`: Total methylation calls on promoters per sample. 
- `data/well_covered_msp1_frags.tsv`: coordinates of msp1 fragments that are covered >= 20 in at least 70% of tumor samples and 70% of normal samples
-  `data/samp_genomic_meth.tsv`: average genomic methylation per sample. 

### Notable generated files

- `data/all_norm_meth.tsv`: TME normalized methylation from ER+/ER- and normal samples for both promoters and genomic regions. 
- `data/TME_features.tsv`: Immune and CAF features per sample (expression and methylation). 
- `data/features_loci_cors.tsv`: Correlation between genomic/promoter loci and the epigenomic features. 
- `data/data/loci_annot_epigenomic_features.tsv`: Epigenomic scores for every genomic/promoter locus in the genome. 
- `data/epigenomic_features.tsv`: epignomic scores for each samples. 
- `data/epigenomic_features_raw_meth.tsv`: raw methylation in epignomic scores regions for each samples. 
- `data/all_meth_summary.tsv"`: Average methylation of all ER+/ER-/normal samples per locus. 

### Spatial methylation distribution 

#### Full chromosomes

`data/tor_clock_chrom_trace_chr1_10000.tsv`
`data/tor_clock_chrom_trace_chr1_100000.tsv`
`data/tor_clock_chrom_trace_chr10_10000.tsv`

see [scripts/clock/chromosomal-traces.R] `get_tor_clock_chrom_trace`. We calculate average methylation in genomic bins for groups of METABRIC samples loss clock (top and bottom 30%).

#### Examples for regulation _in-cis_

`data/cis_promoter_examples_cg_meth.tsv`: average methylation around promoters shown in Figure 3E. 
`data/cis_genomic_examples_cg_meth.tsv`: average methylation around promoters shown in Figure 3H. 

### Copy number abberations

`data/pheno.prom.gene.csv`: CNA status per sample per promoter, together with the epigenomic features. 
`data/driver_gene_list.csv`: a list of driver and tumor suppressor genes. 

### Genomic annotations

See definition at [scripts/init_metadata.R](scripts/init_metadata.R), `define_genomic_regions`. 

- `data/k4me3_peaks_intervals.csv`
- `data/promoter_intervals.csv`
- `data/k27ac_intervals.csv`
- `data/enhancer_intervals.csv`
- `data/enhancer_intervals_tumors.csv`
- `data/k27me3_intervals.csv`
- `data/gene_tss.tsv`
- `data/gene_tss_coords.tsv`

- `data/genes_annot.csv`: manual annotation of genes to functional groups (Cell Cycle, Embryonic TF and other).

### External resources 

- `data/phenoAge.tsv`: phenoAge CpGs from PMID: 29676998
- `data/pheno_age_score.tsv`: phenoAge score for METABRIC samples. 


