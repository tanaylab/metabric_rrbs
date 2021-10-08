packages <- c(
    "here",
    "tanaylab/tglkmeans",
    "tanaylab/tgstat",
    "tanaylab/tgutil",
    "glue",
    "tanaylab/misha",
    "tanaylab/misha.ext",
    "tanaylab/gpatterns",
    "ComplexHeatmap",
    "tidyverse",
    "Matrix",
    "extrafont",
    "tanaylab/tgppt",
    "tanaylab/methylayer",
    "cowplot",
    "patchwork",
    "vegan",
    "ggpubr",
    "umap"
)

suppressPackageStartupMessages(library(here))
source(here("scripts/utils.R"))
suppressPackageStartupMessages(install_and_load_deps(packages))

theme_set(tgppt::theme_arial(6))

set.seed(17)

init_global_defs <- function() {    
    gsetroot(here("db", "trackdb"))
    options(gmax.data.size = 1e9)
    options(tgutil.verbose = FALSE)

    scripts_dir <<- here("scripts")
    config_dir <<- here("scripts", "config")
    data_dir <<- here("data")
    main_config_file <<- glue("{config_dir}/config.yaml")
    main_config <<- yaml::read_yaml(main_config_file)    
}

options(tgutil.cache = TRUE) 

init_global_defs()
source_files(scripts_dir)



# Load sample metadata
init_metadata(main_config_file)
init_colors(main_config$colors)

define_genomic_regions(main_config_file)

options(tgutil.cache = TRUE) # Set to FALSE in order to recreate all files
