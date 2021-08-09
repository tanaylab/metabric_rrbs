get_all_features <- function(){
	fread(here("data/epigenomic_features.tsv")) %>% mutate(ML = -ML, clock = -clock, immune.meth = -immune.meth, caf.meth = -caf.meth) %>% as_tibble()
}

get_all_features_raw <- function(){
    fread(here("data/epigenomic_features_raw_meth.tsv")) %>% as_tibble()
}