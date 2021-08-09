pat_covs_marginal <- function(tracks, pat_len=5){
    add_covs_chrom <- function(tracks, chr){
        print(chr)
        chrom_len <- gintervals.all() %>% filter(chrom == chr) %>% .$end
        poss <- rep(0, chrom_len)        
        for (track in tracks){   
            gvtrack.create('vt', glue("{track}.pat{pat_len}"))
            gvtrack.array.slice(vtrack='vt', func='sum')
            a <- gextract('vt', intervals=gintervals.all() %>% filter(chrom == chr), colnames='cov') %>% as_tibble()       
            poss[a$start] <- poss[a$start] + a$cov
        }
        tibble(start=which(poss != 0), chrom=chr) %>% mutate(cov = poss[start]) %>% select(chrom, start, cov)    
    }
    
    res <- as.character(gintervals.all()$chrom) %>% plyr::alply(1, function(x) add_covs_chrom(tracks, x), .parallel=TRUE)
    pat_cov <- map_df(res, ~ .x) %>% as_tibble() %>% mutate(end = start + 1) %>% select(chrom, start, end, cov)
    return(pat_cov)    
}

get_pat_covs <- function(recalc=FALSE){
	fn <- 'data/pat5_covs.csv'
	if (recalc || !file.exists(fn)){
		pat_cov <-  pat_covs_marginal(tracks=samp_data$track, pat_len=5)	
		fwrite(pat_cov, fn)	
	} 

	fread(fn) %>% as_tibble()
	
}

get_pat_space <- function(recalc=FALSE){
	fn <- 'data/msp1_pat_space.csv'
	if (recalc || !file.exists(fn)){
		pat5_covs <- get_pat_covs(recalc=recalc)

		msp_frags <- gintervals.load('intervs.msp1.fid') %>% as_tibble()
		pat_space <- pat5_covs %>% 
			gintervals.neighbors1(msp_frags) %>% 
			filter(dist == 0) %>% 
			group_by(chrom1, start1, end1, FID) %>% 
			filter(cov == max(cov, na.rm=TRUE)) %>% 
			ungroup() %>% 
			select(chrom:end, fid=FID)

		fwrite(pat_space, fn)
	} 

	pat_space <- fread(fn) %>% as_tibble()	

	return(pat_space)
}

track_to_pats <- function(track, pat_space){
	a <- gtrack.array.extract(glue('{track}.pat5'), intervals=pat_space) %>%
		select(-intervalID) %>% 
		left_join(pat_space) %>% 
		as_tibble()
	colnames(a) <- gsub('^pat_', '', colnames(a))
	a %>% 
		gather('pattern', 'n', -chrom, -start, -end, -fid) %>% 
		select(fid, pattern, n) %>% 
		uncount(n)
}

create_patterns_track <- function(track, pat_space=get_pat_space()){
	print(track)
	patterns_tab <- track_to_pats(track, pat_space)
	gpatterns.create_patterns_track(track, patterns_tab=patterns_tab, pat_space=pat_space, add_read_id=FALSE, description='xxx')
}

create_downsampled_track <- function(track, dsn=30){
	devtools::load_all('/home/aviezerl/repo/gpatterns')
	gpatterns.create_downsampled_track(track, dsn, add_read_id=FALSE, description='xxx')
}

create_all_patterns_tracks <- function(){
	commands <- glue('create_patterns_track("{samp_data$track}")')

	res <- gcluster.run2(command_list = commands, io_saturation=TRUE)	
	gdb.reload()
}

create_all_downsampled_tracks <- function(dsn=30){
	commands <- glue('create_downsampled_track("{samp_data$track}", dsn={dsn})')
	
	res <- gcluster.run2(command_list = commands, io_saturation=TRUE)	
	gdb.reload()
}

