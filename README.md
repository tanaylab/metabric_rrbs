
# DNA methylation landscapes of 1538 breast cancers reveal a replication-linked clock, epigenomic instability and cis-regulation

<!-- badges: start -->
<!-- badges: end -->

This is the code that generates the figures for the [METABRIC RRBS paper](https://doi.org/10.1038/s41467-021-25661-w). The code is splitted to jupyter notebooks that can be found under the analysis folder. 

## Run the notebooks

See below for docker instructions. 

Due to the size of the METABRIC-RRBS dataset (~2.2TB full, 55GB pileup alone), we generated a few smaller processed files to help reproduce the analysis. Even this bundle is quite large (~50GB), and you can download it from: 

https://metabric-rrbs.s3.eu-west-1.amazonaws.com/analysis_files.tar.gz

The above file contains a folder named `db` and a folder named `data`. Please copy them to the main folder before running the notebooks. 

See `files.md` for a description of these files and `pipeline.Rmd` for the code that generated them. 

You can also run all the notebooks at once by running (inside the analysis folder):

```r
bookdown::render_book("index.Rmd", new_session = TRUE, output = bookdown::gitbook())
```

## Dependencies

The initialization script (`scripts/init.R`) installs automatically the necessary R packages to run the notebooks. 

## Notebook order 

It is recommended to run the notebooks in the following order: 

1. coverage-stats
2. TME
3. Epigenomic-scores
4. Loss-clock
5. epigenomic-instability
6. CNA
7. mutations
8. survival
9. cis
10. dosage-compensation

Running in a different order might work for some of the notebooks, but it might fail for others due to dependencies to data that was generated in previous notebooks. 

## Find a specific figure

`figures-key.md` lists where each figure in the paper was generated. 

## Docker 

We also provide a docker image which contains all the needed dependencies, to use it run:

Download the analysis files: 

```bash
wget https://metabric-rrbs.s3.eu-west-1.amazonaws.com/analysis_files.tar.gz
tar xvzf analysis_files.tar.gz
```

Change permissions for the analysis files:

```bash
chmod a+wx data/
chmod a+r -R data/
chmod a+rx -R db/
```

Run the container:

```bash
docker run -v $(pwd)/db:/workdir/db -v $(pwd)/data:/workdir/data -ti -p 8888:8888 tanaylab/metabric-rrbs
```

Connect to the jupyter server running at port 8888.




