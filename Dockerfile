# Note: the dockerfile is assumed to be built from the wiser project directory
FROM jupyter/datascience-notebook:r-4.1.1

# Install rpm dependencies
USER root
RUN apt-get update && apt-get install -y git-core libcurl4-openssl-dev libgit2-dev libicu-dev libssl-dev libxml2-dev make pandoc pandoc-citeproc zlib1g-dev libgtk2.0-dev libhiredis-dev libcairo2-dev libxt-dev xvfb xauth xfonts-base vim && rm -rf /var/lib/apt/lists/*
USER ${NB_UID}

# Install BiocManager and reticulate
RUN R -e 'install.packages("BiocManager", repos = "http://cran.us.r-project.org")'
RUN R -e 'install.packages("reticulate", repos = "http://cran.us.r-project.org")'

# Install bioc dependencies
RUN R -e 'BiocManager::install(c("here", "tgstat", "glue", "tidyverse", "Matrix", "extrafont", "cowplot", "patchwork", "vegan", "bookdown", "ggpubr", "umap", "survminer"))'

# Install github dependencies
RUN R -e 'remotes::install_github("tanaylab/tglkmeans")'
RUN R -e 'remotes::install_github("tanaylab/tgutil")'
RUN R -e 'remotes::install_github("tanaylab/misha")'
RUN R -e 'remotes::install_github("tanaylab/misha.ext")'
RUN R -e 'remotes::install_github("tanaylab/gpatterns")'
RUN R -e 'remotes::install_github("tanaylab/tgppt")'
RUN R -e 'remotes::install_github("tanaylab/methylayer")'
RUN R -e 'BiocManager::install("ComplexHeatmap")'

ENV JUPYTER_ENABLE_LAB=yes

USER root

COPY ./scripts /workdir/scripts
COPY ./analysis /workdir/analysis

RUN chown -R ${NB_UID} /workdir
RUN chgrp -R ${NB_GID} /workdir

USER ${NB_UID}

RUN touch /workdir/.here

WORKDIR /workdir

