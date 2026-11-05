FROM rocker/binder:4.4.1
LABEL maintainer='Francesco Cecere'

# Command-line tools for Module 1 (Ubuntu 22.04 packages: bwa 0.7.17, samtools/bcftools 1.13, fastp 0.20.1, bedtools 2.30).
USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends bwa samtools bcftools tabix fastp bedtools && \
    rm -rf /var/lib/apt/lists/*

COPY --chown=${NB_USER} . ${HOME}
USER ${NB_USER}

RUN R -e "options(repos = list(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/latest')); \
          pkgs <- c('tidyverse', 'jsonlite', 'DT', 'BiocManager', 'knitr', 'rmarkdown', 'R.utils'); \
          install.packages(pkgs); \
          BiocManager::install('maftools', update = FALSE, ask = FALSE); \
          all_pkgs <- c(pkgs, 'maftools'); \
          missing <- all_pkgs[!sapply(all_pkgs, requireNamespace, quietly = TRUE)]; \
          if (length(missing) > 0) stop('Failed to install: ', paste(missing, collapse = ', '))"
