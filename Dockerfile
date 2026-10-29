# processor-cycombine: cyCombine batch correction of arcsinh-transformed FCS files.
# R 4.6.1 (rocker/r-ver, Ubuntu 24.04), binary packages from Posit Public Package Manager
# pinned to one snapshot date, Bioconductor 3.23, cyCombine pinned to a GitHub commit.
FROM --platform=linux/amd64 rocker/r-ver:4.6.1@sha256:bafb3b9a72791eeeb810f0bc0472e74c7f4939b20b46c96cb7ebc6838b6a25ef

ARG P3M_SNAPSHOT=2026-10-26
ARG BIOC_VERSION=3.23
# biosurf/cyCombine main as of 2026-04-29 (merge of PR #73, DESCRIPTION Version 0.3.0)
ARG CYCOMBINE_SHA=fcd42dbbe92bfa6bffa6094bf9b4de536e71c249

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=en_US.UTF-8 \
    OMP_NUM_THREADS=1

# Runtime libraries the binary packages link against (XML for sva's annotate/genefilter chain).
RUN apt-get update \
 && apt-get install -y --no-install-recommends libxml2 libglpk40 libgmp10 \
 && rm -rf /var/lib/apt/lists/*

RUN Rscript -e "\
  options(repos = c(CRAN = 'https://p3m.dev/cran/__linux__/noble/${P3M_SNAPSHOT}'), \
          BioC_mirror = 'https://p3m.dev/bioconductor/__linux__/noble/${P3M_SNAPSHOT}', \
          Ncpus = 4); \
  install.packages(c('BiocManager', 'remotes')); \
  BiocManager::install(c('flowCore', 'sva'), version = '${BIOC_VERSION}', update = FALSE, ask = FALSE); \
  install.packages(c('magrittr', 'kohonen', 'tibble', 'dplyr', 'readxl', 'readr', 'stringr', \
                     'emdist', 'ggplot2', 'ggridges', 'cowplot', 'plyr', 'tidyr', 'jsonlite')); \
  remotes::install_github('biosurf/cyCombine@${CYCOMBINE_SHA}', dependencies = FALSE, upgrade = 'never', \
                          build_vignettes = FALSE); \
  for (p in c('flowCore', 'sva', 'kohonen', 'cyCombine', 'emdist', 'ggridges', 'cowplot', 'jsonlite')) \
    if (!requireNamespace(p, quietly = TRUE)) stop('package missing: ', p)" \
 && rm -rf /tmp/downloaded_packages /tmp/Rtmp* /root/.cache

WORKDIR /app
COPY VERSION /app/VERSION
COPY R/ /app/R/

# No network at run time; everything the run needs is in the image.
ENV INPUT_DIR=/data/input \
    OUTPUT_DIR=/data/output

ENTRYPOINT ["Rscript", "/app/R/main.R"]
