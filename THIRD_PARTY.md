# Third-party software

This app is Apache-2.0 (see `LICENSE`). The image bundles:

| Component | Version | Licence | Source | Role |
|---|---|---|---|---|
| cyCombine | 0.3.0, commit fcd42dbbe92bfa6bffa6094bf9b4de536e71c249 | MIT (LICENSE: © 2020 Lars Rønn Olsen) | https://github.com/biosurf/cyCombine | Normalisation, SOM, ComBat per node, EMD/MAD evaluation, density plots |
| sva | 3.60.0 | Artistic-2.0 | Bioconductor | ComBat |
| kohonen | 3.0.13 | GPL (>= 2) | CRAN | Self-organising map |
| emdist | 0.3-3 | MIT | CRAN | Earth mover's distance |
| flowCore | 2.24.0 | Artistic-2.0 | Bioconductor | FCS reading |
| cytolib | 2.24.0 | AGPL-3.0-only | Bioconductor | flowCore's C++ FCS backend |
| ggplot2, ggridges, cowplot | 4.0.3, 0.5.7, 1.2.0 | MIT, GPL-2, GPL-2 | CRAN | Density plots |
| jsonlite | 2.0.0 | MIT | CRAN | provenance.json |
| R | 4.6.1 | GPL-2 / GPL-3 | rocker/r-ver:4.6.1 | Runtime |
| Ubuntu 24.04 base | | various free licences | rocker/r-ver | OS |

## Licence of the image as a whole

Our R code is Apache-2.0. It loads the packages above into one R process, and several
are GPL (R itself, kohonen, ggridges, cowplot, mgcv and others in the table) and one is AGPL-3.0-only: **cytolib**, the C++ library behind
flowCore's FCS reader. The image as distributed is therefore a combined work under the
GNU GPL v3 / AGPL v3 terms (Apache-2.0 is compatible with GPLv3 in that direction). The
AGPL's network clause concerns software offered to users over a network; this app is a
batch container without a network interface. Source for every component is available
from the repositories listed (CRAN/Bioconductor sources at the versions shown; GitHub
packages at the pinned commits). R packages are redistributed under their own licences.

## Every R package in the image

R 4.6.1 (rocker/r-ver:4.6.1, Ubuntu 24.04), packages from Posit Public Package Manager
snapshot 2026-10-26 (CRAN) and Bioconductor 3.23. Generated from
`installed.packages()` in the built image (base packages omitted).

| Package | Version | Licence |
|---|---|---|
| annotate | 1.90.0 | Artistic-2.0 |
| AnnotationDbi | 1.74.0 | Artistic-2.0 |
| askpass | 1.2.1 | MIT + file LICENSE |
| BH | 1.90.0-1 | BSL-1.0 |
| Biobase | 2.72.0 | Artistic-2.0 |
| BiocGenerics | 0.58.1 | Artistic-2.0 |
| biocmake | 1.4.0 | MIT + file LICENSE |
| BiocManager | 1.30.27 | Artistic-2.0 |
| BiocParallel | 1.46.0 | GPL-2 / GPL-3 / BSL-1.0 |
| BiocVersion | 3.23.1 | Artistic-2.0 |
| Biostrings | 2.80.2 | Artistic-2.0 |
| bit | 4.6.0 | GPL-2 / GPL-3 |
| bit64 | 4.8.6 | GPL-2 / GPL-3 |
| blob | 1.3.0 | MIT + file LICENSE |
| boot | 1.3-32 | Unlimited |
| cachem | 1.1.0 | MIT + file LICENSE |
| cellranger | 1.1.0 | MIT + file LICENSE |
| class | 7.3-23 | GPL-2 / GPL-3 |
| cli | 3.6.6 | MIT + file LICENSE |
| clipr | 0.8.1 | GPL-3 |
| cluster | 2.1.8.2 | GPL (>= 2) |
| codetools | 0.2-20 | GPL |
| cowplot | 1.2.0 | GPL-2 |
| cpp11 | 0.5.5 | MIT + file LICENSE |
| crayon | 1.5.3 | MIT + file LICENSE |
| curl | 8.1.0 | MIT + file LICENSE |
| cyCombine | 0.3.0 | MIT + file LICENSE |
| cytolib | 2.24.0 | AGPL-3.0-only |
| DBI | 1.3.0 | LGPL (>= 2.1) |
| dir.expiry | 1.20.0 | GPL-3 |
| docopt | 0.7.2 | MIT + file LICENSE |
| dplyr | 1.2.1 | MIT + file LICENSE |
| edgeR | 4.10.7 | GPL (>=2) |
| emdist | 0.3-3 | MIT + file LICENSE |
| farver | 2.1.2 | MIT + file LICENSE |
| fastmap | 1.2.0 | MIT + file LICENSE |
| filelock | 1.0.3 | MIT + file LICENSE |
| flowCore | 2.24.0 | Artistic-2.0 |
| foreign | 0.8-91 | GPL (>= 2) |
| formatR | 1.14 | GPL |
| futile.logger | 1.4.9 | LGPL-3 |
| futile.options | 1.0.1 | LGPL-3 |
| genefilter | 1.94.0 | Artistic-2.0 |
| generics | 0.1.4 | MIT + file LICENSE |
| ggplot2 | 4.0.3 | MIT + file LICENSE |
| ggridges | 0.5.7 | GPL-2 / file LICENSE |
| glue | 1.8.1 | MIT + file LICENSE |
| gtable | 0.3.6 | MIT + file LICENSE |
| hms | 1.1.4 | MIT + file LICENSE |
| httr | 1.4.9 | MIT + file LICENSE |
| IRanges | 2.46.0 | Artistic-2.0 |
| isoband | 0.3.0 | MIT + file LICENSE |
| jsonlite | 2.0.0 | MIT + file LICENSE |
| KEGGREST | 1.52.2 | Artistic-2.0 |
| KernSmooth | 2.23-26 | Unlimited |
| kohonen | 3.0.13 | GPL (>= 2) |
| labeling | 0.4.3 | MIT + file LICENSE / Unlimited |
| lambda.r | 1.2.6 | LGPL-3 |
| lattice | 0.22-9 | GPL (>= 2) |
| lifecycle | 1.0.5 | MIT + file LICENSE |
| limma | 3.68.5 | GPL (>=2) |
| littler | 0.3.23 | GPL (>= 2) |
| locfit | 1.5-9.12 | GPL (>= 2) |
| magrittr | 2.0.5 | MIT + file LICENSE |
| MASS | 7.3-65 | GPL-2 / GPL-3 |
| Matrix | 1.7-5 | GPL (>= 2) / file LICENCE |
| MatrixGenerics | 1.24.0 | Artistic-2.0 |
| matrixStats | 1.5.0 | Artistic-2.0 |
| memoise | 2.0.1 | MIT + file LICENSE |
| mgcv | 1.9-4 | GPL (>= 2) |
| mime | 0.13 | GPL |
| nlme | 3.1-169 | GPL (>= 2) |
| nnet | 7.3-20 | GPL-2 / GPL-3 |
| openssl | 2.4.2 | MIT + file LICENSE |
| pillar | 1.11.1 | MIT + file LICENSE |
| pkgconfig | 2.0.3 | MIT + file LICENSE |
| plyr | 1.8.9 | MIT + file LICENSE |
| png | 0.1-9 | GPL-2 / GPL-3 |
| prettyunits | 1.2.0 | MIT + file LICENSE |
| progress | 1.2.3 | MIT + file LICENSE |
| purrr | 1.2.2 | MIT + file LICENSE |
| R6 | 2.6.1 | MIT + file LICENSE |
| RColorBrewer | 1.1-3 | Apache License 2.0 |
| Rcpp | 1.1.2 | GPL (>= 2) |
| readr | 2.2.0 | MIT + file LICENSE |
| readxl | 1.6.0 | MIT + file LICENSE |
| rematch | 2.0.0 | MIT + file LICENSE |
| remotes | 2.5.0 | MIT + file LICENSE |
| Rhdf5lib | 2.0.0 | Artistic-2.0 |
| rlang | 1.3.0 | MIT + file LICENSE |
| rpart | 4.1.27 | GPL-2 / GPL-3 |
| RProtoBufLib | 2.24.0 | BSD_3_clause |
| RSQLite | 3.53.3 | LGPL (>= 2.1) |
| S4Vectors | 0.50.3 | Artistic-2.0 |
| S7 | 0.2.2 | MIT + file LICENSE |
| scales | 1.4.0 | MIT + file LICENSE |
| Seqinfo | 1.2.0 | Artistic-2.0 |
| snow | 0.4-4 | GPL |
| spatial | 7.3-18 | GPL-2 / GPL-3 |
| statmod | 1.5.2 | GPL-2 / GPL-3 |
| stringi | 1.8.9 | file LICENSE |
| stringr | 1.6.0 | MIT + file LICENSE |
| survival | 3.8-6 | LGPL (>= 2) |
| sva | 3.60.0 | Artistic-2.0 |
| sys | 3.4.3 | MIT + file LICENSE |
| tibble | 3.3.1 | MIT + file LICENSE |
| tidyr | 1.3.2 | MIT + file LICENSE |
| tidyselect | 1.2.1 | MIT + file LICENSE |
| tzdb | 0.5.0 | MIT + file LICENSE |
| utf8 | 1.2.6 | Apache License (== 2.0) / file LICENSE |
| vctrs | 0.7.3 | MIT + file LICENSE |
| viridisLite | 0.4.3 | MIT + file LICENSE |
| vroom | 1.7.1 | MIT + file LICENSE |
| withr | 3.0.3 | MIT + file LICENSE |
| XML | 3.99-0.25 | BSD_3_clause + file LICENSE |
| xtable | 1.8-8 | GPL (>= 2) |
| XVector | 0.52.0 | Artistic-2.0 |
