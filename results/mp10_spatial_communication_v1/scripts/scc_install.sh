#!/usr/bin/env bash
# Installation only into this project's isolated Rlib. Official vendored source is retained.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 R_PROFILE_USER=/dev/null
export R_LIBS="$(pwd)/mp10_spatial_communication_v1/Rlib:/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5:/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4:/refdir/Rlib"
/usr/local/bin/Rscript --vanilla -e 'install.packages(c("ks","geometry"),lib="mp10_spatial_communication_v1/Rlib",repos="https://cran.rstudio.com")'
/usr/local/bin/R CMD INSTALL --library=mp10_spatial_communication_v1/Rlib mp10_spatial_communication_v1/vendor/KlugerLab-ALRA-*
/usr/local/bin/R CMD INSTALL --library=mp10_spatial_communication_v1/Rlib mp10_spatial_communication_v1/vendor/MERINGUE-source
/usr/local/bin/R CMD INSTALL --library=mp10_spatial_communication_v1/Rlib mp10_spatial_communication_v1/vendor/jinworks-SpatialCellChat-be88f30
