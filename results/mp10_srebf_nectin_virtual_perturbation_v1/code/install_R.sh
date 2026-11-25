#!/usr/bin/env bash
set -euo pipefail
source code/env.sh
for item in enrichR_3.4 scTenifoldNet_1.4.3 scTenifoldKnk_2.0.0; do
 pkg="${item%_*}"; ver="${item##*_}"
 if ! Rscript --vanilla -e "quit(status=if(requireNamespace('$pkg',quietly=TRUE) && as.character(packageVersion('$pkg'))=='$ver') 0 else 1)"; then
   R CMD INSTALL --library="$PWD/cache/Rlib" "cache/$item.tar.gz"
 fi
done

