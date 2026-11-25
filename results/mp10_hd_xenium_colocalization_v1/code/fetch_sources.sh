#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/run_env.sh"
for name in HALLMARK_CHOLESTEROL_HOMEOSTASIS REACTOME_CHOLESTEROL_BIOSYNTHESIS; do
 p="$OUT/shared/msigdb/$name.json"
 if [[ ! -s "$p" ]]; then curl -fLsS --retry 2 "https://www.gsea-msigdb.org/gsea/msigdb/human/download_geneset.jsp?fileType=json&geneSetName=$name" -o "$p.part"; mv "$p.part" "$p"; fi
done
(cd "$OUT/shared/msigdb" && sha256sum -c SHA256SUMS)
for pair in 9625559_A1R2 9625560_A1R4 9625561_A1R5 9625562_A1R6 9625563_A2R2; do
 for asset in cells.csv.gz cell_boundaries.parquet.gz; do
  gsm=GSM${pair%_*}; f=GSM${pair}_${asset}; p="$OUT/xenium/geo/$f"
  if [[ ! -s "$p" ]]; then curl -fLsS --retry 2 "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9625nnn/$gsm/suppl/$f" -o "$p.part"; gzip -t "$p.part"; mv "$p.part" "$p"; fi
 done
done

