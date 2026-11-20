#!/bin/bash
set -eu
pkg=rust-proseg-3.2.0-hec9b1f2_0.conda
mkdir -p .tools
if [ ! -x .tools/prebuilt/bin/proseg ]; then
 curl -fLsS --connect-timeout 8 --max-time 120 "https://conda.anaconda.org/bioconda/linux-64/$pkg" -o ".tools/$pkg"
 echo "796898b9712790168d13a976dd951d2ec64de54d2631dc703e6f649caae64534  .tools/$pkg" | sha256sum -c -
 cph x ".tools/$pkg" --dest .tools/prebuilt
fi
.tools/prebuilt/bin/proseg --version
