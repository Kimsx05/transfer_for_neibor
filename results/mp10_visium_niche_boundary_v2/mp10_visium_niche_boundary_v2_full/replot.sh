#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
version=${1:-primary}
bash run.sh figures --version "$version"
if [[ -d "$HOME/.agents/skills/jins-routine-plotting" ]]; then bash composition.sh "$version"; fi
bash run.sh diagnostic_plots --version "$version"
bash run.sh core_maps --version "$version"
bash run.sh index --version "$version"

