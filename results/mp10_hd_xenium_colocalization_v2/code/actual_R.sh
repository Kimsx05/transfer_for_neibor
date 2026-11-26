#!/bin/bash
set -eu
export R_HOME=/usr/local/lib/R R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.4" R_ENVIRON_USER=/dev/null R_PROFILE_USER=/dev/null R_DEFAULT_PACKAGES="datasets,utils,grDevices,graphics,stats,methods"
export LD_LIBRARY_PATH="/usr/local/lib/R/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu"
exec /usr/local/lib/R/bin/exec/R --vanilla --slave -f "$1" --args "${@:2}"
