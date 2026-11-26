options(Ncpus=1)
# Install R-only package in task-local library with R's installation engine in-process.
tools:::.install_packages(args=c("--no-test-load","--no-byte-compile","--no-help","--no-demo","--no-docs","--no-html","--library=code/Rlib","code/CRAWDAD"))
