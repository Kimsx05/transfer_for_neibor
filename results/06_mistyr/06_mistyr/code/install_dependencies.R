# Reproducible optional installer. Writes only this module's private R library.
args<-commandArgs(TRUE);root<-normalizePath(args[1]);lib<-file.path(root,'runtime/Rlib')
dir.create(lib,recursive=TRUE,showWarnings=FALSE);.libPaths(c(lib,.libPaths()))
archive<-file.path(root,'sources/mistyR_complete.tar.gz')
expected<-'8e61c5a5c17fdd2ec1a8f8bf56148c054c90fa9dafbcea36c4ed9a6f69e23543'
stopifnot(file.exists(archive),digest::digest(file=archive,algo='sha256')==expected)
if(!requireNamespace('distances',quietly=TRUE))install.packages('distances',lib=lib,repos='https://cloud.r-project.org',Ncpus=1)
install.packages(archive,lib=lib,repos=NULL,type='source',Ncpus=1)
stopifnot(as.character(packageVersion('mistyR'))=='1.18.0')
capture.output(sessionInfo(),file=file.path(root,'installation_sessionInfo.txt'))
