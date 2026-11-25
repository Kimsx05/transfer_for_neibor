library(Matrix)
library(data.table)
setDTthreads(1)
set.seed(123)
out <- Sys.getenv("OUT"); root <- normalizePath(file.path(out,"../.."))
sparse_input <- function(d){
  readint <- function(k) { f <- file.path(d,paste0(k,".bin")); readBin(f,integer(),n=file.size(f)/4,endian="little") }
  genes <- readLines(file.path(d,"genes.txt")); ids <- readLines(file.path(d,"ids.txt"))
  new("dgCMatrix",i=readint("i"),p=readint("p"),x=as.double(readint("x")),Dim=as.integer(c(length(genes),length(ids))),Dimnames=list(genes,ids))
}

