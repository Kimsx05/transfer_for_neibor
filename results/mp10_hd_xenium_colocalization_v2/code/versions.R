library(data.table);setDTthreads(1);packages=c("Seurat","SeuratObject","spatstat.geom","spatstat.explore","spatstat.random","sf","UCell","spacexr","Rcpp")
fwrite(data.frame(package=packages,version=sapply(packages,function(k)if(requireNamespace(k,quietly=TRUE))as.character(packageVersion(k)) else NA_character_)),file.path("audit","R_packages.tsv"),sep="\t")
fwrite(data.frame(package="UCell",version=as.character(packageVersion("UCell",lib.loc=file.path(path.expand("~"),"R/x86_64-pc-linux-gnu-library/4.5"))),usage="v1 cached score provenance and incremental scoring library"),"audit/scoring_R_package.tsv",sep="\t")
writeLines(capture.output(sessionInfo()),file.path("audit","R_session.txt"))
