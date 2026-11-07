suppressPackageStartupMessages({library(data.table);library(Seurat);library(Matrix);library(SingleR);library(BiocParallel);library(jsonlite)})
args<-commandArgs(TRUE);root<-normalizePath(args[1]);old<-file.path(dirname(root),'hd_major_annotation_v1');set.seed(123);future::plan('sequential');options(mc.cores=4)
source(file.path(old,'code/load_result.R'))
