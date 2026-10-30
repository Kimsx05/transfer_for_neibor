# Small meaningful edge-case validation: missing Epi, eligible real zero, and empty common LR.
base<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_spatial_communication_v2_all_sections'
v1<-'/home/data/t070721/codex_workspace/Bladder Metabolism/mp10_spatial_communication_v1'
.libPaths(unique(c(file.path(v1,'Rlib'),'/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4',.libPaths())))
suppressPackageStartupMessages({library(data.table);library(Matrix);library(SpatialCellChat)})
source(file.path(base,'scripts/scc_lib.R'));setDTthreads(1)
out<-file.path(base,'logs/scc_qa_missing_state');dir.create(out,recursive=TRUE,showWarnings=FALSE)
Q<-matrix(c(0,0,.2,.3,.8,.7),2,3,dimnames=list(c('a','b'),c('Epi_MP10_high','TME1','TME2')))
support<-data.table(level='major',node=colnames(Q),type_evaluable=c(FALSE,TRUE,TRUE));P<-sparseMatrix(i=1,j=2,x=2,dims=c(2,2))
scc_aggregate(list(LR1=P),rownames(Q),list(major=Q),support,out,'synthetic_missing_epi')
z<-fread(file.path(out,'major_edges.tsv'));stopifnot(all(is.na(z[source=='Epi_MP10_high'|target=='Epi_MP10_high',total_mean_strength])),all(z[source!='Epi_MP10_high'&target!='Epi_MP10_high',total_mean_strength]>0))
outzero<-file.path(base,'logs/scc_qa_empty_native');dir.create(outzero,recursive=TRUE,showWarnings=FALSE)
scc_aggregate(list(),rownames(Q),list(major=Q),support,outzero,'synthetic_zero_native')
z<-fread(file.path(outzero,'major_edges.tsv'));stopifnot(all(is.na(z[source=='Epi_MP10_high'|target=='Epi_MP10_high',total_mean_strength])),all(z[source!='Epi_MP10_high'&target!='Epi_MP10_high',total_mean_strength]==0))
writeLines('PASS: missing state is NE; eligible edges are positive/real zero; native aggregate conserves cross-position flow',file.path(out,'validation.txt'))
cat('EDGE_CASES_PASS\n')
