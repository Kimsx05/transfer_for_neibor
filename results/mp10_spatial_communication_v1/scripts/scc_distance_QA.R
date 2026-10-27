REQUESTED_R_LIBS<-c('/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4')
.libPaths(unique(c('mp10_spatial_communication_v1/Rlib',REQUESTED_R_LIBS,.libPaths())))
suppressPackageStartupMessages({library(data.table);library(SpatialCellChat)})
for(sec in c('GSE171351_A1','Inhousedata_T4469479')){
 x<-fread(file.path('mp10_spatial_communication_v1/input_review/sections',sec,'coords.tsv'));xy<-as.matrix(x[,.(x_um,y_um)]); d<-as.matrix(dist(xy));actual<-d<=200+1e-6 & d>0
 r<-computeCellDistance(xy,interaction.range=200,contact.range=105,ratio=1,tol=1e-6)
 z<-data.table(section=sec,expected_crosslinks=sum(actual),native_crosslinks=length(r$d.spatial@x),boundary_200_links=sum(abs(d-200)<1e-6),max_abs_error=max(abs(as.matrix(r$d.spatial)-d*actual)),status='PASS')
 stopifnot(z$expected_crosslinks==z$native_crosslinks,z$max_abs_error<1e-6)
 fwrite(z,file.path('mp10_spatial_communication_v1/spatialcellchat',paste0(sec,'_distance_QA.tsv')),sep='\t')
}
