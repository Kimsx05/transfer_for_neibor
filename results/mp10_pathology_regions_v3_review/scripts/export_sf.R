root=normalizePath('.');out=file.path(root,'results/mp10_pathology_regions_v3');.libPaths(unique(c(file.path(root,'mp10_pathology_regions_v2/Rlib'),'/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.5','/home/data/t070721/R/x86_64-pc-linux-gnu-library/4.4','/refdir/Rlib',.libPaths())))
suppressPackageStartupMessages(library(sf));sf_use_s2(FALSE)
files=list.files(file.path(out,'by_dataset'),pattern='^lesion_rois_um.geojson$',recursive=TRUE,full.names=TRUE);records=list()
for(f in files){
 j=jsonlite::read_json(f,simplifyVector=FALSE);if(!length(j$features))next
 roi=st_read(f,quiet=TRUE);st_crs(roi)=NA;stopifnot(all(st_is_valid(roi)));roi$area_um2=as.numeric(st_area(roi));attr(roi,'coordinate_units')='micrometers';attr(roi,'coordinate_system')='local Visium array; no geodetic CRS';attr(roi,'boundary_definition')='union of retained spot Voronoi hexagons; valid outer segments stored independently'
 saveRDS(roi,file.path(dirname(f),'lesion_rois_sf.rds'),compress='xz');records[[length(records)+1]]=data.frame(section=basename(dirname(f)),n_roi=nrow(roi),coordinate_units='um',area_um2=sum(roi$area_um2),status='valid sf polygon adapter; SPATA2 not run')
}
write.table(do.call(rbind,records),file.path(out,'tables/sf_geometry_export.tsv'),sep='\t',quote=FALSE,row.names=FALSE);writeLines(capture.output(sessionInfo()),file.path(out,'logs/sf_sessionInfo.txt'))
