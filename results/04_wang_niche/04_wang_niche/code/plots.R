#!/usr/bin/env Rscript
Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='')
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(jsonlite)})
setDTthreads(1L)
folder<-commandArgs(trailingOnly=TRUE)[[1]];root<-dirname(dirname(folder))
params<-fromJSON(file.path(root,'parameters.json'));ds<-basename(folder)
# Reuse only the existing save_both function AST, never source the executing QC script.
old_code<-'/home/data/t070721/codex_workspace/Bladder Metabolism/ST/code/03_spatial_QC_summary.R'
old_expr<-parse(old_code)
selected<-Filter(function(z)is.call(z)&&identical(z[[1]],as.name('<-'))&&identical(z[[2]],as.name('save_both')),as.list(old_expr))
stopifnot(length(selected)==1L);eval(selected[[1]])
replace_dpi<-function(z){if(is.call(z)||is.pairlist(z)){for(j in seq_along(z))z[[j]]<-replace_dpi(z[[j]])};if(is.numeric(z)&&length(z)==1L&&z==180)z<-600;z}
body(save_both)<-replace_dpi(body(save_both))
capture.output(save_both,file=file.path(folder,'adapted_R_save_both.txt'))
out<-file.path(folder,'figures');dir.create(out,showWarnings=FALSE)
thumb<-file.path(folder,'thumbnails');dir.create(thumb,showWarnings=FALSE)
index<-list();i<-0L
save_view<-function(p,stem,width=7,height=5,view='',section=''){
  p<-p+theme(plot.background=element_rect(fill='white',colour=NA))
  if(!is.null(p$labels$caption))p$labels$caption<-paste(strwrap(p$labels$caption,width=95),collapse='\n')
  save_both(p,out,stem,width,height)
  ggsave(file.path(thumb,paste0(stem,'.png')),p,width=width,height=height,dpi=110,limitsize=FALSE)
  i<<-i+1L;index[[i]]<<-data.table(dataset=ds,section_id=section,view=view,stem=stem,
    pdf=normalizePath(file.path(out,paste0(stem,'.pdf'))),png=normalizePath(file.path(out,paste0(stem,'.png'))),
    thumbnail=normalizePath(file.path(thumb,paste0(stem,'.png'))),png_dpi=600L,status='PASS')
  cat('FIGURE',ds,stem,'\n')
}
c<-fread(file.path(folder,'niche_centroids.tsv'));niches<-sort(unique(c$niche_id));nf<-length(niches)
palette<-setNames(hcl.colors(nf,palette='Dark 3'),niches)
fwrite(data.table(dataset=ds,niche_id=niches,colour=unname(palette)),file.path(folder,'niche_palette.tsv'),sep='\t')
c[,factor:=factor(factor,levels=rev(params$factor_order))];c[,niche_id:=factor(niche_id,levels=niches)]
p<-ggplot(c,aes(niche_id,factor,fill=mean_fraction))+geom_tile()+
  scale_fill_viridis_c(name='Mean fraction',limits=c(0,1),na.value='#bdbdbd')+
  geom_tile(data=c[mean_fraction==0],fill='black')+
  theme_bw(base_size=9)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=60,hjust=1))+
  labs(x=NULL,y=NULL,caption='All 30 original subtype factors; zero=black, NA=grey. Caution flags in factor_interpretation.tsv.')
save_view(p,'niche_subtype_composition',max(6,0.32*nf+3.5),9,view='subtype_mean_fraction_heatmap')
m<-fread(file.path(folder,'niche_major_summary.tsv'));majors<-sort(unique(m$major))
major_palette<-setNames(hcl.colors(length(majors),'Set 2'),majors)
p<-ggplot(m,aes(niche_id,mean_fraction,fill=major))+geom_col()+
  scale_fill_manual(values=major_palette)+scale_y_continuous(labels=scales::percent_format(),limits=c(0,1.0000001))+
  theme_bw(base_size=9)+theme(axis.text.x=element_text(angle=60,hjust=1))+
  labs(x=NULL,y='Mean composition fraction',caption='Major sums explain selected subtype clusters; majors were not clustering features.')
save_view(p,'niche_major_composition',max(7,0.35*nf+3),5.5,view='major_interpretation_stacked_bar')
g<-fread(file.path(folder,'silhouette_grid.tsv'));selected_parameters<-fread(file.path(folder,'selected_parameters.tsv'))
p<-ggplot(g[eligible==TRUE],aes(resolution,mean_silhouette,color=factor(k),group=k))+
  geom_line()+geom_point(size=1.5)+geom_point(data=selected_parameters,shape=21,size=4,fill='white')+
  theme_bw(base_size=10)+labs(x='Louvain resolution',y='Evaluation-spot mean silhouette',color='Composition k',
    caption=paste(unique(g$evaluation_mode),'| Shared evaluation set for all candidates; circle=selected.'))
save_view(p,'silhouette_parameter_grid',7,5,view='silhouette_grid')
coverage<-fread(file.path(folder,'niche_sample_coverage.tsv'))[coverage_scope=='section']
p<-ggplot(coverage,aes(niche_id,section_id,fill=fraction_valid_spots))+geom_tile()+
  scale_fill_viridis_c(name='Fraction of\nvalid spots',limits=c(0,1),na.value='#bdbdbd')+
  geom_tile(data=coverage[n_spots==0],fill='black')+
  theme_bw(base_size=9)+theme(panel.grid=element_blank(),axis.text.x=element_text(angle=60,hjust=1))+
  labs(x=NULL,y=NULL,caption='Section coverage; zero=black, NA=grey. Section IDs are not patient IDs.')
save_view(p,'niche_section_coverage',max(7,0.35*nf+3),max(4,0.3*uniqueN(coverage$section_id)+2),view='section_coverage')
meta<-fread(file.path(folder,'spatial_plot_metadata.tsv'));spots<-fread(file.path(folder,'spatial_plot_spots.tsv'))
geometry<-list()
for(sec in meta$section_id){
  info<-meta[section_id==sec];d<-spots[section_id==sec];im<-png::readPNG(info$image_path)
  # alpha on H&E only; valid categorical spots have fixed opacity 1.
  if(length(dim(im))==3L&&dim(im)[3]==3L)im<-array(c(im,rep(0.25,dim(im)[1]*dim(im)[2])),dim=c(dim(im)[1],dim(im)[2],4L))
  else if(length(dim(im))==3L&&dim(im)[3]==4L)im[,,4]<-0.25
  xy<-as.matrix(d[,.(x,y)]);nn<-BiocNeighbors::findKNN(xy,k=1L,BNPARAM=BiocNeighbors::KmknnParam())$distance[,1]
  spacing<-median(nn);diam<-0.75*spacing
  stopifnot(is.finite(spacing),spacing>0,all(is.finite(xy)),min(d$x)>=0,min(d$y)>=0,
      max(d$x)<=info$image_width,max(d$y)<=info$image_height)
  d[,label:=ifelse(is.na(niche_id),'NA composition',niche_id)]
  palette_full<-c(palette,'NA composition'='#bdbdbd');d[,label:=factor(label,levels=names(palette_full))]
  # Coordinate sign change is display only: image x right / y down, no rescaling.
  p<-ggplot(d,aes(x,-y,color=label))+
    annotation_raster(im,xmin=0,xmax=info$image_width,ymin=-info$image_height,ymax=0)+
    geom_point(size=diam*3.8*25.4/info$image_width,shape=16,alpha=1,stroke=0)+
    scale_color_manual(values=palette_full,drop=FALSE,name='Dataset niche')+
    coord_fixed(xlim=c(0,info$image_width),ylim=c(-info$image_height,0),expand=FALSE)+
    theme_void(base_size=8)+theme(legend.key.height=grid::unit(3,'mm'))+
    labs(caption='H&E alpha=0.25; spots alpha=1. Dataset-scoped labels; NA composition=grey.')
  save_view(p,paste0(sec,'__niche_labels'),6.6,max(5.2,0.20*(nf+1)+1.5),view='spatial_niche_labels',section=sec)
  geometry[[length(geometry)+1L]]<-data.table(dataset=ds,section_id=sec,n_spots=nrow(d),
      n_assigned=sum(!is.na(d$niche_id)),n_NA=sum(is.na(d$niche_id)),image_width=info$image_width,
      image_height=info$image_height,median_spacing_pixels=spacing,spot_diameter_pixels=diam,
      coordinate_orientation='x right, y down; lowres pixels; no extra scaling',status='PASS')
}
fwrite(rbindlist(geometry),file.path(folder,'plot_geometry_validation.tsv'),sep='\t')
fwrite(rbindlist(index),file.path(folder,'figure_index.tsv'),sep='\t')
write_json(list(dataset=ds,status='PASS',n_views=i,png_dpi=600L,HE_alpha=0.25,spot_alpha=1,
    plotting_source=old_code,adapter='Reuse save_both function AST with PNG dpi 180->600; spatial ggplot theme adapted from ST/code/16_run_full_RCTD_one.R',
    NA_colour='#bdbdbd'),file.path(folder,'plot_status.json'),pretty=TRUE,auto_unbox=TRUE)
