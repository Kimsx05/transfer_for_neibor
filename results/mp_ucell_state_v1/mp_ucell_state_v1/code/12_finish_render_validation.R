s <- paste(readLines('/home/data/t070721/codex_workspace/Bladder Metabolism/mp_ucell_state_v1/code/04_figures.R'),collapse='\n')
prefix <- strsplit(s,'for(i in seq_len(nrow(inst))) {',fixed=TRUE)[[1]][1]
eval(parse(text=prefix))
# The prior palette loop finished all 8 branches; verify mappings before final completion.
for(i in seq_len(nrow(inst))) {
 cl<-fread(out(inst$directory[i],'clusters_res0.4.tsv.gz'));lev<-paste0('C',sort(as.integer(sub('C','',unique(cl$cluster)))))
 expected<-setNames(grDevices::hcl(h=(seq(0,360,length.out=length(lev)+1)[seq_along(lev)]+37*i)%%360,c=60,l=65),lev)
 pal<-fread(out(inst$directory[i],'cluster_colors.tsv'));stopifnot(identical(pal$cluster,names(expected)),identical(pal$color,unname(expected)))
}
# Fresh final overview renders after every worker has exited.
dir<-out('04_full_raw');emb<-qread(file.path(dir,'umap.qs'));d<-as.data.table(emb,keep.rownames='cell_id');xy<-d
saveforce<-function(p,base,w,h){ggsave(paste0(base,'.pdf'),p,width=w,height=h,bg='transparent',useDingbats=FALSE);ggsave(paste0(base,'.png'),p,width=w,height=h,dpi=300,bg='transparent')}
saveforce(continuous(d,'EpiMP10',xy),file.path(dir,'figures/FULL_RAW__ALL__res0.4__MP10'),8,6)
saveforce(wrap_plots(lapply(mp,function(m)continuous(d,m,xy)),ncol=4),file.path(dir,'figures/FULL_RAW__ALL__res0.4__all14_raw_UCell'),20,16)
writeLines('complete',out('logs/PALETTE_FINALIZATION_COMPLETE'))
files<-list.files(root,pattern='\\.pdf$',recursive=TRUE,full.names=TRUE)
qc<-rbindlist(lapply(files,function(f) {info<-tryCatch(pdftools::pdf_info(f),error=function(e)e);data.table(path=f,pages=if(inherits(info,'error'))NA_integer_ else info$pages,pass=!inherits(info,'error')&&info$pages>0,error=if(inherits(info,'error'))conditionMessage(info)else'')}))
write_tsv(qc,out('00_input_audit/PDF_parser_validation.tsv'));stopifnot(all(qc$pass));writeLines('complete',out('logs/PDF_VALIDATION_COMPLETE'))
