source(file.path(Sys.getenv('B_OUT'),'code/common.R'))
mp<-readt(out,'sample_mapping.tsv')
hg<-readt(out,'gene_mapping.tsv.gz')
ortho<-hg[status=='one_to_one',.(mouse,human,MGI_class)]
stopifnot(!anyDuplicated(ortho$mouse),!anyDuplicated(ortho$human))
res<-list(); losses<-list()
for(k in c('SCAP','FASN','TCR','CD36','SREBF2')){
 d<-readt(out,paste0('expression_',k,'.tsv.gz'))
 g<-mp[contrast==k];cols<-g$GSM;stopifnot(setequal(cols,names(d)[-1]))
 rna<-k %in%c('CD36','SREBF2')
 if(rna)for(j in cols)set(d,j=j,value=log2(1+d[[j]]))
 dg<-d[,lapply(.SD,median),by=mouse_symbol,.SDcols=cols]
 mat<-as.matrix(dg[,..cols]);rownames(mat)<-dg$mouse_symbol
 ok<-apply(mat,1,function(x)all(is.finite(x))&&sd(x)>0)
 case<-as.numeric(g$arm=='case');fit<-lmFit(mat[ok,,drop=FALSE],model.matrix(~case))
 f<-eBayes(fit,trend=rna)
 z<-as.data.table(topTable(f,coef=2,number=Inf,sort.by='none',confint=TRUE),keep.rownames='mouse')
 z<-merge(z,ortho,by='mouse',all.x=TRUE);z[,contrast:=k]
 z[,direction:=fifelse(logFC>0,if(k=='TCR')'8h_UP' else 'KO_UP',if(k=='TCR')'8h_DOWN' else 'KO_DOWN')]
 z[,oriented_t:=if(k=='TCR')t else -t];setorder(z,-oriented_t)
 z[,rank_oriented:=seq_len(.N)];res[[k]]<-z
 loss<-data.table(mouse=rownames(mat),variable=ok)
 loss<-merge(loss,z[,.(mouse,human,MGI_class)],by='mouse',all.x=TRUE)
 loss[,mapping_status:=fcase(!variable,'constant_or_nonfinite',!is.na(human),'mapped',mouse%in%hg$mouse,'non_unique_or_no_human_in_MGI',default='absent_exact_mouse_symbol')]
 loss[,contrast:=k];losses[[k]]<-loss
}
z<-rbindlist(res);saveT(z,'perturbation_DE.tsv.gz');saveT(z[adj.P.Val<.05],'perturbation_strict_sets.tsv.gz')
saveT(rbindlist(losses),'DE_mapping_and_losses.tsv.gz')
saveT(z[,.(tested=.N,mapped=sum(!is.na(human)),strict_up=sum(adj.P.Val<.05&logFC>0),strict_down=sum(adj.P.Val<.05&logFC<0),contrast_unit=if(contrast=='TCR')'8h minus 0h' else 'KO minus WT'),by=contrast],'reference_summary.tsv')
