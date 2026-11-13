# Exactly the four metrics returned by nichenetr::predict_ligand_activities(single=TRUE).
# Avoids unused correlation tests, Wilcoxon-GST and iRegulon calculations from its internal helper.
nichenet_four_metrics<-function(geneset,background_expressed_genes,ligand_target_matrix,potential_ligands){
 bg<-background_expressed_genes;stopifnot(all(geneset%in%bg),all(bg%in%rownames(ligand_target_matrix)))
 y<-setNames(bg%in%geneset,bg);stopifnot(sd(y)>0);xx<-ligand_target_matrix[bg,potential_ligands,drop=FALSE];stopifnot(all(is.finite(xx)))
 native_split<-get('get_split_auroc',asNamespace('nichenetr'))
 rbindlist(lapply(seq_along(potential_ligands),function(i){x<-xx[,i];ro<-ROCR::prediction(x,y);pc<-ROCR::performance(ro,measure='prec',x.measure='rec');precision<-pc@y.values[[1]];recall<-pc@x.values[[1]];precision[is.na(precision)]<-1;recall[is.na(recall)]<-0;ap<-caTools::trapz(recall,precision);data.table(test_ligand=potential_ligands[i],auroc=native_split(x,y)$auroc,aupr=ap,aupr_corrected=ap-mean(y),pearson=cor(x,y))}))
}
