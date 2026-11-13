select_model<-function(p,pr){
 ix<-if(p$patient=='main')pr$main_patient else pr$replacement_patient;sel<-pr[ix & n_cells>=p$min_cells]
 if(p$omit_dataset!='none')sel<-sel[dataset!=p$omit_dataset];if(p$only_dataset!='ALL')sel<-sel[dataset==p$only_dataset]
 if(p$kind=='Epi'){
 if(p$subset=='common'){ok<-sel[profile_type %in% c('H50','H10','H05','L50'),.N,by=sample_uid][N==4,sample_uid];sel<-sel[sample_uid %in% ok]}
 ok<-sel[profile_type %in% c(p$target,'L50'),.N,by=sample_uid][N==2,sample_uid];sel<-sel[sample_uid %in% ok & profile_type %in% c(p$target,'L50')];setorder(sel,sample_uid,profile_type);sel[,specimen:=factor(sample_uid)];sel[,state:=factor(ifelse(profile_type=='L50','Low','High'),levels=c('Low','High'))]
 }else {sel<-sel[profile_type==p$target];setorder(sel,sample_uid);sel[,exposure_z:=get(paste0(p$exposure,'_z'))];sel[,dataset_factor:=factor(dataset)]}
 if(uniqueN(sel$patient_uid_recorded)<4)return(NULL)
 des<-if(p$kind=='Epi')model.matrix(~specimen+state,sel)else if(uniqueN(sel$dataset)>1)model.matrix(~dataset_factor+exposure_z,sel)else model.matrix(~exposure_z,sel)
 if(qr(des)$rank<ncol(des)||nrow(des)-ncol(des)<3)return(NULL)
 list(profiles=sel,design=des)
}
