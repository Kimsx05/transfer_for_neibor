d=rd('plot_data/project_profiles.tsv.gz')[contrast=='H05' & mask=='M20-masked' & module=='sterol']
d[,group:=factor(ifelse(profile_type=='L50','Low','High'),levels=c('Low','High'))]
fwrite(d[,.(sample_uid,dataset,group,score)],file.path(O,'plot_data/F06_sterol_pairs.tsv'),sep='\t')
p=ggplot(d,aes(group,score,group=sample_uid))+geom_line(colour='grey75',linewidth=.45)+geom_point(aes(colour=group),size=2.2,show.legend=FALSE)+scale_colour_manual(values=c(Low='#32739A',High='#B75230'))+expand_limits(y=0)+labs(x=NULL,y='Sterol program (mean log2 CPM)')
out(p,'F06_sterol_pairs',5.3,4.6)
d=rbindlist(lapply(c('CTRP','DepMap'),function(x)rd(paste0('tables/',x,'_associations.tsv'))[cohort=='BLCA']))
fwrite(d,file.path(O,'plot_data/F07_model_associations.tsv'),sep='\t')
p=ggplot(d,aes(effect,endpoint))+geom_vline(xintercept=0,color='grey60')+geom_segment(aes(x=lo,xend=hi,yend=endpoint))+geom_point(size=2.8,colour='#305B73')+scale_x_continuous(limits=c(-1,1))+labs(x='Spearman correlation with MP10',y=NULL)
out(p,'F07_model_associations',6.4,3.7)

