O=commandArgs(TRUE)[1]
sk=Sys.getenv('JIN_SKILL',unset='/home/data/t070721/.agents/skills/jins-routine-plotting')
for(f in c('theme','palettes','export'))source(file.path(sk,'common',paste0(f,'.R')))
setDTthreads(1);jin_setup()
rd=function(n)fread(file.path(O,n));out=function(p,n,w=6,h=4.6)jin_export(p+jin_theme(),file.path(O,'figures',n),width=w,height=h)
x=rd('tables/project_statistics.tsv')[contrast=='H05' & model=='paired' & mask=='M20-masked']
g=x[!module%in%c('MAPK','PI3K','cycling','stress')]
p=ggplot(g,aes(effect,reorder(module,effect)))+geom_vline(xintercept=0,colour='grey65')+geom_segment(aes(x=lo,xend=hi,yend=reorder(module,effect)),linewidth=.55)+geom_point(size=2.6,color='#305B73')+labs(x='High − Low (mean log2 CPM difference)',y=NULL)
out(p,'F01_project_modules',7,4.8)
p=ggplot(x[module%in%c('MAPK','PI3K')],aes(effect,module))+geom_vline(xintercept=0,colour='grey65')+geom_segment(aes(x=lo,xend=hi,yend=module))+geom_point(size=2.7)+labs(x='High − Low (signed footprint difference)',y=NULL)
out(p,'F02_project_footprints',6,2.8)
for(platform in c('CTRP','DepMap')){
 d=rd(paste0('plot_data/',platform,'_units.tsv'))[eligible==TRUE]
 for(ep in unique(d$endpoint)){
  b=d[endpoint==ep];p=ggplot(b,aes(score,response))+geom_point(size=2.5,color='#305B73')+expand_limits(y=0)+labs(x='MP10 mean expression percentile',y=if(platform=='CTRP')paste(ep,'AUC')else 'GPX4 gene effect')
  out(p,paste0('F03_',gsub('[^A-Za-z0-9]+','_',ep)))
 }
}
d=rd('tables/GSE182638_modules.tsv')[mask=='Full' & module%in%c('MP10','GPX4_GSH','FSP1_CoQ','sterol','stress')]
p=ggplot(d,aes(effect,module))+geom_vline(xintercept=0,color='grey65')+geom_segment(aes(x=lo,xend=hi,yend=module))+geom_point()+facet_grid(line~contrast)+labs(x='Mean paired difference (log2 CPM)',y=NULL)
out(p,'F04_induction_rescue',8,5)
d=rd('plot_data/HMGCS1_source.tsv')[sheet%in%c('Fig. 3a','Fig. 3b') & condition=='SLO']
d[,group:=factor(group,c('shNT','shHMGCS1','shHMGCS1 + Cholesterol'),c('Control','HMGCS1 KD','KD + cholesterol'))]
p=ggplot(d,aes(group,value))+geom_point(position=position_jitter(width=.04,height=0,seed=23),size=2.6,color='#305B73')+stat_summary(fun=mean,geom='point',shape=95,size=9)+facet_wrap(~model)+expand_limits(y=0)+labs(x=NULL,y='Relative membrane CHMP4B')
out(p,'F05_CHMP4B_localization',8,4.4)

source(file.path(O,"scripts/paired_plots.R"))
