source(file.path(Sys.getenv('B_OUT'),'code/common.R'))
skill<-Sys.getenv('JIN_PLOT_SKILL','/home/data/t070721/.agents/skills/jins-routine-plotting')
for(f in c('theme','export'))source(file.path(skill,'common',paste0(f,'.R')))
jin_setup()
theme_set(jin_theme())
fig<-file.path(out,'figures');idx<-list()
draw<-function(p,d,name,w,h,cap){
 fwrite(d,file.path(fig,paste0(name,'.tsv')),sep='\t',na='NA')
 jin_export(p,file.path(fig,name),w,h)
 idx[[name]]<<-data.table(figure=name,PDF=paste0('figures/',name,'.pdf'),PNG=paste0('figures/',name,'.png'),data=paste0('figures/',name,'.tsv'),caption=cap)
}
labs0<-c(SCAP='SCAP dependent',FASN='FASN dependent',TCR='TCR, 8 h',CD36='CD36 dependent',SREBF2='SREBF2 dependent (VAT)',G2M='G2M checkpoint',CHOL_SYNTH='Sterol synthesis',LCFA_IMPORT='Long-chain FA import',CHOL_IMPORT='Cholesterol import')
keys<-names(labs0)[1:5]
pal<-c(CNP0000460='#0072B2',GSE326225='#CC79A7',HRA003620='#D55E00',InhouseData='#009E73',PRJNA662018='#665192')
cv<-readt(out,'coverage.tsv')
short<-function(x) sub('^(PRJNA662018_)+|^(CNP0000460_|GSE326225_|HRA003620_)','',sub('^.*::','',x))
cv[,label:=short(sample_uid)];cv[,class:=fcase(!main_specimen,'Repeat b',n_cells>=50,'Main',n_cells>=30,'30 only',default='Below 30')]
cv[,order_id:=factor(sample_uid,levels=cv[order(dataset,n_cells),sample_uid])]
p<-ggplot(cv,aes(order_id,n_cells,fill=class))+geom_col(width=.75)+geom_hline(yintercept=c(30,50),linetype='dotted',linewidth=.35)+facet_grid(~dataset,scales='free_x',space='free_x')+scale_x_discrete(labels=short)+scale_fill_manual(values=c('Main'='#427D93','Repeat b'='#AA7C46','30 only'='#86B8C4','Below 30'='#D2D9D9'),name=NULL)+labs(x=NULL,y='Treg cells per specimen')+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5))
# Apply custom layout after common theme.
draw(p,cv,'coverage',13,5,'全部40个标本；虚线为30和50细胞。a/b分别展示，b仅用于替换敏感性。')
m<-readt(out,'model_status.tsv')
d<-m[setting=='M0046'&mask=='Full'&adjustment=='total'&program%in%keys];d[,label:=factor(labs0[program],levels=rev(labs0[keys]))]
p<-ggplot(d,aes(beta,label))+geom_vline(xintercept=0,linewidth=.4,linetype=2,color='#777777')+geom_segment(aes(x=lo,xend=hi,yend=label),linewidth=.7,color='#335C67')+geom_point(size=2.7,color='#124E66')+labs(x='Treg program change per SD of Epi MP10 (95% CI)',y=NULL)
draw(p,d,'forest',7.4,4.1,'每行一个预登记全排序投影；原Q50_z、32份标本、数据集固定效应，横线为95%区间。')
en<-readt(out,'directional_enrichment.tsv')[model_id=='M0046'&mask=='Full'&set=='full_rank']
en[,label:=factor(labs0[program],levels=rev(labs0[keys]))]
p<-ggplot(en,aes(stat,label))+geom_vline(xintercept=0,linewidth=.4)+geom_segment(aes(x=0,xend=stat,yend=label),colour='#819EAA',linewidth=1)+geom_point(colour='#23566C',size=3)+labs(x='Correlation of oriented gene ranks',y=NULL)
draw(p,en,'rank_direction',6.7,4,'全基因方向一致性；正值指向WT依赖或8 h。共有基因数逐行提供，基因不构成独立患者重复。')
sc<-readt(out,'specimen_scores.tsv.gz')[mask=='Full'&main_patient==TRUE&n_cells>=50]
for(k in keys){
 d<-sc[program==k];f<-lm(score~dataset+Q50_z,data=d);d[,fitted_common_slope:=as.numeric(predict(f,newdata=d))]
 p<-ggplot(d,aes(Q50_z,score,color=dataset))+geom_line(aes(y=fitted_common_slope),linewidth=.55)+geom_point(size=2.5)+scale_color_manual(values=pal,name='Dataset')+labs(x='Epi MP10 median (reference SD)',y=paste0(labs0[k],' (reference SD)'))
 draw(p,d,paste0('scatter_',k),7.4,4.5,'每点为一份标本。颜色为数据集；线段使用主模型的共同斜率及各数据集截距，仅显示该数据集的暴露范围。')
}
d<-sc[,.(program,sample_uid,dataset,Q50_z,score)];d[,module:=fcase(program%in%c('SCAP','FASN','CD36','SREBF2'),'Lipid perturbations',program%in%c('TCR','G2M'),'Activation / cycle',default='Metabolic processes')];d[,module:=factor(module,levels=c('Lipid perturbations','Activation / cycle','Metabolic processes'))]
d[,label:=factor(labs0[program],levels=rev(unname(labs0)))];lvl<-unique(d[order(dataset,Q50_z),sample_uid]);d[,unit:=factor(sample_uid,levels=lvl)]
p<-ggplot(d,aes(unit,label,fill=score))+geom_tile()+facet_grid(module~dataset,scales='free',space='free')+scale_x_discrete(labels=short)+scale_fill_gradient2(low='#396897',mid='#F7F5F0',high='#AA4B38',midpoint=0,name='Reference SD')+labs(x='Specimen (increasing Epi MP10 within dataset)',y=NULL)+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=.5),strip.text.y=element_text(size=9),strip.text.x=element_text(size=9))
draw(p,d,'program_heatmap',13,6.7,'行按脂质扰动、活化／周期、代谢过程分组；每列为一个实际标本，各数据集内按Epi Q50升序。色标保留全部数值。')
d<-readt(out,'sender_receiver_plot_data.tsv');d[,panel:=ifelse(receiver=='CHOL_IMPORT','Epi efflux / Treg cholesterol import','Epi fatty acyl-CoA / Treg LCFA import')]
p<-ggplot(d,aes(x,y,colour=dataset))+geom_point(size=2.5)+facet_wrap(~panel,ncol=1,scales='free')+scale_color_manual(values=pal,name='Dataset')+labs(x='Epi High program (reference SD)',y='Treg import program (reference SD)')
draw(p,d,'supply_scatter',7.7,7,'同标本H05发送端与单份Treg摄取程序。每条路线18份。主HC3因singleton不可估，图中保留该标本。')
e<-data.table(step=c('Epi synthesis','Epi efflux','Carrier source','Treg import','Treg synthesis','Treg function','Transfer'),
 finding=c('High vs Low: sterol program +2.75 SD; q < 0.001','Efflux program +0.08 SD; q = 0.423. ABCA1 RNA decreased.','C1QC TAM: median 67% of captured APOE RNA.','8/12 cholesterol-import genes covered; main HC3 undefined.','SCAP and VAT SREBF2 dependence: both intervals span zero.','G2M up and OXPHOS down in the original receiver gene ranks.','No donor-labelled cargo or direct transfer measurement.'))
e[,y:=rev(seq_len(.N))]
p<-ggplot(e)+geom_text(aes(x=0,y=y,label=step),hjust=0,fontface='bold',size=3.7)+geom_text(aes(x=2.25,y=y,label=finding),hjust=0,size=3.7)+scale_x_continuous(limits=c(0,12.5),expand=expansion(mult=c(.01,0)))+scale_y_continuous(expand=expansion(add=.5))+labs(x=NULL,y=NULL)+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.line=element_blank())
draw(p,e,'supply_evidence',10.8,3.5,'供给各环节证据。完整中文问题、状态、数值与限制见tables/supply_evidence.tsv。')
fwrite(rbindlist(idx),file.path(out,'FIGURE_INDEX.tsv'),sep='\t')
writeLines(capture.output(sessionInfo()),file.path(out,'logs/figures_session.txt'))
