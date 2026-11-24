source("results/mp10_endothelial_hypothesis_validation_v1/code/common.R")
skill="/home/data/t070721/.agents/skills/jins-routine-plotting"
for(f in c("common/theme.R","common/palettes.R","common/export.R","plot_types/Stacked_composition/plot.R"))source(file.path(skill,f))
jin_setup();theme_set(jin_theme());pdfpath=function(n)file.path(O,"figures",n)
draw=function(n,d,p,w,h){put(d,paste0(n,".tsv"),"plot_data");jin_export(p,pdfpath(n),w,h)}
path=gettab("legacy_full_m20.tsv");path[,name:=sub("^HALLMARK_","",pathway)]
path=path[name%in%c("MYC_TARGETS_V1","MYC_TARGETS_V2","TNFA_SIGNALING_VIA_NFKB","UNFOLDED_PROTEIN_RESPONSE","MTORC1_SIGNALING","G2M_CHECKPOINT","E2F_TARGETS","OXIDATIVE_PHOSPHORYLATION","CHOLESTEROL_HOMEOSTASIS","ANGIOGENESIS","PI3K_AKT_MTOR_SIGNALING","INTERFERON_ALPHA_RESPONSE","INTERFERON_GAMMA_RESPONSE")]
pretty=function(x)gsub("_"," ",x,fixed=TRUE)
p=ggplot(path,aes(NES,reorder(pretty(name),NES,mean),color=mask))+geom_vline(xintercept=0,color="grey65")+geom_point(position=position_dodge(width=.4),size=2.6)+labs(x="Normalized enrichment score",y=NULL,color=NULL)+scale_color_manual(values=c(Full="#196384",M20="#bc6236"))
draw("01_full_m20",path,p,8,5.5)
d=gettab("coverage_model_status.tsv")[exposure=="Q50"];d[,label:=ifelse(status=="completed",as.character(n),paste0(n," (NE)"))]
p=ggplot(d,aes(type,n,fill=factor(min_cells)))+geom_col(position="dodge",width=.7)+geom_text(aes(label=label),position=position_dodge(.7),vjust=-.45,size=3.5)+scale_y_continuous(limits=c(0,21))+scale_fill_manual(values=c("30"="#8fafb9","50"="#196384"))+labs(x=NULL,y="Eligible patient units",fill="Minimum cells")+theme(axis.text.x=element_text(angle=25,hjust=1))
draw("02_coverage",d,p,7.2,4.2)
d=gettab("subtype_composition.tsv")[main_patient==TRUE & total>=50];level=unique(d[order(dataset,Q50_z),sample_uid]);d[,specimen:=sprintf("S%02d",match(sample_uid,level))]
colors=tab(B,"final_label_dictionary.tsv")[parent=="Endothelial",setNames(color,celltype)]
p=plot_stacked_composition(d,"specimen","type","fraction",condition_order=sprintf("S%02d",seq_along(level)),category_order=names(colors),colors=colors,ylab="Fraction of endothelial cells",xlab="Specimen")
draw("03_composition",d,p,8.5,4.6)
r=gettab("patient_effects.tsv.gz");d=r[context=="Capillary EC_50"&mask=="Full"&setting=="Q50"&family=="core"];d[,name:=pretty(program)]
forest=function(d,xlab)ggplot(d,aes(beta,reorder(name,beta)))+geom_vline(xintercept=0,color="grey55",linetype=2)+geom_segment(aes(x=lo,xend=hi,yend=reorder(name,beta)),color="#336981",linewidth=.7)+geom_point(size=2.8,color="#153b52")+labs(x=xlab,y=NULL)
draw("04_patient_effects",d,forest(d,"Rank score change per SD of Epi MP10 Q50 (95% CI)"),8.2,4.1)
d=gettab("perturbation_enrichment.tsv")[type=="Capillary"&mask=="Full"];d[,label:=sub(".*__","",contrast)];d[,label:=gsub("_KD_", " dependence at ", label, fixed=TRUE)];d[,label:=sub("^(RELA|PDHA1|SREBF2|KLF6)_interaction_", "\\1 dependence x TNF ", label)];d[,label:=gsub("minus", " minus ", label, fixed=TRUE)];d[,label:=paste0(sub("__.*","",contrast),"   ",pretty(label))];d[,star:=ifelse(q<.05,"*","")]
p=ggplot(d,aes(arm,factor(label,levels=rev(sort(unique(label)))),fill=aligned_NES))+geom_tile(color="white")+
 geom_text(aes(label=paste0(sprintf("%.2f",aligned_NES),star)),size=3.2)+scale_fill_gradient2(low="#c66c45",mid="#fbf8f0",high="#5f9caa",midpoint=0,name="Aligned NES")+labs(x="External reference arm",y=NULL)
draw("05_perturbation_alignment",d,p,6.9,8.3)
sc=gettab("project_scores.tsv.gz");d=sc[context=="Capillary EC_50"&mask=="Full"&main_patient==TRUE&program=="MYC_TARGETS_V1"]
p=ggplot(d,aes(Q50_z,value,color=dataset))+geom_point(size=3)+expand_limits(y=0)+labs(x="Epi MP10 Q50 (SD)",y="Capillary MYC V1 mean rank percentile",color="Dataset")
draw("06_paired_MYC",d,p,7.3,4.4)
d=r[context=="Capillary EC_50"&mask=="Full"&family=="core"&setting%in%c("CHOL_SYNTH","FA_ELONG","HALLMARK_MTORC1_SIGNALING")];d[,name:=pretty(program)]
p=ggplot(d,aes(beta,name))+facet_wrap(~setting,ncol=1)+geom_vline(xintercept=0,color="grey60",linetype=2)+geom_segment(aes(x=lo,xend=hi,yend=name),color="#336981")+geom_point(size=2)+labs(x="EC rank score change per SD of Epi M20 score (95% CI)",y=NULL)
draw("07_epi_m20_pairs",d,p,8.2,8)
d=r[context=="Capillary EC_50"&program%in%c("ATF4_growth","ATF4_ISR")&((mask=="Full"&setting%in%c("Q50","Q95","depth","cycle"))|(mask!="Full"&setting=="Q50"))];d[,name:=paste(mask,setting)]
p=ggplot(d,aes(beta,name))+facet_wrap(~program,ncol=1)+geom_vline(xintercept=0,color="grey55",linetype=2)+geom_segment(aes(x=lo,xend=hi,yend=name),color="#336981")+geom_point(size=2.5)+labs(x="Rank score change per SD of exposure (95% CI)",y=NULL)
draw("08_atf4_sensitivity",d,p,6.4,6.3)
d=gettab("secondary_prediction_effects.tsv")[type=="Capillary EC"&mask=="Full"&setting=="Q50"];d[,name:=pretty(program)]
draw("09_specific_outputs",d,forest(d,"Rank score change per SD of Epi MP10 Q50 (95% CI)"),7,3.1)

