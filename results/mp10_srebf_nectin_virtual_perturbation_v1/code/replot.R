s=Sys.getenv("JIN_PLOT_SKILL",unset="/home/data/t070721/.agents/skills/jins-routine-plotting")
for(n in c("theme","palettes","export"))source(file.path(s,"common",paste0(n,".R")))
jin_setup();setDTthreads(1);getd=function(id)fread(file.path("plot_data",paste0(id,".tsv")),na.strings=c("","NA"))
exports=list()
emit=function(id,p,w,h,caption,from=id){
 if(id!="F01_coverage")p=p+scale_x_continuous(n.breaks=3,labels=scales::label_number(scale=1000,accuracy=0.1))+labs(x=paste0(p$labels$x," (units of 0.001)"))
 jin_export(p+jin_theme()+theme(strip.background=element_rect(fill="#F1F2F3",colour=NA),panel.spacing=unit(1,"lines"),plot.margin=margin(8,16,8,8)),file.path("figures",id),w,h)
 exports[[length(exports)+1]]<<-data.table(figure=id,pdf=paste0("figures/",id,".pdf"),png=paste0("figures/",id,".png"),plot_data=paste0("plot_data/",from,".tsv"),caption=caption)
}
ink="#195571";orange="#AF5438";ko_cols=c(SREBF1="#195571",SREBF2="#AF5438")
id="F01_coverage";a=getd(id);a[,gene:=factor(gene,levels=rev(c("SREBF1","SREBF2","NECTIN2","NECTIN1","NECTIN3","NECTIN4","PVR","HMGCR","SQLE","LDLR")))]
p=ggplot(a,aes(x=100*detect_frac,y=gene))+geom_segment(aes(x=0,xend=100*detect_frac,yend=gene),colour="grey78")+geom_point(colour=ink,size=3)+geom_text(aes(label=sprintf("%.1f%%",100*detect_frac)),hjust=-.25,size=3.5)+scale_x_continuous(limits=c(0,100))+labs(x="Epithelial cells with detectable RNA (%)",y=NULL)
emit(id,p,6.2,4.6,"全部合格 Epi 的原始非零 UMI 比例。模型覆盖和入选原因见同名图数据与覆盖接口。")
id="F02_patients";a=getd(id);a[,patient_id:=factor(patient_id,levels=rev(sort(unique(patient_id))))]
p=ggplot(a,aes(effect_value,patient_id))+geom_point(data=a[is.finite(effect_value)],colour=ink,size=1.9)+geom_text(data=a[!is.finite(effect_value)],aes(x=0,label="NE"),hjust=0,size=3)+facet_wrap(~perturbed_gene,nrow=1)+labs(x="NECTIN2 regulatory distance",y="Recorded patient")
emit(id,p,9.2,10,"每点为该患者表达统计量在固定主 WT 核上的 NECTIN2 距离。两列各对应一次单基因扰动。患者差异为描述性结果。")
id="F03_states";a=getd(id);a[,state:=factor(state,levels=rev(paste0("C",0:9)))]
p=ggplot(a,aes(effect_value,state))+geom_point(data=a[is.finite(effect_value)],aes(colour=predicted_direction),size=3)+geom_text(data=a[!is.finite(effect_value)],aes(x=0,label="NE"),hjust=0,size=3)+facet_wrap(~perturbed_gene,nrow=1)+scale_colour_manual(values=c(down=ink,up=orange),name="Direction")+expand_limits(x=0)+labs(x="NECTIN2 regulatory distance",y="Existing epithelial cluster")
emit(id,p,8.6,4.9,"每簇独立重建 WT 并模拟两个单 KO。标 NE 的行缺乏建网覆盖；有方向颜色的点取模型原生方向评分。")
id="F04_controls";a=getd(id);a[,label:=factor(perturbed_gene,levels=rev(unique(c("SREBF1","SREBF2",sort(unique(a[role=="Control",perturbed_gene]))))))]
p=ggplot(a,aes(effect_value,label))+geom_segment(aes(x=0,xend=effect_value,yend=label),colour="grey82")+geom_point(aes(colour=role),size=3)+scale_colour_manual(values=c(SREBF=orange,Control=ink),name=NULL)+facet_wrap(~matched_to,scales="free_y")+labs(x="NECTIN2 regulatory distance",y="Perturbed gene")
emit(id,p,8.4,4.2,"每个 SREBF 匹配五个可扰动对照，允许两个匹配集合重叠。均来自同一 WT。匹配变量、距离与读出见图数据。")
id="F05_refits";a=getd(id);a[,label:=ifelse(analysis_type=="LODO",sub("drop_","",analysis_id),ifelse(analysis_type=="alternative_seed","Seed 456",sub("sub","",analysis_id)))];a[,label:=factor(label,levels=rev(c(sprintf("%02d",1:20),sort(unique(a[analysis_type=="LODO",label])),"Seed 456")))]
for(kind in c("direction_score","effect_value")){
 p=ggplot(a,aes(x=.data[[kind]],y=label))+geom_vline(xintercept=0,colour="grey75")+geom_point(data=a[is.finite(get(kind))],aes(colour=analysis_type),size=2.4)+geom_text(data=a[!is.finite(get(kind))],aes(x=0,label="NE"),size=3)+facet_wrap(~perturbed_gene,nrow=1)+scale_colour_manual(values=c(patient_subsample=ink,LODO=orange,alternative_seed="#725B81"),labels=c(alternative_seed="Seed",LODO="Dataset omitted",patient_subsample="Patient subset"),name=NULL)+labs(x=if(kind=="direction_score")"NECTIN2 direction score" else "NECTIN2 regulatory distance",y="Network refit")
 emit(if(kind=="direction_score")id else "F10_refit_distance",p,9.2,7.8,"01–20 为患者无放回重采样；数据集名称代表剔除该来源后重拟合；另列固定种子 456。展示每次估计。",id)
}
id="F06_secondary";a=getd(id)
p=ggplot(a,aes(effect_value,endpoint_gene))+geom_segment(data=a[is.finite(effect_value)],aes(x=0,xend=effect_value,yend=endpoint_gene),colour="grey82")+geom_point(data=a[is.finite(effect_value)],aes(colour=predicted_direction),size=3)+geom_text(data=a[!is.finite(effect_value)],aes(x=0,label="NE"),hjust=0,size=3)+facet_wrap(~perturbed_gene,nrow=1)+scale_colour_manual(values=c(down=ink,up=orange),name="Direction")+labs(x="Regulatory distance",y=NULL)
emit(id,p,8.8,4.4,"上皮网络中的 nectin 家族、PVR 及三个预设校准读出。模型内排序、P、BH q 与定向评分均见图数据。")
id="F07_programs";a=getd(id);a[,label:=gsub("_"," ",program)]
p=ggplot(a,aes(effect_value,label))+geom_point(data=a[is.finite(effect_value)],colour=ink,size=2.8)+geom_text(data=a[!is.finite(effect_value)],aes(x=0,label="NE"),hjust=0,size=3)+facet_grid(mask~perturbed_gene)+expand_limits(x=0)+labs(x="Median regulatory distance of covered program members",y=NULL)
emit(id,p,9,6.4,"Full 与 M20 分别沿用原成员与去 MP10 Top20 的成员。每行至少覆盖三基因及半数冻结成员。MP10 自身 M20 为空。")
id="F08_strata";a=getd(id);a[,label:=factor(state,levels=c("H05","H10","H50","L50"),labels=c("High (5%)","High (10%)","High (50%)","Low (50%)"))]
p=ggplot(a,aes(effect_value,label))+geom_segment(aes(x=uncertainty_low,xend=uncertainty_high,yend=label),colour="grey60")+geom_point(colour=ink,size=3)+facet_wrap(~perturbed_gene,nrow=1)+labs(x="Patient-conditioned NECTIN2 regulatory distance",y=NULL)
emit(id,p,8.5,3.5,"High 分别对应冻结 H05、H10、H50，Low 对应 L50。点为有资格患者的中位数，横线为患者实际最小到最大值。各层可能重叠。")
id="F09_direction";a=getd(id);a[,label:=factor(adjustment,levels=c("pc1_score","depth_score","direction_score"),labels=c("Depth and PC1 removed","Depth removed","Unadjusted"))]
p=ggplot(a,aes(score,label))+geom_vline(xintercept=0,colour="grey65")+geom_point(colour=ink,size=3)+facet_wrap(~perturbed_gene,nrow=1)+labs(x="NECTIN2 direction score",y=NULL)
emit(id,p,8.7,3.1,"方向评分来自 WT 表达。深度为所纳入基因的文库大小；最后一行敏感性去除深度残差表达矩阵的第一主轴。")
fwrite(rbindlist(exports),"FIGURE_INDEX.tsv",sep="\t")
writeLines(capture.output(sessionInfo()),"logs/plot_session.txt")

