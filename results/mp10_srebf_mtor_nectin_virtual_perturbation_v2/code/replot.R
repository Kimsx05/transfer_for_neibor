suppressPackageStartupMessages({library(data.table);library(ggplot2);library(grid)})
s=Sys.getenv("JIN_PLOT_SKILL",unset="/home/data/t070721/.agents/skills/jins-routine-plotting")
for(n in c("theme","palettes","export"))source(file.path(s,"common",paste0(n,".R")))
jin_setup();setDTthreads(1);index=list();K=c("SREBF1","SREBF2","RPTOR","MTOR","RICTOR")
getd=function(id){a=fread(file.path("plot_data",paste0(id,".tsv")));a[,ko:=factor(ko,levels=K)];a}
emit=function(id,p,w,h,caption){
 p=p+jin_theme()+theme(panel.grid.major.y=element_line(colour="#efefef",linewidth=.25),strip.background=element_rect(fill="#f1f2f3",colour=NA),legend.position="bottom",plot.margin=margin(8,14,8,8))
 dest=file.path("figures",id);jin_export(p,dest,w,h);ggsave(paste0(dest,".svg"),p,device=grDevices::svg,width=w,height=h,bg="white")
 index[[length(index)+1]]<<-data.table(id=id,pdf=paste0(dest,".pdf"),svg=paste0(dest,".svg"),png=paste0(dest,".png"),plot_data=paste0("plot_data/",id,".tsv"),caption=caption)
}
id="F01";a=getd(id);a[,analysis:=factor(analysis,levels=c("v1","full","clean"),labels=c("v1","P_full","P_clean"))]
emit(id,ggplot(a,aes(analysis,count,fill=category))+geom_col(width=.7)+scale_fill_manual(values=c("Other"="#a6afb3","MT genome"="#8d4d6d","Mitochondrial ribosome"="#bd7840","Cytosolic ribosome"="#386e82","Unclassified"="#766e62"))+facet_wrap(~ko,nrow=1)+labs(x=NULL,y="Genes among top 50",fill=NULL,subtitle="8,000 cells · 38 patients | Historical v1 has two perturbations"),10,4,"主种子前50名的精确注释类别。无历史模型的列留空。")
id="F02";a=getd(id);a[,gene:=reorder(gene,-rank,median)]
emit(id,ggplot(a,aes(ko,gene,fill=rank))+geom_tile(colour="white")+geom_text(aes(label=round(rank)),size=3.5)+scale_fill_gradient(low="#b3d9d8",high="white")+labs(x="Perturbed source",y=NULL,fill="Rank",subtitle="Union of five top-eight lists | P_clean · 8,000 cells"),7.2,5,"各源前八名的并集，在同一主模型内列出非源基因名次。")
id="F03";a=getd(id)
emit(id,ggplot(a,aes(rank_fraction,gene,colour=direction))+geom_point(size=2.3)+facet_wrap(~ko,nrow=1)+scale_x_continuous(limits=c(0,1),breaks=c(0,.5,1))+scale_colour_manual(values=c(down="#195571",up="#b05c3e"))+labs(x="Distance rank fraction (0 = highest)",y=NULL,colour="Predicted direction",subtitle="P_clean · 2,193 non-source genes per perturbation"),10.2,5.8,"预设脂质、NECTIN 和 mTOR 相关读出，位次与方向分列。")
id="F04";a=getd(id);a[,readout:=gsub("HALLMARK_MTORC1_SIGNALING","Hallmark mTORC1",readout)]
emit(id,ggplot(a,aes(median_rank_fraction,readout))+geom_segment(aes(xend=control_rank_fraction,yend=readout),colour="#999999")+geom_point(aes(colour="Source"),size=2.4)+geom_point(aes(x=control_rank_fraction,colour="Matched median"),shape=1,size=2.5)+facet_wrap(~ko,nrow=1)+scale_colour_manual(values=c("Source"="#195571","Matched median"="#aa5037"))+labs(x="Median rank fraction in a program (0 = highest)",y=NULL,colour=NULL,subtitle="Full programs · common non-source members · 38 patients"),12,5.2,"五个源及各自五对照在共同程序成员上的中位排名，完整定量见图数据。")
id="F05";a=getd(id);a[,label:=factor(refit_type,levels=c("clean","seed","patient_subsample","LODO","state"),labels=c("Main","Other seeds","Patient subsets","Source omitted","Clusters"))]
emit(id,ggplot(a[is.finite(rank_fraction)],aes(100*rank_fraction,label))+geom_point(position=position_jitter(width=0,height=.12,seed=123),size=1.7,alpha=.75,colour="#195571")+facet_wrap(~ko,nrow=1)+labs(x="NECTIN2 distance rank percentile (0 = highest)",y=NULL,subtitle="1 main · 2 alternative seeds · 20 patient subsets · 5 source omissions · 9 eligible clusters"),11,4.3,"全部真实重拟合逐点显示。C9 有29个细胞、2位患者，未达到簇建网门槛。")
id="F06";a=getd(id)
emit(id,ggplot(a,aes(ko,gene,fill=n_controls_lower))+geom_tile(colour="white")+geom_text(aes(label=n_controls_lower),size=4)+scale_fill_gradient(limits=c(0,5),low="white",high="#83b7bd")+labs(x="Perturbed source",y=NULL,fill="Controls below",subtitle="Number of five matched controls with smaller native distance | 38 patients"),8,4.8,"与五个固定对照逐基因比较，仅作描述，不换算经验 P 值。")
id="F07";a=getd(id);a[,symbol:=fcase(comparison=="NE","NE",is.na(signed_delta),"NE",signed_delta>0,"Up",signed_delta<0,"Down",default="Zero")]
emit(id,ggplot(a,aes(ko,gene,fill=comparison))+geom_tile(colour="white")+geom_text(aes(label=symbol),size=3.6)+scale_fill_manual(values=c(NE="#e0e2e3",agree="#9ac6c4",disagree="#e0b09a",zero="#ccccaa"))+labs(x="Perturbed TF",y=NULL,fill="Direction comparison",subtitle="Text = CellOracle signed delta direction | native 3 steps · patient-equal medians"),7.2,6,"同一队列两种模型的预测方向，颜色指示一致性，无法评估保留 NE。")
id="F08";a=getd(id);a[,stratum:=factor(stratum,levels=c("H05","H10","H50"),labels=c("High 5%","High 10%","High 50%"))]
emit(id,ggplot(a,aes(median,stratum))+geom_vline(xintercept=0,colour="grey65")+geom_segment(aes(x=min,xend=max,yend=stratum),colour="#9fadae")+geom_point(colour="#195571",size=2.5)+geom_text(aes(label=paste0("n=",count)),vjust=-.75,size=2.7)+facet_wrap(~ko,nrow=1)+labs(x="Paired NECTIN2 distance difference (High minus Low 50%)",y=NULL,subtitle="Frozen MP10 labels · common patients per pair · median and observed min–max"),12,3.8,"每一对 High/Low 的共同合格患者，线为实际患者范围，n 为患者数。")
fwrite(rbindlist(index),"FIGURE_INDEX.tsv",sep="\t");writeLines(capture.output(sessionInfo()),"logs/plot_session.txt")

