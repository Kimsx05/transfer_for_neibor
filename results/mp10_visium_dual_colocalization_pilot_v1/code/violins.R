source("code/init.R")
skill <- path.expand("~/.agents/skills/jins-routine-plotting")
for (f in c("common/theme.R","common/palettes.R","common/export.R","plot_types/Violin_distribution/plot.R")) source(file.path(skill,f))
jin_setup()
names <- c(Synthesis="Cholesterol synthesis",Homeostasis="Cholesterol homeostasis",Fatty_acid="Fatty acid metabolism",Efflux="Cholesterol efflux")
index <- list()
for(s in fread("tables/join_audit.tsv")$section_id){
 d <- fread(file.path("spots",paste0(s,".tsv.gz")))
 for(x in c("H05","MP10")) for(k in names(names)) for(v in c("full","minus20")){
  col <- paste(k,v,sep="__");reg <- paste0("region_",x)
  use <- d[query==TRUE & get(reg)%in%c("high","low")]
  if(!nrow(use))next
  p <- plot_violin_distribution(use,reg,col,method=names[k],ylab=paste(names[k],if(v=="minus20")"(minus Top20) score" else "score"),order=c("low","high"),colors=c(low="#A6CEE3",high="#D87869"))+ggplot2::labs(x=if(x=="H05")"Epi mapping score" else "MP10 score")
  name <- paste("region",x,k,v,sep="_");path<-file.path("figures",s,name)
  jin_export(p,path,width=5.2,height=4.5)
  index[[length(index)+1L]]<-data.table(section_id=s,kind="functional distribution",figure=name,pdf=paste0(path,".pdf"),png=paste0(path,".png"))
 }
}
fwrite(rbindlist(index),"figure_index_violin.tsv",sep="\t")

